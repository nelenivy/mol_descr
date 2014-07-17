#include "molecule_manager.h"

#include <algorithm>
#include <iostream>

#include "extensions.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"
#include "types_calculation.h"

namespace molecule_descriptor
{

void MoleculeManager::SetCurrFilePrefix(const std::string& file_prefix)
{
	Clear();
	m_curr_file_prefix = file_prefix;
}

void MoleculeManager::FindSingularPoints(const bool calculate, const bool calc_as_average)
{
	if (calculate)
	{
		CalculateSingularPoints(calc_as_average);
		WriteSingularPoints();
		WriteSegmentedSurface();
		WriteSurfaceWithTypes();
	}
	else
	{
		const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
		ReadVector(non_marked_points_file_name, m_non_marked_singular_points.first);
		ReadSingularPointsLevelsTypes(m_levels_num);
	}
}

void MoleculeManager::ReadSingularPointsLevelsTypes(const int levels_num)
{
	//find out how many levels are there
	/*std::ifstream test_stream;
	int levels_num = -1; 
	do
	{	
		test_stream.close();
		++levels_num;
		const std::string non_marked_points_file_name_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypesLevels() + std::to_string(levels_num);
		test_stream.open(non_marked_points_file_name_lev);
	}
	while(test_stream.is_open());*/
	CV_Assert(levels_num > 0);
	m_non_marked_singular_points_levels.first.resize(levels_num);
	//read levels
	for (int level = 0; level < levels_num; ++level)
	{
		const std::string non_marked_points_file_name_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypesLevels() + std::to_string(level);
		ReadVector(non_marked_points_file_name_lev, m_non_marked_singular_points_levels.first[level]);			
	}
}
void MoleculeManager::ReadSingularPointsLevelsLabels(const int levels_num)
{
	////find out how many levels are there
	//std::ifstream test_stream;
	//int levels_num = -1; 
	//do
	//{	
	//	test_stream.close();
	//	++levels_num;
	//	const std::string non_marked_points_file_name_l_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabelsLevels() + std::to_string(levels_num);
	//	test_stream.open(non_marked_points_file_name_l_lev);
	//}
	//while(test_stream.is_open());
	CV_Assert(levels_num > 0);
	m_non_marked_singular_points_levels.second.resize(levels_num);
	//read levels
	for (int level = 0; level < levels_num; ++level)
	{
		const std::string non_marked_points_file_name_l_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabelsLevels() + std::to_string(level);
		ReadVector(non_marked_points_file_name_l_lev,m_non_marked_singular_points_levels.second[level]);
	}
}
void MoleculeManager::ClassifySingularPoints()
{
	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
	ReadVector(non_marked_points_file_name, m_non_marked_singular_points.first);
	ReadSingularPointsLevelsTypes(m_levels_num);

	CalculatePropertiesTypes();
	WriteSingularPoints();
}

void MoleculeManager::FindPairs(const bool calculate)
{
	if (calculate)
	{
		CalculatePairs();
		WritePairs();
	}
	else
	{
		ReadPairs();
	}
}
void MoleculeManager::FindPairsLevels(const bool calculate)
{
	if (calculate)
	{
		CalculatePairsLevels();
		WritePairsLevels();
	}
	else
	{
		ReadPairsLevels(m_levels_num);
	}
}

void MoleculeManager::FindTriples(const bool calculate)
{
	if (calculate)
	{
		CalculateTriples();
		WriteTriples();
	}
	else
	{
		ReadTriples();
	}
}

void MoleculeManager::AppendDistances(std::vector<double>& distances)
{
	distances.reserve(distances.size() + m_pairs.size());
	
	for (auto iter = m_pairs.begin(); iter != m_pairs.end(); ++iter)
	{
		distances.push_back(std::get<2>(*iter));//push distance
	}
}

void MoleculeManager::AppendDistancesLevels(std::vector<std::vector<double>>& distances_levels)
{
	CV_Assert(distances_levels.size() == 0 || distances_levels.size() == m_pairs_levels.size());
	distances_levels.resize(m_pairs_levels.size());

	for (size_t level = 0; level < distances_levels.size(); ++level)
	{
		for (auto iter = m_pairs_levels[level].begin(); iter != m_pairs_levels[level].end(); ++iter)
		{
			distances_levels[level].push_back(std::get<2>(*iter));//push distance
		}
	}
}

void MoleculeManager::AppendCharges(std::vector<double>& charges, const bool levels)
{
	const auto& vec = levels ? m_non_marked_singular_points_levels.first[0] : m_non_marked_singular_points.first;
	charges.reserve(charges.size() + vec.size());

	for (auto iter = vec.begin(); iter != vec.end(); ++iter)
	{
		charges.push_back(iter->Property().ElectricPotential());
	}
}
void MoleculeManager::AppendLennardJones(std::vector<double>& lennard_jones, const bool levels)
{
	const auto& vec = levels ? m_non_marked_singular_points_levels.first[0] : m_non_marked_singular_points.first;
	lennard_jones.reserve(lennard_jones.size() + vec.size());

	for (auto iter = vec.begin(); iter != vec.end(); ++iter)
	{
		lennard_jones.push_back(iter->Property().LennardJones());
	}
}
void  MoleculeManager::AppendAreas(std::vector<double>& areas, const bool levels)
{
	const auto& vec = levels ? m_non_marked_singular_points_levels.first[0] : m_non_marked_singular_points.first;
	areas.reserve(areas.size() + vec.size());

	for (auto iter = vec.begin(); iter != vec.end(); ++iter)
	{
		areas.push_back(iter->Property().Area());
	}
}

void MoleculeManager::CalculatePairsWithTypes()
{
	CV_Assert(!m_dist_threshes.empty());
	m_pairs_with_types.clear();
	m_pairs_with_types.reserve(m_pairs.size());
	m_types.clear();
	m_types.reserve(m_pairs.size());
	const size_t points_types_num = GetPointsTypeNum();
	const size_t pairs_max_type = GetPairsTypeNum();
	m_pairs_histogram.resize(pairs_max_type);
	std::fill(m_pairs_histogram.begin(), m_pairs_histogram.end(), 0);
	for (auto iter = m_pairs.begin(); iter != m_pairs.end(); ++iter)
	{
		const size_t type_1 = std::get<0>(*iter).Property();
		const size_t type_2 = std::get<1>(*iter).Property();
		const size_t type_min = std::min(type_1, type_2);
		const size_t type_max = std::max(type_1, type_2);
		const size_t dist_type = CalculateType(std::get<2>(*iter), m_dist_threshes);
		const size_t pair_type = (m_dist_threshes.size() + 1) * (points_types_num * type_max + type_min) + dist_type;
		m_pairs_with_types.push_back(std::make_tuple(
			std::get<0>(*iter), std::get<1>(*iter), pair_type));
		m_types.push_back(pair_type);
		m_pairs_histogram[pair_type]++;
	}

	m_sing_pts_pairs_with_props_and_types.first.clear();
	m_sing_pts_pairs_with_props_and_types.second.clear();
	for (size_t ind_1 = 0; ind_1 < m_non_marked_singular_points.first.size(); ++ind_1)
	{
		for (size_t ind_2 = ind_1 + 1; ind_2 < m_non_marked_singular_points.first.size(); ++ind_2)
		{
			const size_t type_1 = m_non_marked_singular_points.second[ind_1];
			const size_t type_2 = m_non_marked_singular_points.second[ind_2];
			const size_t type_min = std::min(type_1, type_2);
			const size_t type_max = std::max(type_1, type_2);
			const double curr_dist = cv::norm(m_non_marked_singular_points.first[ind_1].Coord() - m_non_marked_singular_points.first[ind_2].Coord());
			const size_t dist_type = CalculateType(curr_dist, m_dist_threshes);
			const size_t curr_pair_type = (m_dist_threshes.size() + 1) * (points_types_num * type_max + type_min) + dist_type;
			m_sing_pts_pairs_with_props_and_types.second.push_back(curr_pair_type);
			m_sing_pts_pairs_with_props_and_types.first.push_back(SingularPointsPair<PropertiesSet>());
			m_sing_pts_pairs_with_props_and_types.first.back().dist = curr_dist;
			m_sing_pts_pairs_with_props_and_types.first.back().elem1 = m_non_marked_singular_points.first[ind_1].Property();
			m_sing_pts_pairs_with_props_and_types.first.back().elem2 = m_non_marked_singular_points.first[ind_2].Property();
		}
	}

	CV_Assert(m_sing_pts_pairs_with_props_and_types.first.size() == m_sing_pts_pairs_with_props_and_types.second.size());
	CV_Assert(m_sing_pts_pairs_with_props_and_types.first.size() == m_types.size());
	int diff = 0;
	for (size_t ind = 0; ind < m_types.size(); ++ind)
	{
		CV_Assert(abs(static_cast<int>(m_sing_pts_pairs_with_props_and_types.second[ind] - m_types[ind])) <= 1);
	}
}

void MoleculeManager::CalculatePairsWithTypesLevels()
{
	CV_Assert(!m_dist_threshes_levels.empty());
	CV_Assert(m_dist_high_threshes.size() == m_dist_threshes_levels.size());
	const int levels_num = m_pairs_levels.size();
	m_pairs_with_types_levels.resize(levels_num);
	m_types_levels.resize(levels_num);
	const size_t pairs_max_type = GetPairsTypeNumLevels();
	m_pairs_histogram_levels.clear();
	m_pairs_histogram_levels.resize(pairs_max_type, 0);
	m_sing_pts_pairs_with_props_and_types_levels.first.resize(levels_num);
	m_sing_pts_pairs_with_props_and_types_levels.second.resize(levels_num);
	
	for (int level = 0; level < levels_num; ++level)
	{
		m_pairs_with_types_levels[level].clear();
		m_pairs_with_types_levels[level].reserve(m_pairs_levels[level].size());
		m_types_levels[level].clear();
		m_types_levels[level].reserve(m_pairs_levels[level].size());
		const size_t points_types_num = GetPointsTypeNum();

		for (auto iter = m_pairs_levels[level].begin(); iter != m_pairs_levels[level].end(); ++iter)
		{
			const double curr_dist = std::get<2>(*iter);
			if (curr_dist > m_dist_high_threshes[level])
			{
				continue;
			}
			const size_t type_1 = std::get<0>(*iter).Property();
			const size_t type_2 = std::get<1>(*iter).Property();
			const size_t type_min = std::min(type_1, type_2);
			const size_t type_max = std::max(type_1, type_2);
			const size_t dist_type = CalculateType(std::get<2>(*iter), m_dist_threshes_levels[level]);
			const size_t pair_type = levels_num * ((m_dist_threshes_levels[level].size() + 1) * (points_types_num * type_max + type_min) + dist_type) + level;
			m_pairs_with_types_levels[level].push_back(std::make_tuple(
				std::get<0>(*iter), std::get<1>(*iter), pair_type));
			m_types_levels[level].push_back(pair_type);
			m_pairs_histogram_levels[pair_type] += m_levels_scales[level];
		}

		/*m_sing_pts_pairs_with_props_and_types_levels.first[level].clear();
		m_sing_pts_pairs_with_props_and_types_levels.second[level].clear();

		for (size_t ind_1 = 0; ind_1 < m_non_marked_singular_points_levels.first[level].size(); ++ind_1)
		{
			for (size_t ind_2 = ind_1 + 1; ind_2 < m_non_marked_singular_points_levels.first[level].size(); ++ind_2)
			{
				const size_t type_1 = m_non_marked_singular_points_levels.second[level][ind_1];
				const size_t type_2 = m_non_marked_singular_points_levels.second[level][ind_2];
				const size_t type_min = std::min(type_1, type_2);
				const size_t type_max = std::max(type_1, type_2);
				const double curr_dist = cv::norm(
					m_non_marked_singular_points_levels.first[level][ind_1].Coord() - m_non_marked_singular_points_levels.first[level][ind_2].Coord());
				if (curr_dist > m_dist_high_threshes[level])
				{
					continue;
				}
				const size_t dist_type = CalculateType(curr_dist, m_dist_threshes);
				const size_t curr_pair_type = levels_num * ((m_dist_threshes_levels[level].size() + 1) * (points_types_num * type_max + type_min) + dist_type) + level;
				m_sing_pts_pairs_with_props_and_types_levels.second[level].push_back(curr_pair_type);
				m_sing_pts_pairs_with_props_and_types_levels.first[level].push_back(SingularPointsPair<PropertiesSet>());
				m_sing_pts_pairs_with_props_and_types_levels.first[level].back().dist = curr_dist;
				m_sing_pts_pairs_with_props_and_types_levels.first[level].back().elem1 = m_non_marked_singular_points_levels.first[level][ind_1].Property();
				m_sing_pts_pairs_with_props_and_types_levels.first[level].back().elem2 = m_non_marked_singular_points_levels.first[level][ind_2].Property();
			}
		}

		CV_Assert(m_sing_pts_pairs_with_props_and_types_levels.first[level].size() == m_sing_pts_pairs_with_props_and_types_levels.second[level].size());
		CV_Assert(m_sing_pts_pairs_with_props_and_types_levels.first[level].size() == m_types_levels[level].size());
		int diff = 0;
		for (size_t ind = 0; ind < m_types_levels[level].size(); ++ind)
		{
			CV_Assert(abs(static_cast<int>(m_sing_pts_pairs_with_props_and_types_levels.second[level][ind] - m_types_levels[level][ind])) <= 1);
		}*/
	}
}

void MoleculeManager::CalculateTriplesWithTypes()
{
	CV_Assert(!m_dist_threshes.empty());
	m_triples_with_types.clear();
	m_triples_with_types.reserve(m_triples.size());
	const size_t points_max_type = GetPointsTypeNum();
	const size_t triples_max_type = GetTriplesTypeNum();
	m_triples_histogram.resize(triples_max_type);
	std::fill(m_triples_histogram.begin(), m_triples_histogram.end(), 0);
	//containers for types calculation
	std::vector<size_t> triangle_max_types(6);
	triangle_max_types[0] = triangle_max_types[2] = triangle_max_types[4] = points_max_type;
	triangle_max_types[1] = triangle_max_types[3] = triangle_max_types[5] = m_dist_threshes.size() + 1;
	std::vector<size_t> triangle_types(6);

	for (auto iter = m_triples.begin(); iter != m_triples.end(); ++iter)
	{
		const size_t type_1 = triangle_types[0] = std::get<0>(*iter).first.Property();
		const size_t type_2 = triangle_types[2] = std::get<1>(*iter).first.Property();
		const size_t type_3 = triangle_types[4] = std::get<2>(*iter).first.Property();

		const size_t dist_type_1 = triangle_types[1] = CalculateType(std::get<0>(*iter).second, m_dist_threshes);
		const size_t dist_type_2 = triangle_types[3] = CalculateType(std::get<1>(*iter).second, m_dist_threshes);
		const size_t dist_type_3 = triangle_types[5] = CalculateType(std::get<2>(*iter).second, m_dist_threshes);

		const size_t triple_type = CalculateRangeType(triangle_types.begin(), triangle_types.end(), 
			triangle_max_types.begin(), triangle_max_types.end());
		m_triples_with_types.push_back(std::make_tuple(
			std::make_pair(std::get<0>(*iter).first, dist_type_1), 
			std::make_pair(std::get<1>(*iter).first, dist_type_2),
			std::make_pair(std::get<2>(*iter).first, dist_type_3)));
		m_triples_histogram[triple_type]++;
	}
}

void MoleculeManager::GetPairsHistogramm(cv::Mat_<size_t>& histogram)
{
	CV_Assert(histogram.rows == 1);
	CV_Assert(histogram.cols == m_pairs_histogram.size());

	std::copy(m_pairs_histogram.begin(), m_pairs_histogram.end(), histogram.begin());
}
void MoleculeManager::GetPairsHistogrammLevels(cv::Mat_<double>& histogram)
{
	CV_Assert(histogram.rows == 1);
	CV_Assert(histogram.cols == m_pairs_histogram_levels.size());

	std::copy(m_pairs_histogram_levels.begin(), m_pairs_histogram_levels.end(), histogram.begin());
}

void MoleculeManager::GetTriplesHistogramm(cv::Mat_<size_t>& histogram)
{
	CV_Assert(histogram.rows == 1);
	CV_Assert(histogram.cols == m_triples_histogram.size());

	std::copy(m_triples_histogram.begin(), m_triples_histogram.end(), histogram.begin());
}

size_t MoleculeManager::GetPairsTypeNum()
{
	CV_Assert(!m_dist_threshes.empty());
		const size_t points_max_type = GetPointsTypeNum();	
		return (m_dist_threshes.size() + 1) * points_max_type * points_max_type;
}

size_t MoleculeManager::GetPairsTypeNumLevels()
{
	CV_Assert(!m_dist_threshes_levels.empty());
	const size_t points_max_type = GetPointsTypeNum();	
	return m_dist_threshes_levels.size() * (m_dist_threshes_levels[0].size() + 1) * points_max_type * points_max_type;
}

size_t MoleculeManager::GetTriplesTypeNum()
{
	CV_Assert(!m_dist_threshes.empty());
	const size_t points_max_type = GetPointsTypeNum();	
	return pow(m_dist_threshes.size() + 1, 3) * pow(points_max_type, 3);
}

size_t MoleculeManager::GetPointsTypeNum()
{
	CV_Assert(!m_lennard_jones_threshes.empty());
	CV_Assert(!m_charge_threshes.empty());
	return kSaddleType * (m_lennard_jones_threshes.size() + 1) * (m_charge_threshes.size() + 1) * (m_area_threshes.size() + 1);
}

void MoleculeManager::Clear()
{
	m_singular_points.clear();
	m_pairs.clear();
	m_pairs_with_types.clear();
	m_types.clear();
	m_pairs_histogram.clear();
}

void MoleculeManager::CalculateTriples()
{
	m_triples.clear();
	m_triples.reserve(m_singular_points.size() * m_singular_points.size() * m_singular_points.size());

	for (auto iter1 = m_singular_points.begin(); iter1 != m_singular_points.end(); ++iter1)
	{
		for (auto iter2 = iter1 + 1; iter2 != m_singular_points.end(); ++iter2)
		{
			const double dist_1 = cv::norm(iter1->Coord() - iter2->Coord());

			for (auto iter3 = iter2 + 1; iter3 != m_singular_points.end(); ++iter3)
			{
				const double dist_2 = cv::norm(iter1->Coord() - iter3->Coord());
				const double dist_3 = cv::norm(iter2->Coord() - iter3->Coord());
				m_triples.push_back(std::make_tuple(std::make_pair(*iter1, dist_1),
					std::make_pair(*iter2, dist_2), std::make_pair(*iter3, dist_3)));
			}
		}
	}
}

void MoleculeManager::WriteTriples()
{
	const std::string triples_file_name = m_curr_file_prefix + Extensions::Triples();
	WriteInterval(triples_file_name, m_triples.begin(), m_triples.end());
}

void MoleculeManager::ReadTriples()
{
	const std::string triples_file_name = m_curr_file_prefix + Extensions::Triples();
	ReadVector(triples_file_name, m_triples);
}

void MoleculeManager::CalculatePairs()
{
	CV_Assert(!m_singular_points.empty());
	m_pairs.clear();
	m_pairs.reserve(
		(m_singular_points.size() * (m_singular_points.size() - 1)) / 2);

	for (auto iter1 = m_singular_points.begin(); iter1 != m_singular_points.end(); ++iter1)
	{
		for (auto iter2 = iter1 + 1; iter2 != m_singular_points.end(); ++iter2)
		{
			const double curr_dist = cv::norm(iter1->Coord() - iter2->Coord());
			m_pairs.push_back(std::make_tuple(*iter1, *iter2, curr_dist));
		}
	}
}

void MoleculeManager::ReadPairs()
{
	const std::string pairs_file_name = m_curr_file_prefix + Extensions::Pairs();
	ReadVector(pairs_file_name, m_pairs);
}

void MoleculeManager::WritePairs()
{
	const std::string pairs_file_name = m_curr_file_prefix + Extensions::Pairs();
	WriteInterval(pairs_file_name, m_pairs.begin(), m_pairs.end());
}

void MoleculeManager::CalculatePairsLevels()
{
	const int levels_num = m_non_marked_singular_points_levels.first.size();
	m_pairs_levels.resize(levels_num);

	for (int level = 0; level < levels_num; ++level)
	{
		m_pairs_levels[level].clear();
		for (size_t ind1 = 0; ind1 < m_non_marked_singular_points_levels.first[level].size(); ++ind1)
		{
			for (size_t ind2 = ind1 + 1; ind2 < m_non_marked_singular_points_levels.first[level].size(); ++ind2)
			{
				pair_with_distance new_pair;
				std::get<0>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[level][ind1].Coord(), 
					m_non_marked_singular_points_levels.second[level][ind1]);
				std::get<1>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[level][ind2].Coord(), 
					m_non_marked_singular_points_levels.second[level][ind2]);
				std::get<2>(new_pair) = cv::norm(m_non_marked_singular_points_levels.first[level][ind1].Coord() - 
					m_non_marked_singular_points_levels.first[level][ind2].Coord());
				m_pairs_levels[level].push_back(new_pair);
			}
		}
	}
}
void MoleculeManager::WritePairsLevels()
{
	for (size_t level = 0; level < m_pairs_levels.size(); ++level)
	{		
		const std::string pairs_file_name_lev = m_curr_file_prefix + Extensions::PairsLevels() + std::to_string(level);
		WriteInterval(pairs_file_name_lev, m_pairs_levels[level].begin(), m_pairs_levels[level].end());		
	}
}
void MoleculeManager::ReadPairsLevels(const int levels_num)
{
	////find out how many levels are there
	//std::ifstream test_stream;
	//int levels_num = -1; 
	//do
	//{	
	//	test_stream.close();
	//	++levels_num;
	//	const std::string pairs_file_name_lev = m_curr_file_prefix + Extensions::PairsLevels() + std::to_string(levels_num);
	//	test_stream.open(pairs_file_name_lev);
	//}
	//while(test_stream.is_open());

	CV_Assert(levels_num > 0);
	m_pairs_levels.resize(levels_num);
	//read levels
	for (int level = 0; level < levels_num; ++level)
	{
		const std::string pairs_file_name_lev = m_curr_file_prefix + Extensions::PairsLevels() + std::to_string(level);
		ReadVector(pairs_file_name_lev, m_pairs_levels[level]);			
	}
}
int MoleculeManager::CalculateSingularPointsTypes(PropertiesSet& prop)
{
	IntRange curvature_range;
	curvature_range.min_val = kConvexType;
	curvature_range.max_val = kSaddleType;
	curvature_range.step = 1;

	IntRange lennard_jones_range;
	lennard_jones_range.min_val = 0;
	lennard_jones_range.max_val = m_lennard_jones_threshes.size();
	lennard_jones_range.step = 1; 

	IntRange charge_range;
	charge_range.min_val = 0;
	charge_range.max_val = m_charge_threshes.size();
	charge_range.step = 1; 

	IntRange area_range;
	area_range.min_val = 0;
	area_range.max_val = m_area_threshes.size();
	area_range.step = 1; 

	std::vector<TypeWithMax> type_with_max(4);
	enum {SIGN, LENNARD_JONES, CURVATURE, AREA};
	type_with_max[SIGN].range = charge_range;
	type_with_max[LENNARD_JONES].range = lennard_jones_range;
	type_with_max[CURVATURE].range = curvature_range;
	type_with_max[AREA].range = area_range;

	//const size_t curr_sign = Sign(m_vertex_charge_map[*iter]) + 1;
	const int curr_sign = CalculateType(static_cast<double>(prop.ElectricPotential()), m_charge_threshes);
	prop.Charge() = curr_sign;
	type_with_max[SIGN].type = curr_sign;

	const size_t curr_lennard_jones = CalculateType(static_cast<double>(prop.LennardJones()), m_lennard_jones_threshes);
	prop.LennardJonesType() = curr_lennard_jones;
	type_with_max[LENNARD_JONES].type = curr_lennard_jones;

	const size_t curr_area_type = CalculateType(static_cast<double>(prop.Area()), m_area_threshes);
	prop.Area() = curr_area_type;
	type_with_max[AREA].type = curr_area_type;

	type_with_max[CURVATURE].type = prop.SurfaceType();
	
	const size_t all_type = CalculateTypesComposition(type_with_max);
	CV_Assert(all_type < GetPointsTypeNum());
	return all_type;
}

void MoleculeManager::CalculatePropertiesTypes()
{
	m_non_marked_singular_points.second.resize(m_non_marked_singular_points.first.size());
	m_singular_points.resize(m_non_marked_singular_points.first.size());
	for (size_t ind = 0; ind < m_non_marked_singular_points.first.size(); ++ind)
	{
		m_non_marked_singular_points.second[ind] = CalculateSingularPointsTypes(m_non_marked_singular_points.first[ind].Property());
		m_singular_points[ind].Property() = m_non_marked_singular_points.second[ind];
		m_singular_points[ind].Coord() = m_non_marked_singular_points.first[ind].Coord();
	}

	const int levels_num = m_non_marked_singular_points_levels.first.size();
	m_non_marked_singular_points_levels.second.resize(levels_num);
	for (int level = 0; level < levels_num; ++level)
	{
		const size_t curr_points_num = m_non_marked_singular_points_levels.first[level].size();
		m_non_marked_singular_points_levels.second[level].resize(curr_points_num);
		for (size_t ind = 0; ind < m_non_marked_singular_points_levels.first[level]	.size(); ++ind)
		{
			m_non_marked_singular_points_levels.second[level][ind] = CalculateSingularPointsTypes(m_non_marked_singular_points_levels.first[level][ind].Property());
		}
		//возможно добавить вычисление нетипизированных особых точек для кернела
	}
}

void MoleculeManager::CalculateSingularPoints(const bool calc_as_average)
{
	//read surface
	const std::string surface_file_name = m_curr_file_prefix + Extensions::Surface();
	SurfaceReader surface_reader;
	CV_Assert(surface_reader.OpenFile(surface_file_name));
	std::vector<cv::Point3d> vertices;
	std::vector<cv::Point3d> normals;
	std::vector<cv::Point3i> triangles;

	surface_reader.ReadVertices(vertices);
	surface_reader.ReadNormals(normals);
	surface_reader.ReadTriangles(triangles);
	WriteTriangles(triangles);
	//read charges
	const std::string charges_file_name = m_curr_file_prefix + Extensions::Charges();
	std::vector<std::pair<cv::Point3d, double>> charges;
	ReadVector(charges_file_name, charges);
	//read wdv radii
	const std::string radius_file_name = m_curr_file_prefix + Extensions::WDV();
	std::vector<std::pair<cv::Point3d, double>> wdv_radii;
	ReadVector(radius_file_name, wdv_radii);
	
	m_sing_pts_finder->Process(vertices, normals, triangles, charges, wdv_radii, calc_as_average, m_levels_num);
	m_sing_pts_finder->GetMarkedSingularPoints(m_singular_points);
	m_sing_pts_finder->GetNonMarkedSingularPoints(m_non_marked_singular_points);
	m_sing_pts_finder->GetNonMarkedSingularPointsLevels(m_non_marked_singular_points_levels);
	m_sing_pts_finder->GetSingularPointsHisto(m_histogram_singular_points);
}

void MoleculeManager::ReadAllSingularPoints(const int levels_num)
{
	const std::string points_file_name = m_curr_file_prefix + Extensions::SingPts();
	ReadVector(points_file_name, m_singular_points);

	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
	ReadVector(non_marked_points_file_name, m_non_marked_singular_points.first);
	const std::string non_marked_points_file_name_l = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabels();
	ReadVector(non_marked_points_file_name_l, m_non_marked_singular_points.second);

	const std::string hist_points_file_name = m_curr_file_prefix + Extensions::HistSingPts();
	ReadVector(hist_points_file_name, m_histogram_singular_points);

	ReadSingularPointsLevelsLabels(levels_num);
	ReadSingularPointsLevelsTypes(levels_num);
}

void MoleculeManager::WriteSingularPoints()
{
	const std::string points_file_name = m_curr_file_prefix + Extensions::SingPts();
	WriteInterval(points_file_name, m_singular_points.begin(), m_singular_points.end());

	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
	WriteInterval(non_marked_points_file_name, m_non_marked_singular_points.first.begin(), m_non_marked_singular_points.first.end());
	const std::string non_marked_points_file_name_l = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabels();
	WriteInterval(non_marked_points_file_name_l, m_non_marked_singular_points.second.begin(), m_non_marked_singular_points.second.end());

	for (size_t level = 0; level < m_non_marked_singular_points_levels.first.size(); ++level)
	{		
		const std::string non_marked_points_file_name_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypesLevels() + std::to_string(level);
		WriteInterval(non_marked_points_file_name_lev, 
			m_non_marked_singular_points_levels.first[level].begin(), m_non_marked_singular_points_levels.first[level].end());
		const std::string non_marked_points_file_name_l_lev = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabelsLevels() + std::to_string(level);
		WriteInterval(non_marked_points_file_name_l_lev, 
			m_non_marked_singular_points_levels.second[level].begin(), m_non_marked_singular_points_levels.second[level].end());
	}
	const std::string hist_points_file_name = m_curr_file_prefix + Extensions::HistSingPts();
	WriteInterval(hist_points_file_name, m_histogram_singular_points.begin(), m_histogram_singular_points.end());
}

void MoleculeManager::WriteSegmentedSurface()
{
	std::vector<std::pair<cv::Point3d, size_t>> segmented_vertices;
	m_sing_pts_finder->GetSegmentedVertices(segmented_vertices);
	const std::string segm_surface_file = m_curr_file_prefix + Extensions::SegmSurf();
	WriteInterval(segm_surface_file, segmented_vertices.begin(), segmented_vertices.end());
}

void MoleculeManager::WriteTriangles(const std::vector<cv::Point3i>& triangles)
{
	const std::string triangles_file = m_curr_file_prefix + Extensions::MeshTriangles();
	WriteInterval(triangles_file, triangles.begin(), triangles.end());
}

void MoleculeManager::WriteSurfaceWithTypes()
{
	std::vector<std::pair<cv::Point3d, size_t>> vertices_with_types;
	m_sing_pts_finder->GetVerticesWithTypes(vertices_with_types);
	const std::string vert_with_types_file = m_curr_file_prefix + Extensions::SurfWithTypes();
	WriteInterval(vert_with_types_file, vertices_with_types.begin(), vertices_with_types.end());
}

}
