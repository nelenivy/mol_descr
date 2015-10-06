#include "molecule_manager.h"

#include <algorithm>
#include <iostream>

#include "extensions.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"
#include "types_calculation.h"
#include "../SingularPoints/sng_pts_finder_scale_space.h"
#include "SingularPoints/mesh_types.h"
#include "SingularPoints\mesh_operations.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"
#include "GraphLib\graph_filter.h"
#include "GraphLib\curvature_calculator.h"
#include "GraphLib\proxy_property_map.h"
#include "GraphLib\graph_functions.h"
#include "CommonUtilities\common_functions.h"
#include "SingularPoints/points_keeper.h"
#include "SingularPoints/properties_calculator.h"
#include "../scale_space_blur.h"

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
		WriteSurfaceWithTypesLevels();

		for (int prop_type = ISingularPointsFinder::FIRST_SURF_PROP; prop_type < ISingularPointsFinder::SURF_PROPS_NUM; ++prop_type)
		{
			WriteSurfaceWithDblProp(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
		}

		if (!m_use_calculated_scale_space)
		{
			for (int prop_type = ISingularPointsFinder::FIRST_SURF_PROP; prop_type < ISingularPointsFinder::SURF_PROPS_NUM; ++prop_type)
			{
				WriteSurfaceWithDblProp(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
				WriteSurfaceWithDblPropLevels(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
			}
		}

		if (!m_use_calculated_detector_function)
		{
			for (int prop_type = ISingularPointsFinder::FIRST_LOG_PROP; prop_type <= ISingularPointsFinder::LAST_LOG_PROP; ++prop_type)
			{
				WriteSurfaceWithDblPropLevels(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
			}
		}

		if (!m_use_calculated_eig_ratio)
		{
			for (int prop_type = ISingularPointsFinder::FIRST_PCA_PROP; prop_type <= ISingularPointsFinder::LAST_PCA_PROP; ++prop_type)
			{
				WriteSurfaceWithDblPropLevels(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
			}
			for (int prop_type = ISingularPointsFinder::FIRST_EIG_PROP; prop_type <= ISingularPointsFinder::LAST_EIG_PROP; ++prop_type)
			{
				WriteSurfaceWithDblPropLevels(static_cast<ISingularPointsFinder::SurfProperty>(prop_type));
			}
			WriteSurfaceWithDblPropLevels(ISingularPointsFinder::PCA_SCALE_SPACE);
			WriteSurfaceWithDblPropLevels(ISingularPointsFinder::PCA_GRAD);
			WriteSurfaceWithDblPropLevels(ISingularPointsFinder::PCA_EIG);
			WriteSurfaceWithDblPropLevels(ISingularPointsFinder::PCA_EIG_LOG);
		}
		WriteSurfaceWithDblPropLevels(ISingularPointsFinder::PCA_LOG);
		WriteSurfaceWithDblProp(ISingularPointsFinder::UNCOIN_BASIS);
	}
	else
	{
		const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
		ReadVector(non_marked_points_file_name, m_non_marked_singular_points.first);
		ReadSingularPointsLevelsTypes(m_sing_pts_levels_num);
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
	ReadSingularPointsLevelsTypes(m_sing_pts_levels_num);

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
void MoleculeManager::FindPairsLevels(const bool calculate, const bool write)
{
	m_pairs_levels_num = m_sing_pts_levels_num - m_pairs_levels_overlap;

	if (calculate)
	{
		CalculatePairsLevels();

		if (write)
		{
			WritePairsLevels();
		}
	}
	else
	{
		ReadPairsLevels(m_pairs_levels_num);
	}
}

void MoleculeManager::FindTriplesLevels(const bool calculate, const bool write)
{
	CalculateTriplesLevels();
	/*m_pairs_levels_num = m_sing_pts_levels_num - m_pairs_levels_overlap;

	if (calculate)
	{
		CalculatePairsLevels();

		if (write)
		{
			WritePairsLevels();
		}
	}
	else
	{
		ReadPairsLevels(m_pairs_levels_num);
	}*/
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
	CV_Assert(m_pairs_levels.size() == m_pairs_levels_num);

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

void MoleculeManager::CalculatePairsWithTypesLevelsAllMesh()
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
	
	//read the surface and compile pairs of all points
	const std::string surface_file_name = m_curr_file_prefix + Extensions::Surface();
	SurfaceReader surface_reader;
	CV_Assert(surface_reader.OpenFile(surface_file_name));
	std::vector<cv::Point3d> vertices;
	std::vector<cv::Point3d> normals;
	std::vector<cv::Point3i> triangles;

	surface_reader.ReadVertices(vertices);
	surface_reader.ReadNormals(normals);
	surface_reader.ReadTriangles(triangles);
	//read charges
	const std::string charges_file_name = m_curr_file_prefix + Extensions::Charges();
	std::vector<std::pair<cv::Point3d, double>> charges;
	ReadVector(charges_file_name, charges);
	//read wdv radii
	const std::string radius_file_name = m_curr_file_prefix + Extensions::WDV();
	std::vector<std::pair<cv::Point3d, double>> wdv_radii;
	ReadVector(radius_file_name, wdv_radii);
	MeshKeeper mesh_keeper;
	mesh_keeper.ConstructMesh(vertices, normals, triangles);
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	VetrticesChargeMap vertex_charge_map;
	CalculateAllPotentials(charges, mesh_keeper.GetMesh()/*mesh_to_use*/, vertex_charge_map);
	VetrticesChargeMap vertex_lennard_jones_map;
	CalculateLennardJonesPotentials(wdv_radii, mesh_keeper.GetMesh()/*mesh_to_use*/, vertex_lennard_jones_map);
	std::vector<std::pair<cv::Point3d, size_t>> vertices_with_types;
	ReadSurfaceWithTypes(vertices_with_types);
	std::vector<pair_with_distance> pairs_vertices;

	for (size_t ind1 = 0; ind1 < vertices_with_types.size(); ++ind1)
	{
		PropertiesSet prop1;
		prop1.ElectricPotential() = vertex_charge_map[ind1];
		prop1.LennardJones() = vertex_lennard_jones_map[ind1];
		prop1.SurfaceType() = vertices_with_types[ind1].second;
		const size_t p1_type = CalculateSingularPointsTypes(prop1);
		MarkedSingularPoint point1;
		point1.Coord() = vertices_with_types[ind1].first;
		point1.Property() = p1_type;
		for (size_t ind2 = ind1 + 1; ind2 < vertices_with_types.size(); ++ind2)
		{
			PropertiesSet prop2;
			prop2.ElectricPotential() = vertex_charge_map[ind2];
			prop2.LennardJones() = vertex_lennard_jones_map[ind2];
			prop2.SurfaceType() = vertices_with_types[ind2].second;
			const size_t p2_type = CalculateSingularPointsTypes(prop2);

			MarkedSingularPoint point2;
			point2.Coord() = vertices_with_types[ind2].first;
			point2.Property() = p2_type;

			pairs_vertices.push_back(std::make_tuple(point1, point2, cv::norm(point2.Coord() - point1.Coord())));
		}
	}

	for (int level = 0; level < levels_num; ++level)
	{
		m_pairs_with_types_levels[level].clear();
		m_pairs_with_types_levels[level].reserve(m_pairs_levels[level].size());
		m_types_levels[level].clear();
		m_types_levels[level].reserve(m_pairs_levels[level].size());
		const size_t points_types_num = GetPointsTypeNum();

		for (auto iter = pairs_vertices.begin(); iter != pairs_vertices.end(); ++iter)
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

void MoleculeManager::CalculatePairsWithTypesLevelsAllMeshSmoothedCurv()
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
	
	//read the surface and compile pairs of all points
	const std::string surface_file_name = m_curr_file_prefix + Extensions::Surface();
	SurfaceReader surface_reader;
	CV_Assert(surface_reader.OpenFile(surface_file_name));
	std::vector<cv::Point3d> vertices;
	std::vector<cv::Point3d> normals;
	std::vector<cv::Point3i> triangles;

	surface_reader.ReadVertices(vertices);
	surface_reader.ReadNormals(normals);
	surface_reader.ReadTriangles(triangles);
	//read charges
	const std::string charges_file_name = m_curr_file_prefix + Extensions::Charges();
	std::vector<std::pair<cv::Point3d, double>> charges;
	ReadVector(charges_file_name, charges);
	//read wdv radii
	const std::string radius_file_name = m_curr_file_prefix + Extensions::WDV();
	std::vector<std::pair<cv::Point3d, double>> wdv_radii;
	ReadVector(radius_file_name, wdv_radii);
	MeshKeeper mesh_keeper;
	mesh_keeper.ConstructMesh(vertices, normals, triangles);
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	VetrticesChargeMap vertex_charge_map;
	CalculateAllPotentials(charges, mesh_keeper.GetMesh()/*mesh_to_use*/, vertex_charge_map);
	VetrticesChargeMap vertex_lennard_jones_map;
	CalculateLennardJonesPotentials(wdv_radii, mesh_keeper.GetMesh()/*mesh_to_use*/, vertex_lennard_jones_map);
	std::vector<std::vector<std::pair<cv::Point3d, size_t>>> vertices_with_types_lev;
	ReadSurfaceWithTypesLevels(vertices_with_types_lev);

	for (int level = 0; level < levels_num; ++level)
	{
		m_pairs_with_types_levels[level].clear();
		m_pairs_with_types_levels[level].reserve(m_pairs_levels[level].size());
		m_types_levels[level].clear();
		m_types_levels[level].reserve(m_pairs_levels[level].size());
		const size_t points_types_num = GetPointsTypeNum();

		for (size_t ind1 = 0; ind1 < vertices_with_types_lev[level].size(); ++ind1)
		{
			PropertiesSet prop1;
			prop1.ElectricPotential() = vertex_charge_map[ind1];
			prop1.LennardJones() = vertex_lennard_jones_map[ind1];
			prop1.SurfaceType() = vertices_with_types_lev[level][ind1].second;
			const size_t p1_type = CalculateSingularPointsTypes(prop1);
			MarkedSingularPoint point1;
			point1.Coord() = vertices_with_types_lev[level][ind1].first;
			point1.Property() = p1_type;
			for (size_t ind2 = ind1 + 1; ind2 < vertices_with_types_lev[level].size(); ++ind2)
			{
				PropertiesSet prop2;
				prop2.ElectricPotential() = vertex_charge_map[ind2];
				prop2.LennardJones() = vertex_lennard_jones_map[ind2];
				prop2.SurfaceType() = vertices_with_types_lev[level][ind2].second;
				const size_t p2_type = CalculateSingularPointsTypes(prop2);

				MarkedSingularPoint point2;
				point2.Coord() = vertices_with_types_lev[level][ind2].first;
				point2.Property() = p2_type;

				const double curr_dist = cv::norm(point2.Coord() - point1.Coord());
				if (curr_dist > m_dist_high_threshes[level])
				{
					continue;
				}
			
				const size_t type_min = std::min(p2_type, p1_type);
				const size_t type_max = std::max(p2_type, p1_type);
				const size_t dist_type = CalculateType(curr_dist, m_dist_threshes_levels[level]);
				const size_t pair_type = levels_num * ((m_dist_threshes_levels[level].size() + 1) * (points_types_num * type_max + type_min) + dist_type) + level;
				m_pairs_with_types_levels[level].push_back(std::make_tuple(point1, point2, pair_type));
				m_types_levels[level].push_back(pair_type);
				m_pairs_histogram_levels[pair_type] += m_levels_scales[level];
			}
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

void MoleculeManager::CalculatePairsWithTypesLevels()
{
	CV_Assert(!m_dist_threshes_levels.empty());
	CV_Assert(m_dist_high_threshes.size() == m_dist_threshes_levels.size());
	const int levels_num = m_pairs_levels.size();
	CV_Assert(levels_num == m_pairs_levels_num);

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

void MoleculeManager::CalculateTriplesWithTypesLevels()
{
	CV_Assert(!m_dist_threshes_levels.empty());
	CV_Assert(!m_dist_threshes_levels[0].empty());
	CV_Assert(m_dist_threshes_levels.size() == m_pairs_levels_num);

	m_triples_with_types_levels.clear();
	m_triples_with_types_levels.reserve(m_triples_levels.size());
	const size_t points_max_type = GetPointsTypeNum();
	const size_t triples_max_type_levels = GetTriplesTypeNumLevels();
	m_triples_histogram_levels.resize(triples_max_type_levels);
	std::fill(m_triples_histogram_levels.begin(), m_triples_histogram_levels.end(), 0);
	//containers for types calculation
	std::vector<size_t> triangle_max_types_levels(7);
	triangle_max_types_levels[0] = triangle_max_types_levels[2] = triangle_max_types_levels[4] = points_max_type;
	triangle_max_types_levels[1] = triangle_max_types_levels[3] = triangle_max_types_levels[5] = m_dist_threshes_levels[0].size() + 1;
	triangle_max_types_levels[6] = m_pairs_levels_num;
	std::vector<size_t> triangle_types(7);

	m_triples_with_types_levels.resize(m_pairs_levels_num);

	for (size_t curr_level = 0; curr_level < m_pairs_levels_num; ++curr_level)
	{
		for (auto iter = m_triples_levels[curr_level].begin(); iter != m_triples_levels[curr_level].end(); ++iter)
		{
			const size_t type_1 = triangle_types[0] = std::get<0>(*iter).first.Property();
			const size_t type_2 = triangle_types[2] = std::get<1>(*iter).first.Property();
			const size_t type_3 = triangle_types[4] = std::get<2>(*iter).first.Property();

			const size_t dist_type_1 = triangle_types[1] = CalculateType(std::get<0>(*iter).second, m_dist_threshes_levels[curr_level]);
			const size_t dist_type_2 = triangle_types[3] = CalculateType(std::get<1>(*iter).second, m_dist_threshes_levels[curr_level]);
			const size_t dist_type_3 = triangle_types[5] = CalculateType(std::get<2>(*iter).second, m_dist_threshes_levels[curr_level]);
			triangle_types[6] = curr_level;

			const size_t triple_type = CalculateRangeType(triangle_types.begin(), triangle_types.end(), 
				triangle_max_types_levels.begin(), triangle_max_types_levels.end());
			m_triples_with_types_levels[curr_level].push_back(std::make_tuple(
				std::make_pair(std::get<0>(*iter).first, dist_type_1), 
				std::make_pair(std::get<1>(*iter).first, dist_type_2),
				std::make_pair(std::get<2>(*iter).first, dist_type_3)));
			m_triples_histogram_levels[triple_type]++;
		}
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

void MoleculeManager::GetTriplesHistogrammLevels(cv::Mat_<float>& histogram)
{
	CV_Assert(histogram.rows == 1);
	CV_Assert(histogram.cols == m_triples_histogram_levels.size());

	std::copy(m_triples_histogram_levels.begin(), m_triples_histogram_levels.end(), histogram.begin());
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
	CV_Assert(m_dist_threshes_levels.size() == m_pairs_levels_num);
	return m_dist_threshes_levels.size() * (m_dist_threshes_levels[0].size() + 1) * points_max_type * points_max_type;
}

size_t MoleculeManager::GetTriplesTypeNum()
{
	CV_Assert(!m_dist_threshes.empty());
	const size_t points_max_type = GetPointsTypeNum();	
	return pow(m_dist_threshes.size() + 1, 3) * pow(points_max_type, 3);
}
size_t MoleculeManager::GetTriplesTypeNumLevels()
{
	CV_Assert(!m_dist_threshes_levels.empty());
	const size_t points_max_type = GetPointsTypeNum();	
	CV_Assert(m_dist_threshes_levels.size() == m_pairs_levels_num);
	return m_dist_threshes_levels.size() * pow(m_dist_threshes_levels[0].size() + 1, 3) * pow(points_max_type, 3);
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
	CV_Assert(levels_num == m_sing_pts_levels_num);
	CV_Assert(m_pairs_levels_num > 0);
	m_pairs_levels.resize(m_pairs_levels_num);

	for (int level = m_pairs_levels_overlap; level < m_sing_pts_levels_num; ++level)
	{
		const int pair_level = level - m_pairs_levels_overlap;
		m_pairs_levels[pair_level].clear();

		//for simplicity put all points in one container
		std::vector<NonMarkedSingularPoint> points;
		std::vector<size_t> types;
		for (int curr_level = level - m_pairs_levels_overlap; curr_level <= level; ++curr_level)
		{
			points.insert(points.end(), m_non_marked_singular_points_levels.first[curr_level].begin(), 
				m_non_marked_singular_points_levels.first[curr_level].end());
			types.insert(types.end(), m_non_marked_singular_points_levels.second[curr_level].begin(), 
				m_non_marked_singular_points_levels.second[curr_level].end());
		}

		for (size_t ind1 = 0; ind1 < points.size(); ++ind1)
		{
			for (size_t ind2 = ind1 + 1; ind2 < points.size(); ++ind2)
			{
				pair_with_distance new_pair;
				std::get<0>(new_pair) = MarkedSingularPoint(points[ind1].Coord(), types[ind1]);
				std::get<1>(new_pair) = MarkedSingularPoint(points[ind2].Coord(), types[ind2]);
				std::get<2>(new_pair) = cv::norm(points[ind1].Coord() - points[ind2].Coord());
				m_pairs_levels[pair_level].push_back(new_pair);
			}
		}
		
		/*for (int curr_level = level - m_pairs_levels_overlap; curr_level <= level; ++curr_level)
		{
			for (size_t ind1 = 0; ind1 < m_non_marked_singular_points_levels.first[curr_level].size(); ++ind1)
			{
				for (size_t ind2 = ind1 + 1; ind2 < m_non_marked_singular_points_levels.first[curr_level].size(); ++ind2)
				{
					pair_with_distance new_pair;
					std::get<0>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[curr_level][ind1].Coord(), 
						m_non_marked_singular_points_levels.second[curr_level][ind1]);
					std::get<1>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[curr_level][ind2].Coord(), 
						m_non_marked_singular_points_levels.second[curr_level][ind2]);
					std::get<2>(new_pair) = cv::norm(m_non_marked_singular_points_levels.first[curr_level][ind1].Coord() - 
						m_non_marked_singular_points_levels.first[curr_level][ind2].Coord());
					m_pairs_levels[pair_level].push_back(new_pair);
				}

				for (int next_level = curr_level + 1; next_level <= level; ++next_level)
				{
					for (size_t ind2 = 0; ind2 < m_non_marked_singular_points_levels.first[next_level].size(); ++ind2)
					{
						pair_with_distance new_pair;
						std::get<0>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[curr_level][ind1].Coord(), 
							m_non_marked_singular_points_levels.second[curr_level][ind1]);
						std::get<1>(new_pair) = MarkedSingularPoint(m_non_marked_singular_points_levels.first[next_level][ind2].Coord(), 
							m_non_marked_singular_points_levels.second[next_level][ind2]);
						std::get<2>(new_pair) = cv::norm(m_non_marked_singular_points_levels.first[curr_level][ind1].Coord() - 
							m_non_marked_singular_points_levels.first[next_level][ind2].Coord());
						m_pairs_levels[pair_level].push_back(new_pair);
					}
				}
			}
		}*/
	}
}
void MoleculeManager::CalculateTriplesLevels()
{
	const int levels_num = m_non_marked_singular_points_levels.first.size();
	CV_Assert(levels_num == m_sing_pts_levels_num);
	CV_Assert(m_pairs_levels_num > 0);
	m_triples_levels.resize(m_pairs_levels_num);

	for (int level = m_pairs_levels_overlap; level < m_sing_pts_levels_num; ++level)
	{
		const int triple_level = level - m_pairs_levels_overlap;
		m_triples_levels[triple_level].clear();
		//write points to one array for simplicity
		std::vector<NonMarkedSingularPoint> combined_points_from_levels;
		std::vector<size_t> combined_types_from_levels;
		for (int curr_level = level - m_pairs_levels_overlap; curr_level <= level; ++curr_level)
		{
			combined_points_from_levels.insert(combined_points_from_levels.end(),
				m_non_marked_singular_points_levels.first[curr_level].begin(), 
				m_non_marked_singular_points_levels.first[curr_level].end());

			combined_types_from_levels.insert(combined_types_from_levels.end(), 
				m_non_marked_singular_points_levels.second[curr_level].begin(), 
				m_non_marked_singular_points_levels.second[curr_level].end());
		}

		
		for (size_t ind1 = 0; ind1 < combined_points_from_levels.size(); ++ind1)
		{
			for (size_t ind2 = ind1 + 1; ind2 < combined_points_from_levels.size(); ++ind2)
			{
				for (size_t ind3 = ind2 + 1; ind3 < combined_points_from_levels.size(); ++ind3)
				{
					MarkedSingularPoint points[3];
					points[0] = MarkedSingularPoint(combined_points_from_levels[ind1].Coord(),
						combined_types_from_levels[ind1]);
					points[1] = MarkedSingularPoint(combined_points_from_levels[ind2].Coord(),
						combined_types_from_levels[ind2]);
					points[2] = MarkedSingularPoint(combined_points_from_levels[ind3].Coord(),
						combined_types_from_levels[ind3]);

					struct PointsLess {
						bool operator()(const MarkedSingularPoint& p1, const MarkedSingularPoint& p2) {
							return p1.Property() < p2.Property();} };
					std::sort(std::begin(points), std::end(points), PointsLess());

					triangle new_triangle;
					std::get<0>(new_triangle).first = points[0];
					std::get<1>(new_triangle).first = points[1];
					std::get<2>(new_triangle).first = points[2];

					std::get<0>(new_triangle).second = cv::norm(points[0].Coord() - points[1].Coord());
					std::get<1>(new_triangle).second = cv::norm(points[2].Coord() - points[1].Coord());
					std::get<2>(new_triangle).second = cv::norm(points[0].Coord() - points[2].Coord());
					m_triples_levels[triple_level].push_back(new_triangle);
				}
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
	
	if (!m_mean_and_sigma.empty())
	{
		m_sing_pts_finder->SetMeanAndSigma(m_mean_and_sigma);
	}
	m_sing_pts_finder->CalcOnlyProps(vertices, normals, triangles, charges, wdv_radii);

	if (m_use_calculated_scale_space)
	{
		std::vector<std::vector<std::vector<double>>> scale_space;
		 ReadScaleSpaceFunctions(m_sing_pts_finder->GetScaleSpaceLevelsNum(),
			ISingularPointsFinder::FIRST_SURF_PROP, 
			static_cast<ISingularPointsFinder::SurfProperty>(ISingularPointsFinder::SURF_PROPS_NUM - 1),
			scale_space); 
		m_sing_pts_finder->SetScaleSpace(scale_space);
	}
	if (m_use_calculated_detector_function)
	{
		std::vector<std::vector<std::vector<double>>> detector_function;
		ReadScaleSpaceFunctions(m_sing_pts_finder->GetDetectorFuncLevelsNum(),
			ISingularPointsFinder::FIRST_LOG_PROP, 
			ISingularPointsFinder::LAST_LOG_PROP,
			detector_function); 
		m_sing_pts_finder->SetDetectorFunction(detector_function);
	}
	if (m_use_calculated_eig_ratio)
	{
		std::vector<std::vector<std::vector<double>>> eig_ratio_lev;
		ReadScaleSpaceFunctions(m_sing_pts_finder->GetDetectorFuncLevelsNum(),
			ISingularPointsFinder::PCA_EIG, 
			ISingularPointsFinder::PCA_EIG,
			eig_ratio_lev); 
		std::vector<std::vector<double>> eig_ratio(eig_ratio_lev.size());

		for (size_t ind = 0; ind < eig_ratio_lev.size(); ++ind)
		{
			eig_ratio[ind] = eig_ratio_lev[ind][0];
		}
		m_sing_pts_finder->SetEigRatio(eig_ratio);
	}
	m_sing_pts_finder->CalcSingPtsFromCalculatedProperties(vertices, normals, triangles, charges, wdv_radii, calc_as_average);
	m_sing_pts_finder->GetNonMarkedSingularPoints(m_non_marked_singular_points);
	m_sing_pts_finder->GetNonMarkedSingularPointsLevels(m_non_marked_singular_points_levels);
}
void MoleculeManager::CollectProperties(std::vector<std::vector<double>>& props)
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

	m_sing_pts_finder->CalcOnlyProps(vertices, normals, triangles, charges, wdv_radii);
	props.resize(ISingularPointsFinder::SURF_PROPS_NUM);
	for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
	{
		m_sing_pts_finder->AppendProp(props[curr_prop], static_cast<ISingularPointsFinder::SurfProperty>(curr_prop));
	}
}

void MoleculeManager::ReadAllSingularPoints(const int levels_num)
{
	const std::string points_file_name = m_curr_file_prefix + Extensions::SingPts();
	ReadVector(points_file_name, m_singular_points);

	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPtsTypes();
	ReadVector(non_marked_points_file_name, m_non_marked_singular_points.first);
	const std::string non_marked_points_file_name_l = m_curr_file_prefix + Extensions::NonMarkedSingPtsLabels();
	ReadVector(non_marked_points_file_name_l, m_non_marked_singular_points.second);

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

void MoleculeManager::WriteSurfaceWithTypesLevels()
{
	std::vector<std::vector<std::pair<cv::Point3d, size_t>>> vertices_with_types_lev;
	m_sing_pts_finder->GetVerticesWithTypesLevels(vertices_with_types_lev);

	for (size_t lev = 0; lev < vertices_with_types_lev.size(); ++lev)
	{
		std::stringstream s;
		s << lev;
		const std::string vert_with_types_file = m_curr_file_prefix + Extensions::SurfWithTypesLev() + s.str();
		WriteInterval(vert_with_types_file, vertices_with_types_lev[lev].begin(), vertices_with_types_lev[lev].end());
	}	
}

void MoleculeManager::WriteSurfaceWithDblProp(const ISingularPointsFinder::SurfProperty prop_type)
{
	std::vector<std::pair<cv::Point3d, double>> vertices_with_dbl_prop;
	m_sing_pts_finder->GetVerticesWithDblProp(vertices_with_dbl_prop, prop_type);
	const std::string vert_with_prop_file = m_curr_file_prefix + 
		SurfPropertyName(prop_type) + Extensions::SurfWithDblProp();
	WriteInterval(vert_with_prop_file, vertices_with_dbl_prop.begin(), vertices_with_dbl_prop.end());
}

void MoleculeManager::WriteSurfaceWithDblPropLevels(const ISingularPointsFinder::SurfProperty prop_type)
{
	std::vector<std::vector<std::pair<cv::Point3d, double>>> vertices_with_dbl_prop_lev;
	m_sing_pts_finder->GetVerticesWithDblPropLevels(vertices_with_dbl_prop_lev, prop_type);

	for (size_t lev = 0; lev < vertices_with_dbl_prop_lev.size(); ++lev)
	{
		std::stringstream s;
		s << lev;
		const std::string vert_with_prop_file = m_curr_file_prefix + SurfPropertyName(prop_type)+ 
			Extensions::SurfWithDblPropLev() + s.str();
		WriteInterval(vert_with_prop_file, vertices_with_dbl_prop_lev[lev].begin(), vertices_with_dbl_prop_lev[lev].end());
	}	
}

void MoleculeManager::ReadScaleSpaceFunctions(const size_t levels_num,
	const ISingularPointsFinder::SurfProperty first_prop, 
	const ISingularPointsFinder::SurfProperty last_prop,
	std::vector<std::vector<std::vector<double>>>& scale_space)
{
	std::vector<std::vector<std::vector<std::pair<cv::Point3d, double>>>> vertices_with_dbl_prop_lev(levels_num);
	scale_space.resize(levels_num);
	for (size_t lev = 0; lev < levels_num; ++lev)
	{
		vertices_with_dbl_prop_lev[lev].resize(last_prop - first_prop + 1);
		scale_space[lev].resize(last_prop - first_prop + 1);

		for (int curr_prop = first_prop; curr_prop <= last_prop; ++curr_prop)
		{
			auto& output = vertices_with_dbl_prop_lev[lev][curr_prop - first_prop];
			std::stringstream s;
			s << lev;
			const std::string vert_with_prop_file = m_curr_file_prefix + 
				SurfPropertyName(static_cast<ISingularPointsFinder::SurfProperty>(curr_prop))+ 
				Extensions::SurfWithDblPropLev() + s.str();
			ReadVector(vert_with_prop_file, output);

			auto& input = scale_space[lev][curr_prop - first_prop];
			input.resize(output.size());

			for (size_t ind = 0; ind < output.size(); ++ind)
			{
				input[ind] = output[ind].second;
			}
		}
	} 
}

void MoleculeManager::ReadSurfaceWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
{
	
	const std::string vert_with_types_file = m_curr_file_prefix + Extensions::SurfWithTypes();
	ReadVector(vert_with_types_file, vertices_with_types);
}

void MoleculeManager::ReadSurfaceWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types_lev)
{
	vertices_with_types_lev.resize(m_sing_pts_levels_num);
	for (size_t lev = 0; lev < m_sing_pts_levels_num; ++lev)
	{
		std::stringstream s;
		s << lev;
		const std::string vert_with_types_file = m_curr_file_prefix + Extensions::SurfWithTypesLev() + s.str();
		ReadVector(vert_with_types_file, vertices_with_types_lev[lev]);
	}
}

}
