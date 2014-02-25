#include "molecule_manager.h"

#include <algorithm>

#include "extensions.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"

namespace molecule_descriptor
{

void MoleculeManager::SetCurrFilePrefix(const std::string& file_prefix)
{
	Clear();
	m_curr_file_prefix = file_prefix;
}

void MoleculeManager::FindSingularPoints(const bool calculate)
{
	if (calculate)
	{
		CalculateSingularPoints();
		WriteSingularPoints();
		WriteSegmentedSurface();
		WriteSurfaceWithTypes();
	}
	else
	{
		ReadSingularPoints();
	}
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

void MoleculeManager::CalculatePairsWithTypes()
{
	CV_Assert(!m_dist_threshes.empty());
	m_pairs_with_types.clear();
	m_pairs_with_types.reserve(m_pairs.size());
	const size_t points_max_type = m_sing_pts_finder->GetTypesNum();
	const size_t pairs_max_type = GetPairsTypeNum();
	m_pairs_histogram.resize(pairs_max_type);
	std::fill(m_pairs_histogram.begin(), m_pairs_histogram.end(), 0);

	for (auto iter = m_pairs.begin(); iter != m_pairs.end(); ++iter)
	{
		const size_t type_1 = std::get<0>(*iter).Property();
		const size_t type_2 = std::get<1>(*iter).Property();
		const size_t dist_type = CalculateDistanceType(std::get<2>(*iter));
		const size_t pair_type = (m_dist_threshes.size() + 1) * (points_max_type * type_1 + type_2) + dist_type;
		m_pairs_with_types.push_back(std::make_tuple(
			std::get<0>(*iter), std::get<1>(*iter), pair_type));
		m_pairs_histogram[pair_type]++;
	}

}

void MoleculeManager::CalculateTriplesWithTypes()
{
	CV_Assert(!m_dist_threshes.empty());
	m_triples_with_types.clear();
	m_triples_with_types.reserve(m_triples.size());
	const size_t points_max_type = m_sing_pts_finder->GetTypesNum();
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

		const size_t dist_type_1 = triangle_types[1] = CalculateDistanceType(std::get<0>(*iter).second);
		const size_t dist_type_2 = triangle_types[3] = CalculateDistanceType(std::get<1>(*iter).second);
		const size_t dist_type_3 = triangle_types[5] = CalculateDistanceType(std::get<2>(*iter).second);

		const size_t triple_type = CalculateType(triangle_types.begin(), triangle_types.end(), 
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

void MoleculeManager::GetTriplesHistogramm(cv::Mat_<size_t>& histogram)
{
	CV_Assert(histogram.rows == 1);
	CV_Assert(histogram.cols == m_triples_histogram.size());

	std::copy(m_triples_histogram.begin(), m_triples_histogram.end(), histogram.begin());
}

size_t MoleculeManager::GetPairsTypeNum()
{
	CV_Assert(!m_dist_threshes.empty());
		const size_t points_max_type = m_sing_pts_finder->GetTypesNum();	
		return (m_dist_threshes.size() + 1) * points_max_type * points_max_type;
}

size_t MoleculeManager::GetTriplesTypeNum()
{
	CV_Assert(!m_dist_threshes.empty());
	const size_t points_max_type = m_sing_pts_finder->GetTypesNum();	
	return pow(m_dist_threshes.size() + 1, 3) * pow(points_max_type, 3);
}

void MoleculeManager::Clear()
{
	m_singular_points.clear();
	m_pairs.clear();
	m_pairs_with_types.clear();
	m_pairs_histogram.clear();
}

size_t MoleculeManager::CalculateDistanceType(const double distance)
{
	CV_Assert(!m_dist_threshes.empty());
	int type = -1;

	for (auto iter = m_dist_threshes.begin(); iter != m_dist_threshes.end() - 1; ++iter)
	{
		if (distance >= *iter && distance < *(iter + 1))
		{
			type = static_cast<int>((iter - m_dist_threshes.begin()) + 1);
			break;
		}
	}

	if (type == -1)
	{
		if (distance < m_dist_threshes.front())
		{
			type = 0;
		}
		else if (distance >= m_dist_threshes.back())
		{
			type = static_cast<int>(m_dist_threshes.size());
		}
	}

	return static_cast<size_t>(type);
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

void MoleculeManager::CalculateSingularPoints()
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

	m_sing_pts_finder->Process(vertices, normals, triangles, charges);
	m_sing_pts_finder->GetMarkedSingularPoints(m_singular_points);
	m_sing_pts_finder->GetNonMarkedSingularPoints(m_non_marked_singular_points);
	m_sing_pts_finder->GetSingularPointsHisto(m_histogram_singular_points);
}

void MoleculeManager::ReadSingularPoints()
{
	const std::string points_file_name = m_curr_file_prefix + Extensions::SingPts();
	ReadVector(points_file_name, m_singular_points);

	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPts();
	ReadVector(non_marked_points_file_name, m_non_marked_singular_points);

	const std::string hist_points_file_name = m_curr_file_prefix + Extensions::HistSingPts();
	ReadVector(hist_points_file_name, m_histogram_singular_points);
}

void MoleculeManager::WriteSingularPoints()
{
	const std::string points_file_name = m_curr_file_prefix + Extensions::SingPts();
	WriteInterval(points_file_name, m_singular_points.begin(), m_singular_points.end());

	const std::string non_marked_points_file_name = m_curr_file_prefix + Extensions::NonMarkedSingPts();
	WriteInterval(non_marked_points_file_name, m_non_marked_singular_points.begin(), m_non_marked_singular_points.end());

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
