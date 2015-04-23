#pragma once
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <array>
#include <memory>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "SingularPoints/i_singular_points_finder.h"

namespace molecule_descriptor
{
//maybe temporal class
class MoleculeManager
{
public:	
	MoleculeManager() : m_sing_pts_levels_num(0)
	{ }
	void Init(int argc, char** argv)
	{
		m_sing_pts_finder->InitParams(argc, argv);
	}
	void SetSingPtsAlgorithm(const SingularPointsAlgorithm alg)
	{
		m_sing_pts_finder = CreateSingularPointsFinder(alg);
	}
	void SetLevelsNum(const int levels_num)
	{
		m_sing_pts_levels_num = levels_num;
	}
	void SetPairsLevelsOverlap(const int levels_overlap)
	{
		m_pairs_levels_overlap = levels_overlap;
	}
	void SetLevelsScales(const std::vector<double>& levels_scales)
	{
		m_levels_scales = levels_scales;
	}
	void SetCurrFilePrefix(const std::string& file_prefix);
	void SetDistThresholds(const std::vector<double>& dist_threshes_)
	{
		m_dist_threshes = dist_threshes_;
	}
	void SetDistLevelsThresholds(const std::vector<std::vector<double>>& dist_threshes_level, const std::vector<double>& dist_high_thresh)
	{
		CV_Assert(m_pairs_levels_num == dist_threshes_level.size());
		CV_Assert(m_pairs_levels_num == dist_high_thresh.size());
		m_dist_threshes_levels = dist_threshes_level;
		m_dist_high_threshes = dist_high_thresh;
	}
	void SetChargesThresholds(const std::vector<double>& charges_threshes_)
	{
		m_charge_threshes = charges_threshes_;
	}
	void SetLennardJonesThresholds(const std::vector<double>& lennard_jones_threshes_)
	{
		m_lennard_jones_threshes = lennard_jones_threshes_;
	}
	void SetAreaThresholds(const std::vector<double>& area_threshes_)
	{
		m_area_threshes = area_threshes_;
	}
	void SetMeansAndSigma(std::vector<std::vector<double>> mean_and_sigma)
	{
		m_mean_and_sigma = mean_and_sigma;
	}
	void CollectProperties(std::vector<std::vector<double>>& props);
	void FindSingularPoints(const bool calculate, const bool calc_as_average);
	void ReadAllSingularPoints(const int levels_num);
	void ReadSurfaceWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	void ReadSurfaceWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types_lev);

	size_t GetPointsTypeNum();
	void ClassifySingularPoints();
	void GetSingularPoints(std::vector<MarkedSingularPoint>& sing_points) 
	{
		sing_points = m_singular_points;
	}

	void GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& sing_points) 
	{
		sing_points = m_non_marked_singular_points.first;
	}

	void AppendDistances(std::vector<double>& distances);
	void AppendDistancesLevels(std::vector<std::vector<double>>& distances);
	void AppendCharges(std::vector<double>& charges, const bool levels = false);
	void AppendLennardJones(std::vector<double>& lennard_jones, const bool levels = false);
	void AppendAreas(std::vector<double>& areas, const bool levels = false);
	void CalculatePropertiesTypes();
	void CalculatePairsWithTypes();
	void CalculatePairsWithTypesLevels();
	void CalculatePairsWithTypesLevelsAllMesh();
	void CalculatePairsWithTypesLevelsAllMeshSmoothedCurv();
	void FindPairs(const bool calculate);
	void FindPairsLevels(const bool calculate, const bool write);
	size_t GetPairsTypeNum();
	size_t GetPairsTypeNumLevels();
	void GetPairsHistogramm(cv::Mat_<size_t>& histogram);
	void GetPairsHistogrammLevels(cv::Mat_<double>& histogram);
	void GetNonMarkedPairs(std::vector<SingularPointsPair<PropertiesSet>>& non_marked_pairs)
	{
		non_marked_pairs = m_sing_pts_pairs_with_props_and_types.first;
	}
	void GetPairsTypes(std::vector<size_t>& pairs_types)
	{
		pairs_types = m_sing_pts_pairs_with_props_and_types.second;
	}
	void CalculateTriplesWithTypes();
	void FindTriples(const bool calculate);
	size_t GetTriplesTypeNum();
	void GetTriplesHistogramm(cv::Mat_<size_t>& histogram);

	void FindTriplesLevels(const bool calculate, const bool write);
	void CalculateTriplesWithTypesLevels();
	size_t GetTriplesTypeNumLevels();
	void GetTriplesHistogrammLevels(cv::Mat_<float>& histogram);
private:
	int CalculateSingularPointsTypes(PropertiesSet& prop);

	void Clear();
	void CalculateTriples();
	void WriteTriples();
	void ReadTriples();

	void CalculatePairs();
	void WritePairs();
	void ReadPairs();

	void CalculatePairsLevels();
	void WritePairsLevels();
	void ReadPairsLevels(const int levels_num);

	void CalculateTriplesLevels();

	void CalculateSingularPoints(const bool calc_as_average);
	void WriteSingularPoints();
	void WriteSegmentedSurface();
	void WriteTriangles(const std::vector<cv::Point3i>& triangles);
	void WriteSurfaceWithTypes();
	void WriteSurfaceWithTypesLevels();
	void ReadSingularPointsLevelsTypes(const int levels_num);
	void ReadSingularPointsLevelsLabels(const int levels_num);
private:
	std::shared_ptr<ISingularPointsFinder> m_sing_pts_finder;
	static const size_t kHistSize = 9;
	typedef MarkedSingularPoint point_with_type;
	//typedef std::pair<cv::Point3d, size_t> point_with_type;
	typedef std::tuple<point_with_type, point_with_type, double> pair_with_distance;
	typedef std::tuple<point_with_type, point_with_type, size_t> pair_with_distance_type;

	typedef std::pair<point_with_type, double> triangle_side;
	typedef std::pair<point_with_type, size_t> triangle_side_with_type;
	typedef std::tuple<triangle_side, triangle_side, triangle_side> triangle;
	typedef std::tuple<triangle_side_with_type, triangle_side_with_type, triangle_side_with_type> triangle_with_type;

	std::vector<point_with_type> m_singular_points;
	std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>> m_non_marked_singular_points;
	std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, std::vector<std::vector<size_t>>> m_non_marked_singular_points_levels;
	std::pair<std::vector<SingularPointsPair<PropertiesSet>>, std::vector<size_t>> m_sing_pts_pairs_with_props_and_types;
	std::pair<std::vector<std::vector<SingularPointsPair<PropertiesSet>>>, 
		std::vector<std::vector<size_t>>> m_sing_pts_pairs_with_props_and_types_levels;

	std::vector<pair_with_distance> m_pairs;	
	std::vector<pair_with_distance_type> m_pairs_with_types;
	std::vector<size_t> m_types;
	std::vector<size_t> m_pairs_histogram;

	std::vector<std::vector<pair_with_distance>> m_pairs_levels;
	std::vector<std::vector<pair_with_distance_type>> m_pairs_with_types_levels;
	std::vector<std::vector<size_t>> m_types_levels;
	std::vector<double> m_pairs_histogram_levels;
	std::vector<double> m_levels_scales;

	std::vector<triangle> m_triples;
	std::vector<triangle_with_type> m_triples_with_types;
	std::vector<size_t> m_triples_histogram;

	std::vector<std::vector<triangle>> m_triples_levels;
	std::vector<std::vector<triangle_with_type>> m_triples_with_types_levels;
	std::vector<size_t> m_triples_histogram_levels;

	std::string m_curr_file_prefix;
	std::vector<double> m_dist_threshes;
	std::vector<std::vector<double>> m_dist_threshes_levels;
	std::vector<double> m_dist_high_threshes;
	std::vector<double> m_charge_threshes;
	std::vector<double> m_lennard_jones_threshes;
	std::vector<double> m_area_threshes;
	std::vector<std::vector<double>> m_mean_and_sigma;
	int m_sing_pts_levels_num;
	int m_pairs_levels_num;
	int m_pairs_levels_overlap;
	std::vector<pair_with_distance> m_pairs_vertices;
};

}