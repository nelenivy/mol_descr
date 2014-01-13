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
	MoleculeManager() : 
		m_sing_pts_finder(CreateSingularPointsFinder()) 
	{ }
	void SetCurrFilePrefix(const std::string& file_prefix);
	void SetDistThresholds(const std::vector<double>& dist_threshes_)
	{
		m_dist_threshes = dist_threshes_;
	}
	void FindSingularPoints(const bool calculate);
	void GetSingularPoints(std::vector<MarkedSingularPoint>& sing_points) 
	{
		sing_points = m_singular_points;
	}

	void GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& sing_points) 
	{
		sing_points = m_non_marked_singular_points;
	}

	void AppendDistances(std::vector<double>& distances);
	void CalculatePairsWithTypes();
	void FindPairs(const bool calculate);
	size_t GetPairsTypeNum();
	void GetPairsHistogramm(cv::Mat_<size_t>& histogram);

	void CalculateTriplesWithTypes();
	void FindTriples(const bool calculate);
	size_t GetTriplesTypeNum();
	void GetTriplesHistogramm(cv::Mat_<size_t>& histogram);
private:
	void Clear();
	void CalculateTriples();
	void WriteTriples();
	void ReadTriples();

	void CalculatePairs();
	void WritePairs();
	void ReadPairs();
	size_t CalculateDistanceType(const double distance);
	void CalculateSingularPoints();
	void ReadSingularPoints();
	void WriteSingularPoints();
	void WriteSegmentedSurface();
	void WriteTriangles(const std::vector<cv::Point3i>& triangles);
	void WriteSurfaceWithTypes();
private:
	std::shared_ptr<ISingularPointsFinder> m_sing_pts_finder;
	typedef MarkedSingularPoint point_with_type;
	//typedef std::pair<cv::Point3d, size_t> point_with_type;
	typedef std::tuple<point_with_type, point_with_type, double> pair_with_distance;
	typedef std::tuple<point_with_type, point_with_type, size_t> pair_with_distance_type;

	typedef std::pair<point_with_type, double> triangle_side;
	typedef std::pair<point_with_type, size_t> triangle_side_with_type;
	typedef std::tuple<triangle_side, triangle_side, triangle_side> triangle;
	typedef std::tuple<triangle_side_with_type, triangle_side_with_type, triangle_side_with_type> triangle_with_type;

	std::vector<point_with_type> m_singular_points;
	std::vector<NonMarkedSingularPoint> m_non_marked_singular_points;
	std::vector<pair_with_distance> m_pairs;
	std::vector<pair_with_distance_type> m_pairs_with_types;
	std::vector<size_t> m_pairs_histogram;

	std::vector<triangle> m_triples;
	std::vector<triangle_with_type> m_triples_with_types;
	std::vector<size_t> m_triples_histogram;

	std::string m_curr_file_prefix;
	std::vector<double> m_dist_threshes;
};

}