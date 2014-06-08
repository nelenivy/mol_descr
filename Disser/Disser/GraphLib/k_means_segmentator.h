#pragma once
#include <vector>
#include <algorithm>

#include "opencv2/core/core.hpp"

#include "GraphLib/dijkstra_dist_matr.h"
#include "GraphLib/array_property_map.h"
#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

template <class GraphType>
class KMeansSegmentator
{//graph as 
public:
	template <typename SegmentsMap, typename IndexMap>
	void Segment(const GraphType& graph, const IndexMap& index_map, const size_t expected_clust_num, 
		SegmentsMap& segmented_graph, size_t& actual_clust_num);
private:
	static const int kMaxIterations = 3000;
	void InitForNewData(const GraphType& graph);
	template <typename IndexMap>
	void KMeansImpl(const GraphType& graph, const IndexMap& index_map, const size_t clust_num);

	DijkstraDistMapCalculator<GraphType> m_dist_mat_calculator;
	cv::Mat_<float> m_distance_mat;
	ContPropMap<GraphType, vector<bool>, VERTEX> m_visited;
	ContPropMap<GraphType, vector<bool>, VERTEX> m_curr_subgraph_mask;
	ContPropMap<GraphType, vector<size_t>, VERTEX>  m_local_segmented_graph;
};

template <class GraphType>
template <typename SegmentsMap, typename IndexMap>
void KMeansSegmentator<GraphType>::Segment(const GraphType& graph, const IndexMap& index_map,
					const size_t expected_clust_num, SegmentsMap& segmented_graph, size_t& actual_clust_num)
{
	if (expected_clust_num == 0)
	{
		return;
	}
	InitForNewData(graph);
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
		curr_vert != end_vert; ++curr_vert)
	{
		segmented_graph[*curr_vert] = 0;
	}

	actual_clust_num = 0;

	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
		curr_vert != end_vert; ++curr_vert)
	{//for each connectivity component run k-means
		if (m_visited[*curr_vert])
		{
			continue;
		}

		std::fill(m_curr_subgraph_mask.begin(), m_curr_subgraph_mask.end(), false);
		std::fill(m_local_segmented_graph.begin(), m_local_segmented_graph.end(), 0);
		size_t curr_vertices_num = 0;
		//count vertices in current connectivity component
		//and mark nodes
		for (auto curr_vert_node = vertices(graph).first, end_vert_node = vertices(graph).second; 
			curr_vert_node != end_vert_node; ++curr_vert_node)
		{
			if (m_distance_mat(index_map[*curr_vert], index_map[*curr_vert_node]) >= 0.0f)
			{
				m_visited[*curr_vert_node] = true;
				m_curr_subgraph_mask[*curr_vert_node] = true;
				curr_vertices_num++;
			}
		}

		const size_t clust_num = static_cast<size_t>(
			Round(expected_clust_num * static_cast<double>(curr_vertices_num) / static_cast<double>(num_vertices(graph))));
		KMeansImpl(graph, get(boost::vertex_index, graph), clust_num);

		for (auto curr_vert_node = vertices(graph).first, end_vert_node = vertices(graph).second; 
			curr_vert_node != end_vert_node; ++curr_vert_node)
		{//write from local segments vector to the whole 
			if (m_local_segmented_graph[*curr_vert_node] > 0)
			{
				segmented_graph[*curr_vert_node] = m_local_segmented_graph[*curr_vert_node] + actual_clust_num;
			}
		}

		actual_clust_num += clust_num;
	}
}

template <typename GraphType>
void KMeansSegmentator<GraphType>::InitForNewData(const GraphType& graph)
{
	m_dist_mat_calculator.Calc(graph, m_distance_mat);
	m_visited.SetGraph(graph);	
	m_curr_subgraph_mask.SetGraph(graph);
	m_local_segmented_graph.SetGraph(graph);

	std::fill(m_local_segmented_graph.begin(), m_local_segmented_graph.end(), 0);
	std::fill(m_visited.begin(), m_visited.end(), false);
	std::fill(m_curr_subgraph_mask.begin(), m_curr_subgraph_mask.end(), false);
}

template <typename GraphType>
template <typename IndexMap>
void KMeansSegmentator<GraphType>::KMeansImpl(const GraphType& graph, const IndexMap& index_map, const size_t clust_num)
{
	if (clust_num == 0)
	{
		return;
	}

	std::fill(m_local_segmented_graph.begin(), m_local_segmented_graph.end(), 0);

	const size_t all_elems_num = std::count(m_curr_subgraph_mask.begin(), m_curr_subgraph_mask.end(), true);
	const size_t mean_elem_num_in_segm = all_elems_num / clust_num;

#pragma region init seeds
	vector<size_t> centers_indexes(clust_num + 1, 0);//clust + 1 because clusters numbers go from 1
	//indexes because we obtain distance through distance matrix
	{
		int ind_in_subgraph = 0, curr_clust = 1;
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
			curr_vert != end_vert && curr_clust <= clust_num; ++curr_vert)
		{
			if (!m_curr_subgraph_mask[*curr_vert])
			{
				continue;
			}

			if (ind_in_subgraph % mean_elem_num_in_segm == 0)
			{
				centers_indexes[curr_clust] = index_map[*curr_vert];
				curr_clust++;
			}
			ind_in_subgraph++;
		}
	}
#pragma endregion init seeds
	
#pragma region main cycle
	size_t changes_number = 1;
	size_t curr_iteration = 0;
	vector<cv::Point3d> actual_centers(clust_num + 1);//clust + 1 because clusters numbers go from 1
	vector<size_t> curr_clust_sizes(clust_num + 1, 0);
	vector<std::pair<double, size_t>> min_dist_to_center(clust_num + 1);//minimum distance on the first place, index on the second

	while (changes_number > 0 && curr_iteration < kMaxIterations)
	{
		curr_iteration++;
		//assign elements to clusters
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
			curr_vert != end_vert; ++curr_vert)
		{
			if (!m_curr_subgraph_mask[*curr_vert])
			{
				continue;
			}

			float min_dist = DijkstraDistMapCalculator<GraphType>::kMaxDist();

			for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
			{
				const float curr_dist = m_distance_mat(index_map[*curr_vert], static_cast<int>(centers_indexes[curr_clust]));
				if (curr_dist < min_dist)
				{
					min_dist = curr_dist;
					m_local_segmented_graph[*curr_vert] = curr_clust;
				}
			}
		}

		//calculate new centers as mean of the cluster coordinates
		std::fill(actual_centers.begin(), actual_centers.end(), cv::Point3d(0.0, 0.0, 0.0));
		std::fill(curr_clust_sizes.begin(), curr_clust_sizes.end(), 0);

		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
			curr_vert != end_vert; ++curr_vert)
		{
			const size_t curr_clust_num = m_local_segmented_graph[*curr_vert];

			if (curr_clust_num > 0)
			{
				actual_centers[curr_clust_num] += get(boost::vertex_info_3d, graph, *curr_vert).Center();
				curr_clust_sizes[curr_clust_num]++;
			}
		}

		for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
		{
			actual_centers[curr_clust].x /= curr_clust_sizes[curr_clust];
			actual_centers[curr_clust].y /= curr_clust_sizes[curr_clust];
			actual_centers[curr_clust].z /= curr_clust_sizes[curr_clust];
		}

		//find the nearest point from the set to the actual center
		std::fill(min_dist_to_center.begin(), min_dist_to_center.end(), 
			make_pair(static_cast<float>(DijkstraDistMapCalculator<GraphType>::kMaxDist()), 0));

		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; 
			curr_vert != end_vert; ++curr_vert)
		{
			const size_t curr_clust_num = m_local_segmented_graph[*curr_vert];

			if (curr_clust_num > 0)
			{
				const double curr_dist_to_center = cv::norm(get(boost::vertex_info_3d, graph, *curr_vert).Center() - actual_centers[curr_clust_num]);

				if (curr_dist_to_center < min_dist_to_center[curr_clust_num].first)
				{
					min_dist_to_center[curr_clust_num] = make_pair(curr_dist_to_center, index_map[*curr_vert]);
				}
			}
		}

		//find the number of changes in centers coordinates
		changes_number = 0;

		for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
		{
			if (centers_indexes[curr_clust] != min_dist_to_center[curr_clust].second)
			{
				changes_number++;
			}
		}
	}
}
#pragma endregion main cycle
}