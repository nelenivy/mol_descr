#pragma once
#include <vector>
#include <algorithm>

#include "opencv2/core/core.hpp"

#include "CommonUtilities/attributes_container.h"
#include "graph_structures.h"
#include "segmentation_types.h"
#include "CommonUtilities/common_functions.h"
#include "CommonUtilities/priority_heap.h"

namespace molecule_descriptor
{

template <class NodeType>
class KMeansSegmentator
{//graph as 
public:
	void Segment(const vector<GraphNode<NodeType>>& graph,	const size_t expected_clust_num, vector<size_t>& segmented_graph, size_t& actual_clust_num);
private:
	static float kMaxDist() { return 1000000.0f; }
	static float kNoConnectionMark() { return -1.0f; }
	static const int kMaxIterations = 1000;
	void InitForNewData(const vector<GraphNode<NodeType>>& graph);
	void CalculateDistanceMatrix(const vector<GraphNode<NodeType>>& graph);
	void KMeansImpl(const vector<GraphNode<NodeType>>& graph, const size_t clust_num);
	cv::Mat_<float> m_distance_mat;
	vector<bool> m_visited;
	vector<bool> m_curr_subgraph_mask;
	vector<size_t>  m_local_segmented_graph;
};

template <class NodeType>
void KMeansSegmentator<NodeType>::Segment(const vector<GraphNode<NodeType>>& graph, 
					const size_t expected_clust_num, vector<size_t>& segmented_graph, size_t& actual_clust_num)
{
	CV_Assert(expected_clust_num);
	InitForNewData(graph);
	segmented_graph.resize(graph.size());
	std::fill(segmented_graph.begin(), segmented_graph.end(), 0);
	actual_clust_num = 0;

	for (size_t ind = 0; ind < graph.size(); ind++)
	{//for each connectivity component run k-means
		if (m_visited[ind])
		{
			continue;
		}

		std::fill(m_curr_subgraph_mask.begin(), m_curr_subgraph_mask.end(), false);
		std::fill(m_local_segmented_graph.begin(), m_local_segmented_graph.end(), 0);
		size_t curr_vertices_num = 0;
		//count vertices in current connectivity component
		//and mark nodes

		for (size_t ind_node = 0; ind_node < graph.size(); ind_node++)
		{
			if (m_distance_mat(static_cast<int>(ind), static_cast<int>(ind_node)) >= 0.0f)
			{
				m_visited[ind_node] = true;
				m_curr_subgraph_mask[ind_node] = true;
				curr_vertices_num++;
			}
		}

		const size_t clust_num = static_cast<size_t>(
			Round(expected_clust_num * static_cast<double>(curr_vertices_num) / static_cast<double>(graph.size())));
		KMeansImpl(graph, clust_num);

		for (size_t ind_node = 0; ind_node < graph.size(); ind_node++)
		{//write from local segments vector to the whole 
			if (m_local_segmented_graph[ind_node] > 0)
			{
				segmented_graph[ind_node] = m_local_segmented_graph[ind_node] + actual_clust_num;
			}
		}

		actual_clust_num += clust_num;
	}
}

template <typename NodeType>
void KMeansSegmentator<NodeType>::InitForNewData(const vector<GraphNode<NodeType>>& graph)
{
	CalculateDistanceMatrix(graph);
	m_visited.resize(graph.size(), false);
	m_curr_subgraph_mask.resize(graph.size(), false);
	m_local_segmented_graph.resize(graph.size(), 0);
}

INT_PROP(IndexProp);

template <typename NodeType>
void KMeansSegmentator<NodeType>::CalculateDistanceMatrix(const vector<GraphNode<NodeType>>& graph)
{
	const int graph_size = static_cast<int>(graph.size());
	m_distance_mat.create(graph_size, graph_size);
	cv::Mat_<float> adjacency_mat = cv::Mat::zeros(graph_size, graph_size, CV_32FC1);
	//add indexes property
	for (int ind = 0; ind < graph_size; ind++)
	{
		graph[ind].attr.Add<IndexProp>();
		graph[ind].attr.Get<IndexProp>() = ind;
	}
	//fill adjacency matrix
	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{
		const int curr_ind = node_iter->attr.Get<IndexProp>();

		for (auto neighb_iter = node_iter->neighbours.begin(); neighb_iter != node_iter->neighbours.end(); ++neighb_iter)
		{
			const int neighb_ind = (*neighb_iter)->attr.Get<IndexProp>();
			adjacency_mat(curr_ind, neighb_ind) = adjacency_mat(neighb_ind, curr_ind)
				= static_cast<float>(cv::norm(node_iter->element->Center() - (*neighb_iter)->element->Center()));
		}
	}

	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{//for each node find shortest distances to other nodes
		auto distance_queue = MakeRisingProrityQueue(graph.begin(), graph.end(), kMaxDist());
		const int base_ind = node_iter->attr.Get<IndexProp>();
		distance_queue.ChangeKey(base_ind, 0.0f, eThrow);

		while (distance_queue.Size() > 0)
		{
			//extract top
			const GraphNode<NodeType>& curr_top_elem = *(distance_queue.Top().Pointer());
			float curr_min_dist =  distance_queue.Top().Key();
			distance_queue.Pop();
			const int curr_ind = curr_top_elem.attr.Get<IndexProp>();

			for (auto neighb_iter = curr_top_elem.neighbours.begin(); neighb_iter != curr_top_elem.neighbours.end(); ++neighb_iter)
			{//update distances for neighbours
				const int neighb_ind = (*neighb_iter)->attr.Get<IndexProp>();
				distance_queue.ChangeKey(static_cast<size_t>(neighb_ind), curr_min_dist + adjacency_mat(curr_ind, neighb_ind), eNotThrow);
			}
			//write distance to matrix
			if (curr_min_dist >= kMaxDist())
			{
				curr_min_dist = kNoConnectionMark();
			}

			m_distance_mat(base_ind, curr_ind) = m_distance_mat(curr_ind, base_ind) = curr_min_dist;
		}
	}
}


template <typename NodeType>
void KMeansSegmentator<NodeType>::KMeansImpl(const vector<GraphNode<NodeType>>& graph, const size_t clust_num)
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
	for (size_t ind = 0, ind_in_subgraph = 0, curr_clust = 1; ind < graph.size() && curr_clust <= clust_num; ind++)
	{
		if (!m_curr_subgraph_mask[ind])
		{
			continue;
		}

		if (ind_in_subgraph % mean_elem_num_in_segm == 0)
		{
			centers_indexes[curr_clust] = ind;
			curr_clust++;
		}
		ind_in_subgraph++;
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
		for (size_t ind = 0; ind < graph.size(); ++ind)
		{
			if (!m_curr_subgraph_mask[ind])
			{
				continue;
			}

			float min_dist = kMaxDist();

			for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
			{
				const float curr_dist = m_distance_mat(static_cast<int>(ind), static_cast<int>(centers_indexes[curr_clust]));
				if (curr_dist < min_dist)
				{
					min_dist = curr_dist;
					m_local_segmented_graph[ind] = curr_clust;
				}
			}
		}

		//calculate new centers as mean of the cluster coordinates
		std::fill(actual_centers.begin(), actual_centers.end(), cv::Point3d(0.0, 0.0, 0.0));
		std::fill(curr_clust_sizes.begin(), curr_clust_sizes.end(), 0);

		for (size_t ind = 0; ind < graph.size(); ++ind)
		{
			const size_t curr_clust_num = m_local_segmented_graph[ind];

			if (curr_clust_num > 0)
			{
				actual_centers[curr_clust_num] += graph[ind].element->Center();
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
		std::fill(min_dist_to_center.begin(), min_dist_to_center.end(), make_pair(static_cast<float>(kMaxDist()), 0));

		for (size_t ind = 0; ind < graph.size(); ++ind)
		{
			const size_t curr_clust_num = m_local_segmented_graph[ind];

			if (curr_clust_num > 0)
			{
				const double curr_dist_to_center = cv::norm(graph[ind].element->Center() - actual_centers[curr_clust_num]);

				if (curr_dist_to_center < min_dist_to_center[curr_clust_num].first)
				{
					min_dist_to_center[curr_clust_num] = make_pair(curr_dist_to_center, ind);
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