#pragma once

#include "boost/graph/properties.hpp"
#include "boost/graph/subgraph.hpp"
#include "boost/graph/graph_traits.hpp"
#include "GraphLib/array_property_map.h"
#include "k_means_segmentator.h"
#include "connected_components_segmentator.h"

namespace molecule_descriptor
{
template <typename GraphType, typename SubgraphMaskMap>
void CreateSubgraph(boost::subgraph<GraphType>& graph, const SubgraphMaskMap& subgraph_mask, 
					boost::subgraph<GraphType>& subgraph);
template <typename GraphType, typename SegmentsMap>
void SplitSegment(const boost::subgraph<GraphType>& graph, const size_t segment_to_split, const size_t parts_num, 
				  SegmentsMap& segm_map, size_t& segments_num);
template <typename GraphType, typename SegmentsMap, typename CentersGraph>
void CalculateSegmentGraphAndCenters(const GraphType& graph, const SegmentsMap& segments_map, 
									 std::vector<typename boost::graph_traits<GraphType>::vertex_descriptor>& centers, 
									 CentersGraph& centers_graph);

template <typename GraphType, typename SubgraphMaskMap>
void CreateSubgraph(boost::subgraph<GraphType>& graph, const SubgraphMaskMap& subgraph_mask, 
					boost::subgraph<GraphType>& subgraph)
{
	subgraph = graph.create_subgraph();
	typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_iterator curr_vert, end_vert;
	for (boost::tie(curr_vert, end_vert) = vertices(graph); curr_vert != end_vert; ++curr_vert)
	{
		if (static_cast<bool>(subgraph_mask[*curr_vert]) == true)
		{
			add_vertex(*curr_vert, subgraph);
		}
	}	
}

template <typename GraphType, typename SegmentsMap>
void SplitSegment(const boost::subgraph<GraphType>& graph, const size_t segment_to_split, const size_t parts_num, 
				  SegmentsMap& segm_map, size_t& segments_num)
{
	ContPropMap<boost::subgraph<GraphType>, std::vector<bool>, VERTEX> subgraph_mask(graph);

	typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_iterator curr_vert, end_vert;
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{//mark necessary nodes
		subgraph_mask[*curr_vert] = segm_map[*curr_vert] == segment_to_split ? true : false;
	}

	//segment using k_means
	boost::subgraph<GraphType> graph_copy = graph;
	boost::subgraph<GraphType> subgraph;
	CreateSubgraph(graph_copy, subgraph_mask, subgraph);

	if (num_vertices(subgraph) == 0)
	{
		return;
	}

	KMeansSegmentator<boost::subgraph<GraphType>> kmeans_segmentator;
	ContPropMap<boost::subgraph<GraphType>, std::vector<size_t>, VERTEX> segmented_subgraph(subgraph);
	size_t new_segments_num = 0;
	kmeans_segmentator.Segment(subgraph, get(boost::vertex_index, subgraph), parts_num, segmented_subgraph, new_segments_num);
	//write segments back
	//renum segments
	std::vector <size_t> segments_map(new_segments_num + 1, 0);
	segments_map[1] = segment_to_split;
	for (size_t ind = 2; ind < segments_map.size(); ++ind)
	{
		segments_map[ind] = segments_num + ind - 1;
	}
	//write new value of segments number
	segments_num += new_segments_num - 1;
	for (auto curr_vert = vertices(subgraph).first, end_vert = vertices(subgraph).second; curr_vert != end_vert; ++curr_vert)
	{
		const size_t curr_segm_num = segments_map[segmented_subgraph[*curr_vert]];
		const auto parent_vert = subgraph.local_to_global(*curr_vert);
		segm_map[parent_vert] = curr_segm_num;		
	}
}

template <typename GraphType, typename SegmentsMap, typename CentersGraph>
void CalculateSegmentGraphAndCenters(const GraphType& graph, const SegmentsMap& segments_map, 
	std::vector<typename boost::graph_traits<GraphType>::vertex_descriptor>& centers, 
	CentersGraph& centers_graph)
{
	//find clusters number
	size_t clust_num = 0;

	for (auto node_iter = vertices(graph).first, graph_end = vertices(graph).second; node_iter != graph_end; ++node_iter)
	{
		const size_t curr_segm_num = segments_map[*node_iter];

		if (curr_segm_num > clust_num)
		{
			clust_num = curr_segm_num;
		}
	}

	//calculate new centers as mean of the cluster coordinates
	std::vector<cv::Point3d> actual_centers(clust_num + 1, cv::Point3d(0.0, 0.0, 0.0));//clust + 1 because clusters numbers go from 1
	std::vector<size_t> curr_clust_sizes(clust_num + 1, 0);

	for (auto node_iter = vertices(graph).first, graph_end = vertices(graph).second; node_iter != graph_end; ++node_iter)
	{
		const size_t curr_clust_num = segments_map[*node_iter];

		if (curr_clust_num > 0)
		{
			actual_centers[curr_clust_num] += get(boost::vertex_info_3d, graph, *node_iter).Center();
			curr_clust_sizes[curr_clust_num]++;
		}
	}

	for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
	{
		if (curr_clust_sizes[curr_clust])
		{
			actual_centers[curr_clust] *= 1.0 / curr_clust_sizes[curr_clust];
		}
	}

	//find the nearest point from the set to the actual center
	const double kMaxDist = 1000000.0;
	typedef typename boost::graph_traits<GraphType>::vertex_descriptor vertex_descriptor;
	std::vector<std::pair<double, vertex_descriptor>> min_dist_to_center(clust_num + 1, make_pair(kMaxDist, vertex_descriptor()));
	//minimum distance on the first place, index on the second

	for (auto node_iter = vertices(graph).first, graph_end = vertices(graph).second; node_iter != graph_end; ++node_iter)
	{
		const size_t curr_clust_num = segments_map[*node_iter];

		if (curr_clust_num > 0)
		{
			const double curr_dist_to_center = cv::norm(get(boost::vertex_info_3d, graph, *node_iter).Center() - actual_centers[curr_clust_num]);

			if (curr_dist_to_center < min_dist_to_center[curr_clust_num].first)
			{
				min_dist_to_center[curr_clust_num] = make_pair(curr_dist_to_center, *node_iter);
			}
		}
	}

	//write found nodes to the output vector
	centers.resize(clust_num);
	std::vector<typename boost::graph_traits<CentersGraph>::vertex_descriptor> centers_nums(clust_num + 1);
	centers_graph = CentersGraph();
	size_t non_zero_cluster = 0;
	for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
	{
		if (curr_clust_sizes[curr_clust])
		{
			centers[non_zero_cluster] = min_dist_to_center[curr_clust].second;
			const auto added = add_vertex(centers_graph);
			centers_nums[curr_clust] = added;
			put(boost::vertex_info_3d, centers_graph, added, 
				get(boost::vertex_info_3d, graph, centers[non_zero_cluster]));
			put(boost::vertex_parent, centers_graph, added, centers[non_zero_cluster]);		
			++non_zero_cluster;
		}
	}

	centers.resize(non_zero_cluster);
	//calculate segments graph
	for (auto node_iter = vertices(graph).first, graph_end = vertices(graph).second; node_iter != graph_end; ++node_iter)
	{
		const size_t curr_clust_num = segments_map[*node_iter];

		if (curr_clust_num == 0)
		{
			continue;
		}

		for (auto neighb_iter = adjacent_vertices(*node_iter, graph).first, 
			neighb_end = adjacent_vertices(*node_iter, graph).second; 
			neighb_iter != neighb_end; ++neighb_iter)
		{
			const size_t neighb_clust_num = segments_map[*neighb_iter];
			if (neighb_clust_num == 0)
			{
				continue;
			}
			if (neighb_clust_num != curr_clust_num)
			{
				if (!edge(centers_nums[curr_clust_num], centers_nums[neighb_clust_num], centers_graph).second)
				{
					add_edge(centers_nums[curr_clust_num], centers_nums[neighb_clust_num], centers_graph);
				}
			}
		}
	}
}

}