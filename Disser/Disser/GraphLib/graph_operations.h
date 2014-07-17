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
void CreateSubgraphFromMask(boost::subgraph<GraphType>& graph, const SubgraphMaskMap& subgraph_mask, 
	const typename boost::property_traits<SubgraphMaskMap>::value_type subgraph_mark,
	boost::subgraph<GraphType>& subgraph);
template <typename GraphType, typename SegmentsMap>
void SplitSegment(const boost::subgraph<GraphType>& graph, const size_t segment_to_split, const size_t parts_num, 
				  SegmentsMap& segm_map, size_t& segments_num);
template <typename GraphType, typename SegmentsMap, typename CentersGraph>
void CalculateSegmentGraphAndCenters(const boost::subgraph<GraphType>& graph, const SegmentsMap& segments_map, 
									 std::vector<typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_descriptor>& centers, 
									 CentersGraph& centers_graph);

template <typename GraphType, typename SubgraphMaskMap>
void CreateSubgraphFromMask(boost::subgraph<GraphType>& graph, const SubgraphMaskMap& subgraph_mask, 
							const typename boost::property_traits<SubgraphMaskMap>::value_type subgraph_mark,
					boost::subgraph<GraphType>& subgraph)
{
	subgraph = graph.create_subgraph();
	typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_iterator curr_vert, end_vert;
	for (boost::tie(curr_vert, end_vert) = vertices(graph); curr_vert != end_vert; ++curr_vert)
	{
		if (subgraph_mask[*curr_vert] == subgraph_mark)
		{
			add_vertex(*curr_vert, subgraph);
		}
	}	
}

template <typename GraphType, typename SegmentsMap, typename CoordMap>
std::pair<typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_descriptor, bool>
	CalculateSegmentCenter(boost::subgraph<GraphType>& graph, const SegmentsMap& segm_map, 
							const typename boost::property_traits<SegmentsMap>::value_type curr_segm,
							const CoordMap& coord_map)
{
	typedef typename boost::graph_traits<boost::subgraph<GraphType>>::vertex_descriptor vertex_descriptor;
	boost::subgraph<GraphType> subgraph;
	CreateSubgraphFromMask(graph, segm_map, curr_segm, subgraph);

	if (num_vertices(subgraph) == 0)
	{
		return make_pair(vertex_descriptor(), false);
	}
	ContPropMap<boost::subgraph<GraphType>, std::vector<double>, EDGE> edge_distances(subgraph);
	std::fill(edge_distances.begin(), edge_distances.end(), 0.0);

	for (auto curr_edge = edges(subgraph).first, edge_end = edges(subgraph).second; curr_edge != edge_end; ++curr_edge)
	{
		const vertex_descriptor v1 = source(*curr_edge, subgraph);
		const vertex_descriptor v2 = target(*curr_edge, subgraph);
		const double curr_dist = Distance(coord_map[subgraph.local_to_global(v1)], coord_map[subgraph.local_to_global(v2)]);
		edge_distances[*curr_edge] = curr_dist;
	}
	cv::Mat_<float> dist_mat;
	CalcDistMatr(subgraph, get(boost::vertex_index, subgraph), edge_distances, 
		std::numeric_limits<float>::max(), 0.0, dist_mat);

	int min_pos = 0;
	double min_dist = std::numeric_limits<double>::max();
	for (int y = 0; y < dist_mat.rows; ++y)
	{
		double sum_dist = 0;
		for (int x = 0; x < dist_mat.cols;++x)
		{
			sum_dist += dist_mat(y, x);
		}
		if (sum_dist < min_dist)
		{
			min_pos = y;
			min_dist = sum_dist;
		}
	}

	const vertex_descriptor local_cent = vertex(min_pos, subgraph);
	return make_pair(subgraph.local_to_global(local_cent), true);
}

template <typename GraphType, typename SegmentsMap>
void SplitSegment(const boost::subgraph<GraphType>& graph, const size_t segment_to_split, const size_t parts_num, 
				  SegmentsMap& segm_map, size_t& segments_num)
{
	CV_Assert(parts_num > 0);
	ContPropMap<boost::subgraph<GraphType>, std::vector<bool>, VERTEX> subgraph_mask(graph);

	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{//mark necessary nodes
		subgraph_mask[*curr_vert] = segm_map[*curr_vert] == segment_to_split ? true : false;
	}

	//segment using k_means
	boost::subgraph<GraphType> graph_copy = graph;
	boost::subgraph<GraphType> subgraph;
	CreateSubgraphFromMask(graph_copy, subgraph_mask, true, subgraph);

	if (num_vertices(subgraph) == 0)
	{
		return;
	}

	KMeansSegmentator<boost::subgraph<GraphType>> kmeans_segmentator;
	ContPropMap<boost::subgraph<GraphType>, std::vector<size_t>, VERTEX> segmented_subgraph(subgraph);
	size_t new_segments_num = 0;
	kmeans_segmentator.Segment(subgraph, get(boost::vertex_index, subgraph), parts_num, segmented_subgraph, new_segments_num);

	if (new_segments_num == 0)
	{
		return;
	}
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
void CalculateSegmentGraphAndCenters(const boost::subgraph<GraphType>& graph, const SegmentsMap& segments_map,
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

	//calculate new centers
	std::vector<typename boost::graph_traits<CentersGraph>::vertex_descriptor> centers_nums(clust_num + 1);
	centers_graph = CentersGraph();
	size_t non_zero_cluster = 0;
	boost::subgraph<GraphType> graph_copy = graph;

	for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
	{
		const auto res = CalculateSegmentCenter(graph_copy, segments_map, curr_clust,
			get(boost::vertex_info_3d, graph_copy));

		if (res.second == true)
		{
			const auto center_descr = res.first;
			const auto added = add_vertex(centers_graph);
			centers_nums[curr_clust] = added;
			put(boost::vertex_info_3d, centers_graph, added, 
				get(boost::vertex_info_3d, graph, center_descr));
			put(boost::vertex_parent, centers_graph, added, center_descr);		
			++non_zero_cluster;
		}
	}

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