#pragma once

#include "k_means_segmentator.h"
#include "connected_components_segmentator.h"
#include "segmentation_types.h"

namespace molecule_descriptor
{
template <typename NodeType>
void CreateSubgraph(const std::vector<GraphNode<NodeType>>& graph, std::vector<GraphNode<NodeType>>& subgraph);
template <typename NodeType>
void SplitSegment(const std::vector<GraphNode<NodeType>>& graph, const size_t segment_to_split, const size_t parts_num, size_t& segments_num);
template <typename NodeType>
void CalculateSegmentCenters(const std::vector<GraphNode<NodeType>>& graph, std::vector<NodeType>& centers);

template <typename NodeType>
void CreateSubgraph(const std::vector<GraphNode<NodeType>>& graph, std::vector<GraphNode<NodeType>>& subgraph)
{
	typedef PtrToSubgraph<NodeType> PtrToSubgraphCurr;
	typedef PtrToGraph<NodeType> PtrToGraphCurr;

	//count vertices which hare marked with mask
	struct IfNodeIsMarked
	{
		bool operator()(const GraphNode<NodeType>& node)
		{
			return (node.attr.Get<SubgraphMask>() == true);
		}
	} if_node_is_marked;

	int nodes_to_copy_num = static_cast<int>(count_if(graph.begin(), graph.end(), if_node_is_marked));
	subgraph.resize(nodes_to_copy_num);

	//copy marked nodes	
	{
		auto subgraph_iter = subgraph.begin();
		for (auto graph_iter = graph.begin(); graph_iter != graph.end(); ++graph_iter)
		{
			if (if_node_is_marked(*graph_iter))
			{
				subgraph_iter->element = graph_iter->element;
				subgraph_iter->attr.Add<PtrToGraphCurr>();
				subgraph_iter->attr.Get<PtrToGraphCurr>() = &*graph_iter;
				graph_iter->attr.Add<PtrToSubgraphCurr>();
				graph_iter->attr.Get<PtrToSubgraphCurr>() = &*subgraph_iter;
				++subgraph_iter;
			}
			else
			{
				graph_iter->attr.Add<PtrToSubgraphCurr>();
				graph_iter->attr.Get<PtrToSubgraphCurr>() = nullptr;
			}
		}
	}
	//fill neighbours	
	for (auto subgraph_iter = subgraph.begin(); subgraph_iter != subgraph.end(); ++subgraph_iter)
	{
		const GraphNode<NodeType>* graph_item = subgraph_iter->attr.Get<PtrToGraphCurr>();

		for (auto neighb_iter = graph_item->neighbours.begin(); neighb_iter != graph_item->neighbours.end(); ++neighb_iter)
		{				
			GraphNode<NodeType>* subgraph_neighb = (*neighb_iter)->attr.Get<PtrToSubgraphCurr>();

			if (subgraph_neighb == nullptr)
			{
				continue;
			}

			subgraph_iter->neighbours.push_back(subgraph_neighb);
		}
	}

	//delete property
	for (auto graph_iter = graph.begin(); graph_iter != graph.end(); ++graph_iter)
	{
		graph_iter->attr.Delete<PtrToSubgraphCurr>();
	}
}

template <typename NodeType>
void SplitSegment(const std::vector<GraphNode<NodeType>>& graph, const size_t segment_to_split, const size_t parts_num, size_t& segments_num)
{
	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{//mark necessary nodes
		node_iter->attr.Add<SubgraphMask>();
		node_iter->attr.Get<SubgraphMask>() = 
			node_iter->attr.Get<SegmentNumProp>() == segment_to_split ? true : false;
	}

	struct SubgraphMaskDeleter//struct to automaticaly delete subgraph mask
	{
		SubgraphMaskDeleter(const std::vector<GraphNode<NodeType>>& graph) : m_graph(graph) { }
		~SubgraphMaskDeleter()
		{
			for (auto node_iter = m_graph.begin(); node_iter != m_graph.end(); ++node_iter)
			{//
				node_iter->attr.Delete<SubgraphMask>();
			}
		}

		const std::vector<GraphNode<NodeType>>& m_graph;
	};

	SubgraphMaskDeleter deleter(graph);
	//segment using k_means
	std::vector <GraphNode<NodeType>> subgraph;
	CreateSubgraph(graph, subgraph);

	if (subgraph.size() == 0)
	{
		return;
	}

	KMeansSegmentator<NodeType> kmeans_segmentator;
	std::vector<size_t> segmented_subgraph;
	size_t new_segments_num = 0;
	kmeans_segmentator.Segment(subgraph, parts_num, segmented_subgraph, new_segments_num);
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

	for (size_t ind = 0; ind < subgraph.size(); ++ind)
	{
		const size_t curr_segm_num = segments_map[segmented_subgraph[ind]];
		const GraphNode<NodeType>& parent = *(subgraph[ind].attr.Get<PtrToGraph<NodeType>>());
		parent.attr.Get<SegmentNumProp>() = curr_segm_num;
	}
}

template <typename NodeType>
void CalculateSegmentCenters(const std::vector<GraphNode<NodeType>>& graph, std::vector<NodeType>& centers)
{
	//find clusters number
	size_t clust_num = 0;

	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{
		const size_t curr_segm_num = node_iter->attr.Get<SegmentNumProp>();

		if (curr_segm_num > clust_num)
		{
			clust_num = curr_segm_num;
		}
	}

	//calculate new centers as mean of the cluster coordinates
	std::vector<cv::Point3d> actual_centers(clust_num + 1, cv::Point3d(0.0, 0.0, 0.0));//clust + 1 because clusters numbers go from 1
	std::vector<size_t> curr_clust_sizes(clust_num + 1, 0);

	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{
		const size_t curr_clust_num = node_iter->attr.Get<SegmentNumProp>();

		if (curr_clust_num > 0)
		{
			actual_centers[curr_clust_num] += node_iter->element->Center();
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
	std::vector<std::pair<double, NodeType*>> min_dist_to_center(clust_num + 1, make_pair(kMaxDist, nullptr));
	//minimum distance on the first place, index on the second

	for (auto node_iter = graph.begin(); node_iter != graph.end(); ++node_iter)
	{
		const size_t curr_clust_num = node_iter->attr.Get<SegmentNumProp>();

		if (curr_clust_num > 0)
		{
			const double curr_dist_to_center = cv::norm(node_iter->element->Center() - actual_centers[curr_clust_num]);

			if (curr_dist_to_center < min_dist_to_center[curr_clust_num].first)
			{
				min_dist_to_center[curr_clust_num] = make_pair(curr_dist_to_center, node_iter->element);
			}
		}
	}

	//write found nodes to the output vector
	centers.resize(clust_num);
	size_t non_zero_cluster = 0;
	for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
	{
		if (curr_clust_sizes[curr_clust])
		{
			centers[non_zero_cluster++] = *(min_dist_to_center[curr_clust].second);
		}
	}

	centers.resize(non_zero_cluster);
}

}