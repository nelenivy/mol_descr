#pragma once

#include <vector>
#include "boost/graph/graph_traits.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/graph/properties.hpp"
#include "opencv2/core/core.hpp"
#include "CommonUtilities/priority_heap.h"
#include "GraphLib/array_property_map.h"
#include "SingularPoints/mesh_types.h"

namespace molecule_descriptor
{

template <typename GraphType, typename IndexMap, typename EdgeLengthMap>
void CalcDistMatr(const GraphType& graph, const IndexMap& index_map, const EdgeLengthMap& edge_length_map, 
				  const float max_dist_mark, const float no_connection_mark, cv::Mat_<float>& dist_mat);

template <typename GraphType>
class DijkstraDistMapCalculator
{
public:
	static float kMaxDist() { return 1000000.0f; }
	static float kNoConnectionMark() { return -1.0f; }
	void Calc(const GraphType& graph, cv::Mat_<float>& dist_mat)
	{
		m_edge_weights.SetGraph(graph);
		const auto info_3d_map = get(boost::vertex_info_3d, graph);
		typename boost::graph_traits<GraphType>::edge_iterator curr_edge, end_edge;

		for (auto curr_edge = edges(graph).first, end_edge = edges(graph).second; curr_edge != end_edge; ++curr_edge)
		{
			const auto first_vert = source(*curr_edge, graph);
			const auto second_vert = target(*curr_edge, graph);
			const Point3d first_coord = info_3d_map[first_vert].Center();
			const Point3d second_coord = info_3d_map[second_vert].Center();
			m_edge_weights[*curr_edge] = cv::norm(first_coord - second_coord);
		}

		const auto index_map = get(boost::vertex_index, graph);
		CalcDistMatr(graph, index_map, m_edge_weights, kMaxDist(), kNoConnectionMark(), dist_mat);
	}
private:
	ContPropMap<GraphType, std::vector<float>, EDGE> m_edge_weights;
};

template <typename GraphType, typename IndexMap, typename EdgeLengthMap>
void CalcDistMatr(const GraphType& graph, const IndexMap& index_map, const EdgeLengthMap& edge_length_map, 
	const float max_dist_mark, const float no_connection_mark, cv::Mat_<float>& dist_mat)
{
	const int graph_size = static_cast<int>(num_vertices(graph));
	dist_mat.create(graph_size, graph_size);
		
	typedef typename boost::graph_traits<GraphType>::vertex_descriptor vertex_descriptor;
	typename boost::graph_traits<GraphType>::vertex_iterator beg_node, end_node;
	boost::tie(beg_node, end_node) = vertices(graph);
	std::vector<vertex_descriptor> graph_vertices(beg_node, end_node);

	for (auto node_iter = graph_vertices.begin(); node_iter != graph_vertices.end(); ++node_iter)
	{//for each node find shortest distances to other nodes
		auto distance_queue = MakeRisingProrityQueue(graph_vertices.begin(), graph_vertices.end(), max_dist_mark);
		const int base_ind = index_map[*node_iter];
		distance_queue.ChangeKey(base_ind, 0.0f, eThrow);

		while (distance_queue.Size() > 0)
		{
			//extract top
			const vertex_descriptor& curr_top_elem = *(distance_queue.Top().Pointer());
			float curr_min_dist =  distance_queue.Top().Key();
			distance_queue.Pop();

			for (auto neighb_iter = adjacent_vertices(curr_top_elem, graph).first, 
				end_neighb = adjacent_vertices(curr_top_elem, graph).second;
				neighb_iter != end_neighb; ++neighb_iter)
			{//update distances for neighbours
				const int neighb_ind = index_map[*neighb_iter];
				const auto curr_edge = edge(curr_top_elem, *neighb_iter, graph).first;
				distance_queue.ChangeKey(neighb_ind, curr_min_dist + edge_length_map[curr_edge], eNotThrow);
			}
			//write distance to matrix
			if (curr_min_dist >= max_dist_mark)
			{
				curr_min_dist = no_connection_mark;
			}

			const int curr_ind = index_map[curr_top_elem];
			dist_mat(base_ind, curr_ind) = dist_mat(curr_ind, base_ind) = curr_min_dist;
		}
	}
}

}