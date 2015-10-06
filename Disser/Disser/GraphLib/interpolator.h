#pragma once

#include <stdint.h>
#include <deque>
#include "boost/graph/properties.hpp"
#include "boost/graph/adjacency_iterator.hpp"
#include "boost/graph/adjacency_list.hpp"
#include <boost/graph/graph_traits.hpp>

#include "GraphLib/array_property_map.h"

namespace molecule_descriptor
{

template <typename GraphT, class PropT>
class InterpolatorFromBorders
{
public:
	/** \brief interpolates values in the pixels, which color values are less or equal than bad_depth.
	Pixels only in area_to_process is processed*/ 
	template <class BadVerticesMask, class PropMap>
	void Interpolate(const GraphT& graph, const BadVerticesMask& bad_vert_mask, 
		const uint8_t bad_vertice_mark, PropMap& interp_prop_map);
private:
	void Init(const GraphT& graph, const PropT& zero_prop);
	/** Buffers for the algorithm*/
	typedef ContPropMap<GraphT, std::vector<uint8_t>, VERTEX> UInt8PropMap;
	UInt8PropMap m_vertices_to_process_mask;
	UInt8PropMap m_counter_of_visits; // because its value can't be more than 8, because pixel have maximum 8 neighbour pixels
	typedef ContPropMap<GraphT, std::vector<PropT>, VERTEX> PropTPropMap;
	PropTPropMap m_sum_of_props;
	typedef typename boost::graph_traits<GraphT>::vertex_descriptor VertexDescr;
	std::deque<VertexDescr> m_vertices_queue;
};

template <class GraphT, class PropT>
template <class BadVerticesMask, class PropMap>
void InterpolatorFromBorders<GraphT, PropT>::Interpolate(const GraphT& graph, const BadVerticesMask& bad_vert_mask, 
									   const uint8_t bad_vertice_mark, PropMap& interp_prop_map)
{
	//check if size changed
	const PropT zero_prop = interp_prop_map[*(vertices(graph).first)]- interp_prop_map[*(vertices(graph).first)];
	Init(graph, zero_prop);

	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert;
		++curr_vert)
	{
		m_vertices_to_process_mask[*curr_vert] = bad_vert_mask[*curr_vert];

		if (bad_vert_mask[*curr_vert] == bad_vertice_mark)				
		{
			continue;
		}
		//find if it's border vertices
		bool is_border_vert = false;

		for (auto curr_neighb = adjacent_vertices(*curr_vert, graph).first, end_neighb = adjacent_vertices(*curr_vert, graph).second;
			curr_neighb != end_neighb; ++curr_neighb)
		{
			if (bad_vert_mask[*curr_neighb] == bad_vertice_mark)
			{
				is_border_vert = true;
				break;
			}
		}

		if (is_border_vert)
		{
			m_vertices_queue.push_back(*curr_vert);	
		}
	}

	//border of parts of processed pixels
	const VertexDescr kBorderInQueue = VertexDescr(-1);
	m_vertices_queue.push_back(kBorderInQueue);

	//main loop, where all bad pixels are processed
	while (1)
	{
		const VertexDescriptor vert_to_process = m_vertices_queue.front();
		m_vertices_queue.pop_front();

		if (vert_to_process == kBorderInQueue)
		{			
			if (m_vertices_queue.empty())
			{
				break;
			}

			//mark all pixels in queue as processed and assign a value to the processed pixels
			for (size_t queue_index = 0; queue_index < m_vertices_queue.size(); ++queue_index)
			{
				const VertexDescr curr_vert = m_vertices_queue[queue_index];
				m_vertices_to_process_mask[curr_vert] = bad_vertice_mark + 1;
				interp_prop_map[curr_vert] = m_sum_of_props[curr_vert] * (1.0 / m_counter_of_visits[curr_vert]);
			}

			m_vertices_queue.push_back(kBorderInQueue);
			continue;
		}

		for (auto curr_neighb = adjacent_vertices(vert_to_process, graph).first, 
			end_neighb = adjacent_vertices(vert_to_process, graph).second;
			curr_neighb != end_neighb; ++curr_neighb)
		{
			if (m_vertices_to_process_mask[*curr_neighb] != bad_vertice_mark)
			{
				continue;
			}				
			//if we haven't processed this pixel yet
			if (m_counter_of_visits[*curr_neighb] == 0)
			{
				m_vertices_queue.push_back(*curr_neighb);						
			}

			m_sum_of_props[*curr_neighb] += interp_prop_map[vert_to_process];
			++m_counter_of_visits[*curr_neighb];
		}
	}
}

template <class GraphT, class PropT>
void InterpolatorFromBorders<GraphT, PropT>::Init(const GraphT& graph, const PropT& zero_prop)
{
	m_vertices_to_process_mask.SetGraph(graph);
	m_counter_of_visits.SetGraph(graph);
	m_sum_of_props.SetGraph(graph);

	std::fill(m_sum_of_props.begin(), m_sum_of_props.end(), zero_prop);
}

}
