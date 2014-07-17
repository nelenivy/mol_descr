#pragma once

#include <vector>
#include <deque>

#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"

namespace molecule_descriptor
{
template <class Graph, class PropMap, class CompareFunc>
void FindLocalMaximums(const Graph& graph, const PropMap& prop_map, 
		std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& maximums, CompareFunc comp_func)
{
	maximums.clear();

	for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		bool maximum = true;

		for (auto curr_neighb = adjacent_vertices(*curr_vertice, graph).first, 
			end_neighb = adjacent_vertices(*curr_vertice, graph).second; 
			curr_neighb != end_neighb; ++curr_neighb)
		{
			if (! comp_func(abs(prop_map[*curr_vertice]), abs(prop_map[*curr_neighb])))
			{
				maximum = false;
				break;
			}
		}

		if (maximum)
		{
			maximums.push_back(*curr_vertice);
		}
	}
}

//here we search not only through local members but through neighbours on levels
template <class Graph, class PropMap>
void FindLocalMaximumsOnLevels(const Graph& graph, const std::vector<PropMap>& prop_map_levels, 
					   std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& maximums)
{
	maximums.clear();
	std::vector<size_t> max_level_count(prop_map_levels.size(), 0);

	for (size_t level = 1; level < prop_map_levels.size(); ++level)
	{
		for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			bool maximum = true;
			const size_t prev_level = level > 0? level - 1 : 0;
			const size_t next_level = std::min(level + 1, prop_map_levels.size() - 1);

			for (size_t neighb_level = prev_level; neighb_level <= next_level; ++neighb_level)
			{
				if (neighb_level != level && abs(prop_map_levels[neighb_level][*curr_vertice]) >= abs(prop_map_levels[level][*curr_vertice]))
				{
					maximum = false;
					break;
				}
				for (auto curr_neighb = adjacent_vertices(*curr_vertice, graph).first, 
					end_neighb = adjacent_vertices(*curr_vertice, graph).second; 
					curr_neighb != end_neighb; ++curr_neighb)
				{
					if (abs(prop_map_levels[neighb_level][*curr_neighb]) >= abs(prop_map_levels[level][*curr_vertice]))
					{
						maximum = false;
						break;
					}
				}
			}

			if (maximum)
			{
				max_level_count[level]++;
				maximums.push_back(*curr_vertice);
			}
		}

		std::cout << max_level_count[level] << " ";
	}
	std::cout << "\n";
}

template <class Graph, class WaveMap>
void WaveAlgorithm(const Graph& graph, const typename boost::property_traits<WaveMap>::value_type radius, 
				   const typename boost::graph_traits<Graph>::vertex_descriptor seed, WaveMap& prop_map)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef typename boost::property_traits<WaveMap>::value_type PropType;
	for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		prop_map[*curr_vertice] = 0;
	}

	std::deque<vertex_descriptor> vertices_queue;
	vertices_queue.push_back(seed);

	while (!vertices_queue.empty())
	{
		vertex_descriptor curr_vertice = vertices_queue.front();
		vertices_queue.pop_front();
		const PropType curr_prop = prop_map[curr_vertice];
		CV_Assert(curr_prop <= radius);

		if (curr_prop == radius)
		{
			continue;
		}

		for (auto neighb_it = adjacent_vertices(curr_vertice, graph).first,
			end_neighb = adjacent_vertices(curr_vertice, graph).second; neighb_it != end_neighb; ++neighb_it)
		{
			if (prop_map[*neighb_it] == 0)
			{
				prop_map[*neighb_it] = curr_prop + 1;
				vertices_queue.push_back(*neighb_it);
			}
		}
	}	
}

}
