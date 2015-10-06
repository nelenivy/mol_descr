#pragma once

#include <vector>
#include <utility>
#include <deque>

#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"

namespace molecule_descriptor
{
//find vertices from graph_2 which are within distance threshold from vertice vert_1 of graph_1
template <class Graph1Descr, class Graph2, class Graph12DistMap, class MaskMap>
size_t GetVerticesWithinDistMask(const Graph1Descr vert_1, const Graph2& graph_2, const Graph12DistMap& dist_map, 
						   const double dist_thresh, MaskMap& mask_map)
{
	if (mask_map.Size() != boost::num_vertices(graph_2))
	{
		mask_map.SetGraph(graph_2);
	}
	size_t vert_num_in_mask = 0;
	for (auto neighb_vertice = vertices(graph_2).first, end_neighb_vertices = vertices(graph_2).second; 
		neighb_vertice != end_neighb_vertices; ++neighb_vertice)
	{

		if (dist_map(vert_1, *neighb_vertice) <= dist_thresh)
		{
			mask_map[*neighb_vertice] = 255;
			++vert_num_in_mask;
		}
		else
		{
			mask_map[*neighb_vertice] = 0;
		}
	}

	return vert_num_in_mask;
}
template <class GraphDescr, class Graph, class GraphDistMap, class MaskMap>
size_t GetVerticesWithinDistMaskPlusAdjacent(const GraphDescr vert, const Graph& graph, const GraphDistMap& dist_map, 
								 const double dist_thresh, MaskMap& mask_map)
{
	size_t vert_num_in_mask = GetVerticesWithinDistMask( vert, graph, dist_map, dist_thresh, mask_map);

	for (auto neighb_it = adjacent_vertices(vert, graph).first,
		end_neighb = adjacent_vertices(vert, graph).second; 
		neighb_it != end_neighb; ++neighb_it)
	{
		if (mask_map[*neighb_it] == 0)
		{
			++vert_num_in_mask;
			mask_map[*neighb_it] == 255;
		}
	}

	return vert_num_in_mask;
}

//find vertices from graph_2 which are within distance threshold from vertice vert_1 of graph_1
template <class Graph1Descr, class Graph2, class Graph12DistMap, class Graph2Descr>
void GetVerticesWithinDist(const Graph1Descr vert_1, const Graph2& graph_2, const Graph12DistMap& dist_map, 
								 const double dist_thresh, std::vector<Graph2Descr>& vertices_in_mask)
{
	vertices_in_mask.clear();

	for (auto neighb_vertice = vertices(graph_2).first, end_neighb_vertices = vertices(graph_2).second; 
		neighb_vertice != end_neighb_vertices; ++neighb_vertice)
	{

		if (dist_map(vert_1, *neighb_vertice) <= dist_thresh)
		{
			vertices_in_mask.push_back(*neighb_vertice);
		}
	}
}

//find vertices from graph_2 which are within distance threshold from vertice vert_1 of graph_1
template <class GraphDescr, class Graph, class GraphDistMap>
void GetVerticesWithinDistPlusAdjacent(const GraphDescr vert, const Graph& graph, const GraphDistMap& dist_map, 
						   const double dist_thresh, std::vector<GraphDescr>& vertices_in_mask)
{
	GetVerticesWithinDist(vert, graph, dist_map, dist_thresh, vertices_in_mask);

	for (auto neighb_it = adjacent_vertices(vert, graph).first,
		end_neighb = adjacent_vertices(vert, graph).second; 
		neighb_it != end_neighb; ++neighb_it)
	{
		if (vertices_in_mask.end() == std::find(vertices_in_mask.begin(), vertices_in_mask.end(), *neighb_it))
		{
			vertices_in_mask.push_back(*neighb_it);
		}
	}
}

template <class Graph, class PropMap, class CompareFuncMax, class CompareFuncMin>
void FindLocalMaxAndMin(const Graph& graph, const PropMap& prop_map, 
		std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& maximums, 
		CompareFuncMax comp_func_max, CompareFuncMin comp_func_min)
{
	maximums.clear();
	int low_maximums = 0;
	for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		bool maximum = true;
		bool minimum = true;
		double diff = 0;
		for (auto curr_neighb = adjacent_vertices(*curr_vertice, graph).first, 
			end_neighb = adjacent_vertices(*curr_vertice, graph).second; 
			curr_neighb != end_neighb; ++curr_neighb)
		{
			diff += abs(prop_map[*curr_vertice]) - abs(prop_map[*curr_neighb]);
			if (! comp_func_max(prop_map[*curr_vertice], prop_map[*curr_neighb]))
			{
				maximum = false;
			}
			if (! comp_func_min(prop_map[*curr_vertice], prop_map[*curr_neighb]))
			{
				minimum = false;
			}

			/*for (auto curr_neighb_neighb = adjacent_vertices(*curr_neighb, graph).first, 
				end_neighb_neighb = adjacent_vertices(*curr_neighb, graph).second; 
				curr_neighb_neighb != end_neighb_neighb; ++curr_neighb_neighb)
			{
				if (*curr_neighb_neighb == *curr_vertice)
				{
					continue;
				}
				diff += abs(prop_map[*curr_vertice]) - abs(prop_map[*curr_neighb]);
				if (! comp_func_max(prop_map[*curr_vertice], prop_map[*curr_neighb_neighb]))
				{
					maximum = false;
				}
				if (! comp_func_min(prop_map[*curr_vertice], prop_map[*curr_neighb_neighb]))
				{
					minimum = false;
				}
			}*/
		}

		diff /= (adjacent_vertices(*curr_vertice, graph).second - adjacent_vertices(*curr_vertice, graph).first);
		if (maximum || minimum/*&& diff > 0.2*/)
		{
//			std::cout << diff /  << " ";
			maximums.push_back(*curr_vertice);
		}
		else if (maximum)
		{
			low_maximums++;
		}
	}
}

//here we search not only through local members but through neighbours on levels
template <class Graph, class PropMap, class CompareFuncGreater, class CompareFuncLess, class VertDistMap>
void FindLocalMaximumsOnLevels(const Graph& graph, const std::vector<PropMap>& prop_map_levels, 
					   CompareFuncGreater comp_func_greater, const CompareFuncLess comp_func_less, 
					   const VertDistMap& vert_dist_map, const double dist_thresh,
					   const bool use_levels, 
					   std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>>& maximums)
{
	maximums.assign(prop_map_levels.size() - 2, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>());
	std::vector<size_t> max_level_count(prop_map_levels.size(), 0);
	std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vertices_within_dist;

	for (size_t level = 1; level < prop_map_levels.size() - 1; ++level)
	{
		for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			bool maximum = true;
			bool minimum = true;
			const size_t prev_level = use_levels ? level - 1 : level;
			const size_t next_level = use_levels ? level + 1 : level;
			GetVerticesWithinDistPlusAdjacent(*curr_vertice, graph, vert_dist_map, dist_thresh, vertices_within_dist);

			for (size_t neighb_level = prev_level; neighb_level <= next_level; ++neighb_level)
			{
				for (auto curr_neighb = vertices_within_dist.begin(), end_neighb = vertices_within_dist.end(); 
					curr_neighb != end_neighb; ++curr_neighb)
				{
					if ((neighb_level == level && *curr_neighb == *curr_vertice) || 
						(neighb_level != level && *curr_neighb != *curr_vertice))
					{
						continue;
					}
					if (! comp_func_greater(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb]))
					{
						maximum = false;
					}
					if (! comp_func_less(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb]))
					{
						minimum = false;
					}
				}
			}

			if (maximum || minimum)
			{
				max_level_count[level]++;
				maximums[level - 1].push_back(*curr_vertice);
			}
		}

		std::cout << max_level_count[level] << " ";
	}
	std::cout << "\n";
}

//here we search not only through local members but through neighbours on levels
template <class Graph, class PropMap, class CoordMap, class PropMapProj, 
class CompareFuncMax, class CompareFuncMin, class VertDistMap>
void FindLocalMaximumsOnLevelsVect(const Graph& graph, const CoordMap& coord_map, const std::vector<std::vector<PropMap>>& prop_map_levels, 
								const std::vector<PropMapProj>& prop_vect_project_levels, 
								const VertDistMap& vert_dist_map, const double dist_thresh,
								CompareFuncMax comp_func_max, CompareFuncMin comp_func_min, 
								const bool use_levels, const bool use_projecter_from_center,
							   std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>>& maximums)
{
	maximums.assign(prop_map_levels.size() - 2, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>());
	std::vector<size_t> max_level_count(prop_map_levels.size(), 0);
	std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vertices_within_dist;

	for (size_t level = 1; level < prop_map_levels.size() - 1; ++level)
	{
		for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			bool maximum = true;
			bool minimum = true;
			const size_t prev_level = use_levels ? level - 1 : level;
			const size_t next_level = use_levels ? level + 1 : level;
			//vector on which we will project
			const double curr_scalar_prod = ProjectPCADiff(prop_map_levels[level], coord_map, prop_vect_project_levels[level],*curr_vertice);

			GetVerticesWithinDistPlusAdjacent(*curr_vertice, graph, vert_dist_map, dist_thresh, vertices_within_dist);

			double mean_diff_neighb = 0;
			int neighb_num = 2;
			for (size_t neighb_level = prev_level; neighb_level <= next_level; ++neighb_level)
			{
				if (neighb_level != level)
				{ 
					const double curr_neighb_prod =ProjectPCADiff(prop_map_levels[neighb_level], coord_map, 
						prop_vect_project_levels[neighb_level],*curr_vertice);

					mean_diff_neighb += curr_scalar_prod - curr_neighb_prod;
					if(! comp_func_max(curr_scalar_prod, curr_neighb_prod))
					{
						maximum = false;
					}
					if(! comp_func_min(curr_scalar_prod, curr_neighb_prod))
					{
						minimum = false;
					}
				}

				for (auto curr_neighb = vertices_within_dist.begin(), end_neighb = vertices_within_dist.end(); 
					curr_neighb != end_neighb; ++curr_neighb)
				{
					if (*curr_neighb == *curr_vertice)
					{
						continue;
					}

					const auto& used_projecter = use_projecter_from_center ? 
						prop_vect_project_levels[neighb_level][*curr_vertice] : prop_vect_project_levels[neighb_level][*curr_neighb];
					const double curr_neighb_prod = ProjectPCADiffDifferentVert(prop_map_levels[neighb_level], coord_map, 
						*curr_neighb, used_projecter);
						//ScalarProductOneVert(prop_vect_project_levels[neighb_level], prop_map_levels[neighb_level], *curr_neighb);
						//ScalarProduct(curr_vect_to_project, prop_map_levels[neighb_level], *curr_neighb);
					mean_diff_neighb += curr_scalar_prod - curr_neighb_prod;
					neighb_num++;
					if (! comp_func_max(curr_scalar_prod, curr_neighb_prod))
					{
						maximum = false;
					}
					if (! comp_func_min(curr_scalar_prod, curr_neighb_prod))
					{
						minimum = false;
					}					
				}
			}
			
			const double mean_diff_to_neighb = mean_diff_neighb / neighb_num / curr_scalar_prod;
			if ((maximum || minimum))
			{
				max_level_count[level]++;
				maximums[level - 1].push_back(*curr_vertice);
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
