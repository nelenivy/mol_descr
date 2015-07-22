#pragma once

#include <vector>
#include <deque>

#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"

namespace molecule_descriptor
{
template <class Graph, class PropMap, class CompareFuncMax, class CompareFuncMin>
void FindLocalMaximumsOfAbsVal(const Graph& graph, const PropMap& prop_map, 
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
template <class Graph, class PropMap, class CompareFuncGreater, class CompareFuncLess>
void FindLocalMaximumsOnLevels(const Graph& graph, const std::vector<PropMap>& prop_map_levels, 
					   std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>>& maximums, 
					   CompareFuncGreater comp_func_greater, CompareFuncLess comp_func_less, const bool one_ring_neighb)
{
	maximums.assign(prop_map_levels.size() - 2, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>());
	std::vector<size_t> max_level_count(prop_map_levels.size(), 0);

	for (size_t level = 1; level < prop_map_levels.size() - 1; ++level)
	{
		for (auto curr_vertice = vertices(graph).first, end_vertices = vertices(graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			bool maximum = true;
			bool minimum = true;
			const size_t prev_level = level - 1;
			const size_t next_level = level + 1;

			double diff = 0;
			double num_neighb = 0;

			for (size_t neighb_level = prev_level; neighb_level <= next_level; ++neighb_level)
			{
				if (neighb_level != level)
				{
					num_neighb++;
					diff += prop_map_levels[neighb_level][*curr_vertice];
					if (! comp_func_greater(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_vertice]))
					{
						maximum = false;
					}
					if (! comp_func_less(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_vertice]))
					{
						minimum = false;
					}
				}

				for (auto curr_neighb = adjacent_vertices(*curr_vertice, graph).first, 
					end_neighb = adjacent_vertices(*curr_vertice, graph).second; 
					curr_neighb != end_neighb; ++curr_neighb)
				{
					num_neighb++;
					diff += prop_map_levels[neighb_level][*curr_neighb];

					if (! comp_func_greater(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb]))
					{
						maximum = false;
					}
					if (! comp_func_less(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb]))
					{
						minimum = false;
					}

					if (! one_ring_neighb)
					{
						for (auto curr_neighb_neighb = adjacent_vertices(*curr_neighb, graph).first, 
							end_neighb_neighb = adjacent_vertices(*curr_neighb, graph).second; 
							curr_neighb_neighb != end_neighb_neighb; ++curr_neighb_neighb)
						{
							if (*curr_neighb_neighb == *curr_vertice)
							{
								continue;
							}

							if (! comp_func_greater(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb_neighb]))
							{
								maximum = false;
							}
							if (! comp_func_less(prop_map_levels[level][*curr_vertice], prop_map_levels[neighb_level][*curr_neighb_neighb]))
							{
								minimum = false;
							}
						}
					}
				}
			}

			//const double av_diff = 
			if (maximum || minimum)
			{
				//std::cout << diff / num_neighb / prop_map_levels[level][*curr_vertice] << " ";
				max_level_count[level]++;
				maximums[level - 1].push_back(*curr_vertice);
			}
		}

		//std::cout << max_level_count[level] << "\n";
	}
	std::cout << "\n";
}

template <class PropMap, class VertDescr>
double ScalarProduct(const std::vector<double>& vect, const std::vector<PropMap>& prop, const VertDescr vert)
{
	CV_Assert(vect.size() == prop.size());
	double sum = 0, norm = 0;
	for (size_t curr_prop = 0; curr_prop < prop.size(); ++curr_prop)
	{
		sum += vect[curr_prop] * prop[curr_prop][vert];
		norm += vect[curr_prop] * vect[curr_prop];
	}
	
	return sum / sqrt(norm);
}

template <class PropMap, class VertDescr>
double Norm(const std::vector<PropMap>& prop, const VertDescr vert)
{
	double norm = 0;
	for (size_t curr_prop = 0; curr_prop < prop.size(); ++curr_prop)
	{
		norm += prop[curr_prop][vert] * prop[curr_prop][vert];
	}

	return norm;
}

template <class PropMap, class VertDescr>
double ScalarProductOneVert(const std::vector<PropMap>& prop_norm, const std::vector<PropMap>& prop, const VertDescr vert)
{
	CV_Assert(prop.size() == prop_norm.size());

	double norm = 0;
	double sum = 0;

	for (size_t curr_prop = 0; curr_prop < prop.size(); ++curr_prop)
	{
		norm += prop_norm[curr_prop][vert] * prop_norm[curr_prop][vert];
		sum += prop[curr_prop][vert] * prop_norm[curr_prop][vert];
	}

	return sum / sqrt(norm);
}
//this is not for original vectors but for the difference of original vector
template <class PropMap, class CoordMap, class PCAProjMap, class VertDescr>
double ProjectPCADiff(const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const PCAProjMap& proj_map, const VertDescr vert)
{
	cv::Mat_<double> prop_vect(1, prop_map.size() + 3);
	prop_vect(0, 0)=coord_map[vert].x;
	prop_vect(0, 1)=coord_map[vert].y;
	prop_vect(0, 2)=coord_map[vert].z;
	for (int ind = 3; ind < 3 + prop_map.size(); ++ind)
	{
		prop_vect(0, ind) = prop_map[ind - 3][vert];
	}
	prop_vect += proj_map[vert].pca.mean;
	cv::Mat_<double> proj_prop_vect;
	proj_map[vert].pca.project(prop_vect, proj_prop_vect);
	double res = proj_prop_vect.dot(proj_map[vert].vect_to_project_on)/ cv::norm(proj_map[vert].vect_to_project_on);
	return res;
}
template <class PropMap, class CoordMap, class PCAProjMap, class VertDescr>
double ProjectPCA(const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const PCAProjMap& proj_map, const VertDescr vert)
{
	cv::Mat_<double> prop_vect(1, prop_map.size() + 3);
	prop_vect(0, 0)=coord_map[vert].x;
	prop_vect(0, 1)=coord_map[vert].y;
	prop_vect(0, 2)=coord_map[vert].z;
	for (int ind = 3; ind < 3 + prop_map.size(); ++ind)
	{
		prop_vect(0, ind) = prop_map[ind - 3][vert];
	}
	cv::Mat_<double> proj_prop_vect;
	proj_map[vert].pca.project(prop_vect, proj_prop_vect);
	double res = proj_prop_vect.dot(proj_map[vert].vect_to_project_on)/ cv::norm(proj_map[vert].vect_to_project_on);
	return res;
}
//here we search not only through local members but through neighbours on levels
template <class Graph, class PropMap, class CoordMap, class PropMapProj, class CompareFuncMax, class CompareFuncMin>
void FindLocalMaximumsOnLevelsVect(const Graph& graph, const CoordMap& coord_map, const std::vector<std::vector<PropMap>>& prop_map_levels, 
								const std::vector<PropMapProj>& prop_vect_project_levels, 
							   std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>>& maximums, 
							   CompareFuncMax comp_func_max, CompareFuncMin comp_func_min, const bool use_levels, const bool one_ring_neighb)
{
	maximums.assign(prop_map_levels.size() - 2, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>());
	std::vector<size_t> max_level_count(prop_map_levels.size(), 0);
	//std::vector<double> curr_vect_to_project(prop_map_levels[0].size());

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
			/*for (size_t curr_prop = 0; curr_prop < prop_vect_project_levels[level].size(); ++curr_prop)
			{
				curr_vect_to_project[curr_prop] = prop_vect_project_levels[level][curr_prop][*curr_vertice];
			}*/
			//const double curr_norm = Norm(prop_vect_project_levels[level], *curr_vertice);

			const double curr_scalar_prod = ProjectPCADiff(prop_map_levels[level], coord_map, prop_vect_project_levels[level],*curr_vertice);
			//ScalarProduct(curr_vect_to_project, prop_map_levels[level], *curr_vertice);
			double mean_diff_neighb = 0;
			int neighb_num = 2;
			for (size_t neighb_level = prev_level; neighb_level <= next_level; ++neighb_level)
			{
				//for (size_t curr_prop = 0; curr_prop < prop_vect_project_levels[level].size(); ++curr_prop)
				{
					//curr_vect_to_project[curr_prop] = prop_vect_project_levels[neighb_level /*level*/][curr_prop][*curr_vertice];
				}
				if (neighb_level != level)
				{ 
					const double curr_neighb_prod =ProjectPCADiff(prop_map_levels[neighb_level], coord_map, 
						prop_vect_project_levels[neighb_level],*curr_vertice);
						//ScalarProductOneVert(prop_vect_project_levels[neighb_level], prop_map_levels[neighb_level], *curr_vertice);
						//ScalarProduct(curr_vect_to_project, prop_map_levels[neighb_level], *curr_vertice);

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

				for (auto curr_neighb = adjacent_vertices(*curr_vertice, graph).first, 
					end_neighb = adjacent_vertices(*curr_vertice, graph).second; 
					curr_neighb != end_neighb; ++curr_neighb)
				{
					const double curr_neighb_prod = ProjectPCADiff(prop_map_levels[neighb_level], coord_map, 
						prop_vect_project_levels[neighb_level],*curr_neighb);
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

					if (!one_ring_neighb)
					{
						for (auto curr_neighb_neighb = adjacent_vertices(*curr_neighb, graph).first, 
							end_neighb_neighb = adjacent_vertices(*curr_neighb, graph).second; 
							curr_neighb_neighb != end_neighb_neighb; ++curr_neighb_neighb)
						{
							if (*curr_neighb_neighb == *curr_vertice)
							{
								continue;
							}
							const double curr_neighb_prod =  ProjectPCADiff(prop_map_levels[neighb_level], coord_map, 
								prop_vect_project_levels[neighb_level],*curr_neighb_neighb);
								//ScalarProductOneVert(prop_vect_project_levels[neighb_level], prop_map_levels[neighb_level], *curr_neighb_neighb);
								//ScalarProduct(curr_vect_to_project, prop_map_levels[neighb_level], *curr_neighb_neighb);
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
				}
			}
			
			const double mean_diff_to_neighb = mean_diff_neighb / neighb_num / curr_scalar_prod;
			if ((maximum || minimum)/* && std::abs(mean_diff_to_neighb) > 0.4*/)
			{
				max_level_count[level]++;
				maximums[level - 1].push_back(*curr_vertice);
				//std::cout << "lev " << level - 1 << " val " << mean_diff_neighb / neighb_num / curr_scalar_prod << " ";
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
