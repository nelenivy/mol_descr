#pragma once 

#include <opencv2\core\core.hpp>
#include "boost/graph/graph_traits.hpp"
#include "GraphLib/connected_components_segmentator.h"
#include "GraphLib/array_property_map.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

template <typename Graph_1, typename CoordMap_1, typename Graph_2, typename CoordMap_2>
void CalcDistBetweenGraphs(const Graph_1& graph_1, const CoordMap_1& coord_map_1, 
						   const Graph_2& graph_2, const CoordMap_2& coord_map_2, cv::Mat_<double> dist_mat)
{
	dist_mat.create(num_vertices(graph_1), num_vertices(graph_2));
	for (auto vert_1 = vertices(graph_1).first, end_vert_1 = vertices(graph_1).second;
		vert_1 != end_vert_1; ++vert_1)
	{
		for (auto vert_2 = vertices(graph_2).first, end_vert_2 = vertices(graph_2).second;
			vert_2 != end_vert_2; ++vert_2)
		{
			cv::Point3d coord1 = coord_map_1[*vert_1];
			cv::Point3d coord2 = coord_map_2[*vert_2];
			dist_mat(*vert_1, *vert_2) = cv::norm(coord1 - coord2);
		}
	}
}
}