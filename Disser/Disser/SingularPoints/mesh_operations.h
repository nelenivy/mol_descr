#pragma once

#include "SingularPoints/mesh_types.h"
#include "GraphLib/graph_filter.h"
#include "GraphLib/proxy_property_map.h"

namespace molecule_descriptor
{

template <class FilterKernel>
void FilterVerticesGraph(const VerticesGraph& graph_in, const FilterKernel& kernel, VerticesGraph& graph_out)
{
	graph_out = graph_in;

	auto map_3d_in = get(boost::vertex_info_3d, graph_in);
	auto map_3d_out = get(boost::vertex_info_3d, graph_out);

	auto map_3d_in_coord = GetProxyPropMap(map_3d_in, GetCoord<Vertice>());
	auto map_3d_out_coord = GetProxyPropMap(map_3d_out, GetCoord<Vertice>());
	FilterGraphDist(kernel, graph_in, map_3d_in_coord, map_3d_in_coord, map_3d_out_coord);

	auto map_3d_in_norm = GetProxyPropMap(map_3d_in, GetNormal<Vertice>());
	auto map_3d_out_norm = GetProxyPropMap(map_3d_out, GetNormal<Vertice>());
	FilterGraphDist(kernel, graph_in, map_3d_in_coord, map_3d_in_norm, map_3d_out_norm);

	//normalize normal vectors
	for (auto curr_vertice = vertices(graph_out).first, end_vertices = vertices(graph_out).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const double curr_norm = norm(map_3d_out[*curr_vertice].normal);

		if (curr_norm > 0.0)
		{
			map_3d_out[*curr_vertice].normal *= 1.0 / curr_norm;
		}
		else
		{
			map_3d_out[*curr_vertice].normal = cv::Point3d(1, 0, 0);
		}
	}
}

template <class FilterKernel>
void FilterMesh(const Mesh& mesh_in, const FilterKernel& kernel, Mesh& mesh_out)
{
	FilterVerticesGraph(mesh_in.vertices, kernel, mesh_out.vertices);
	mesh_out.triangles = mesh_in.triangles;
	auto triangle_map_out = get(boost::vertex_info_3d, mesh_out.triangles);

	for (auto tr_iter = vertices(mesh_out.triangles).first, tr_end = vertices(mesh_out.triangles).second;
		tr_iter != tr_end; ++tr_iter)
	{
		auto& curr_triangle = triangle_map_out[*tr_iter];
		curr_triangle.m_graph = &mesh_out.vertices;
	}
}
}