#include "mesh_constructor.h"
#include <algorithm>
#include <iostream>

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

void MeshKeeper::ConstructMesh(const vector<cv::Point3d>& vertices, 
							   const vector<cv::Point3d>& normals, 
							   const vector<Point3i>& triangles)
{
	CV_Assert(vertices.size() == normals.size());
	Clear();
	FillVertices(vertices, normals);
	FillTriangles(triangles);
	FillMeshVertices();	
	FillMeshTriangles();
}

void MeshKeeper::FillVertices(const vector<cv::Point3d>& vertices, const vector<cv::Point3d>& normals)
{
	CV_Assert(vertices.size() == normals.size());
	//fill vertices
	m_vertices.resize(vertices.size(), Vertice(cv::Point3d(), cv::Point3d()));
	auto pmap = get(boost::vertex_info_3d_t(), m_vertices_graph);
	for (int ind = 0; ind < m_vertices.size(); ind++)
	{
		auto added_vertex = boost::add_vertex(m_vertices_graph);
		pmap[added_vertex].coord = vertices[ind];
		pmap[added_vertex].normal = normals[ind];		
	}
}

void MeshKeeper::FillTriangles(const vector<Point3i>& triangles)
{
	CV_Assert(m_vertices.size());
	m_triangles_with_indexes = triangles;
	auto pmap = get(boost::vertex_info_3d_t(), m_triangles_graph);

	for (int ind = 0; ind < m_triangles_with_indexes.size(); ind++)
	{
		auto added_triangle = add_vertex(m_triangles_graph);
		const Point3i& curr_triangle = triangles[ind];
		const VertexDescriptor a = vertex(curr_triangle.x, m_vertices_graph);
		const VertexDescriptor b = vertex(curr_triangle.y, m_vertices_graph);
		const VertexDescriptor c = vertex(curr_triangle.z, m_vertices_graph);
		pmap[added_triangle] = MeshTriangle(a, b, c, m_vertices_graph);
	}
}

void MeshKeeper::FindTriangleNeighboursOfEachVertex()
{
	CV_Assert(m_vertices.size() && m_triangles_with_indexes.size());

	m_neighbours_numbers.resize(m_vertices.size());
	fill(m_neighbours_numbers.begin(), m_neighbours_numbers.end(), 0);
	//find neighbours triangles number for each vertex	
	for (auto iter = m_triangles_with_indexes.begin(); iter != m_triangles_with_indexes.end(); iter++)
	{
		const Point3i& curr_triangle = *iter;
		m_neighbours_numbers[curr_triangle.x]++;
		m_neighbours_numbers[curr_triangle.y]++;
		m_neighbours_numbers[curr_triangle.z]++;
	}	

	m_triangle_neighb_of_vertex.resize(m_vertices.size());

	for (int ind = 0; ind < m_triangle_neighb_of_vertex.size(); ind++)
	{
		m_triangle_neighb_of_vertex[ind].reserve(m_neighbours_numbers[ind]);
	}
	//find neighbours triangles for each vertex
	for (int ind = 0; ind < m_triangles_with_indexes.size(); ind++)
	{
		const Point3i& curr_triangle_indexes = m_triangles_with_indexes[ind];

		m_triangle_neighb_of_vertex[curr_triangle_indexes.x].push_back(ind);
		m_triangle_neighb_of_vertex[curr_triangle_indexes.y].push_back(ind);
		m_triangle_neighb_of_vertex[curr_triangle_indexes.z].push_back(ind);
	}
}

void MeshKeeper::FillMeshVertices()
{
	CV_Assert(m_vertices.size() && m_triangles_with_indexes.size());
	
	for (auto iter = m_triangles_with_indexes.begin(); iter != m_triangles_with_indexes.end(); iter++)
	{
		const size_t vert_1_ind = iter->x;
		const int vert_2_ind = iter->y;
		const int vert_3_ind = iter->z;

		VertexDescriptor vert_1 = boost::vertex(vert_1_ind, m_vertices_graph);
		VertexDescriptor vert_2 = boost::vertex(vert_2_ind, m_vertices_graph);
		VertexDescriptor vert_3 = boost::vertex(vert_3_ind, m_vertices_graph);

		if (!edge(vert_1, vert_2, m_vertices_graph).second)
		{
			add_edge(vert_1, vert_2, m_vertices_graph);
		}
		if (!edge(vert_1, vert_3, m_vertices_graph).second)
		{
			add_edge(vert_1, vert_3, m_vertices_graph);
		}
		if (!edge(vert_2, vert_3, m_vertices_graph).second)
		{
			add_edge(vert_2, vert_3, m_vertices_graph);
		}
	}	
}

void MeshKeeper::FillMeshTriangles()
{
	CV_Assert(m_vertices.size() && m_triangles_with_indexes.size());

	FindTriangleNeighboursOfEachVertex();
	const auto pmap = get(boost::vertex_info_3d, m_triangles_graph);
	for (int ind_vect = 0; ind_vect < m_triangle_neighb_of_vertex.size(); ind_vect++)
	{
		const vector<int>& curr_neighbours = m_triangle_neighb_of_vertex[ind_vect];

		for (int ind_curr = 0; ind_curr < curr_neighbours.size(); ind_curr++)
		{//find neighbour triangles for each triangle
			const int curr_triangle_ind = curr_neighbours[ind_curr];
			const TriangleDescriptor curr_triangle = boost::vertex(curr_triangle_ind, m_triangles_graph);

			for (int ind_neighb = 0; ind_neighb < curr_neighbours.size(); ind_neighb++)
			{
				if (ind_neighb == ind_curr) 
				{ 
					continue; 
				}

				const int neighb_triangle_ind = curr_neighbours[ind_neighb];
				const TriangleDescriptor neighb_triangle = boost::vertex(neighb_triangle_ind, m_triangles_graph);

				if (CountSameVertices(pmap[neighb_triangle], pmap[curr_triangle]) == 2)
				{
					if (!edge(neighb_triangle, curr_triangle, m_triangles_graph).second)
					{
						add_edge(neighb_triangle, curr_triangle, m_triangles_graph);
					}
				}
			}
		}
	}
}

void MeshKeeper::Clear()
{
	m_vertices_graph = VerticesGraph();
	m_triangles_graph = TrianglesGraph();
	m_vertices.clear();
	m_triangles_with_indexes.clear();
	m_triangle_neighb_of_vertex.clear();
	m_neighbours_numbers.clear();
}

}