#include "mesh_constructor.h"
#include <algorithm>
#include <iostream>

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

void MeshKeeper::ConstructMesh(const vector<cv::Point3d>& vertices, const vector<cv::Point3d>& normals, const vector<Point3i>& triangles)
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
	m_mesh_vertices.resize(vertices.size());

	for (int ind = 0; ind < m_vertices.size(); ind++)
	{
		m_vertices[ind].coord = vertices[ind];
		m_vertices[ind].normal = normals[ind];
		m_mesh_vertices[ind].element = &m_vertices[ind];
		m_mesh_vertices[ind].neighbours.reserve(20);
	}
}

void MeshKeeper::FillTriangles(const vector<Point3i>& triangles)
{
	CV_Assert(m_vertices.size());
	m_triangles.resize(triangles.size());
	m_mesh_triangles.resize(triangles.size());
	m_triangles_with_indexes = triangles;

	for (int ind = 0; ind < m_triangles.size(); ind++)
	{
		const Point3i& curr_triangle = triangles[ind];
		m_triangles[ind] = Triangle(&m_vertices[curr_triangle.x],
									&m_vertices[curr_triangle.y],
									&m_vertices[curr_triangle.z] );
		m_mesh_triangles[ind].element = &m_triangles[ind];
		m_mesh_triangles[ind].neighbours.reserve(3);
	}
}

void MeshKeeper::FindTriangleNeighboursOfEachVertex()
{
	CV_Assert(m_vertices.size() && m_triangles.size() && m_triangles_with_indexes.size());

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
	CV_Assert(m_vertices.size() && m_triangles.size() && m_triangles_with_indexes.size() && m_mesh_vertices.size());
	
	for (auto iter = m_triangles_with_indexes.begin(); iter != m_triangles_with_indexes.end(); iter++)
	{
		const int vert_1_ind = iter->x;
		const int vert_2_ind = iter->y;
		const int vert_3_ind = iter->z;

		MeshVertice* vert_1_ptr = &m_mesh_vertices[vert_1_ind];
		MeshVertice* vert_2_ptr = &m_mesh_vertices[vert_2_ind];
		MeshVertice* vert_3_ptr = &m_mesh_vertices[vert_3_ind];

		m_mesh_vertices[vert_1_ind].neighbours.push_back(vert_2_ptr);
		m_mesh_vertices[vert_1_ind].neighbours.push_back(vert_3_ptr);

		m_mesh_vertices[vert_2_ind].neighbours.push_back(vert_1_ptr);
		m_mesh_vertices[vert_2_ind].neighbours.push_back(vert_3_ptr);

		m_mesh_vertices[vert_3_ind].neighbours.push_back(vert_1_ptr);
		m_mesh_vertices[vert_3_ind].neighbours.push_back(vert_2_ptr);
	}	

	//remove the same pointers
	for (auto iter = m_mesh_vertices.begin(); iter != m_mesh_vertices.end(); iter++)
	{
		std::vector<MeshVertice*>& curr_neighbours = iter->neighbours;
		std::sort(curr_neighbours.begin(), curr_neighbours.end());
		curr_neighbours.erase(
			std::unique(curr_neighbours.begin(), curr_neighbours.end()),
			curr_neighbours.end());
	}	
}

void MeshKeeper::FillMeshTriangles()
{
	CV_Assert(m_vertices.size() && m_triangles.size() && m_triangles_with_indexes.size() 
		&& m_mesh_vertices.size() && m_mesh_triangles.size());

	FindTriangleNeighboursOfEachVertex();

	for (int ind_vect = 0; ind_vect < m_triangle_neighb_of_vertex.size(); ind_vect++)
	{
		const vector<int>& curr_neighbours = m_triangle_neighb_of_vertex[ind_vect];

		for (int ind_curr = 0; ind_curr < curr_neighbours.size(); ind_curr++)
		{//find neighbour triangles for each triangle
			const int curr_triangle_ind = curr_neighbours[ind_curr];
			const Triangle& curr_triangle = m_triangles[curr_triangle_ind];

			for (int ind_neighb = 0; ind_neighb < curr_neighbours.size(); ind_neighb++)
			{
				if (ind_neighb == ind_curr) 
				{ 
					continue; 
				}

				const int neighb_triangle_ind = curr_neighbours[ind_neighb];
				const Triangle& neighb_triangle = m_triangles[neighb_triangle_ind];

				if (CountSameVertices(curr_triangle, neighb_triangle) == 2)
				{
					MeshTriangle& curr_mesh_triangle = m_mesh_triangles[curr_triangle_ind];
					MeshTriangle& neighb_mesh_triangle = m_mesh_triangles[neighb_triangle_ind];
					curr_mesh_triangle.neighbours.push_back(&neighb_mesh_triangle);
				}
			}
		}
	}
}

void MeshKeeper::Clear()
{
	m_vertices.clear();
	m_triangles.clear();
	m_triangles_with_indexes.clear();
	m_mesh_vertices.clear();
	m_mesh_triangles.clear();
	m_triangle_neighb_of_vertex.clear();
	m_neighbours_numbers.clear();
}

}