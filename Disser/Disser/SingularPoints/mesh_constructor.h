#pragma once
#ifndef MESH_CONSTRUCTOR_H
#define MESH_CONSTRUCTOR_H

#include <vector>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"
#include "opencv2/core/core.hpp"
#include "mesh_types.h"

namespace molecule_descriptor
{

//class for constructing and keeping graph
class MeshKeeper
{
public:
	MeshKeeper() { }
	void ConstructMesh(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, const std::vector<cv::Point3i>& triangles);
	const VerticesGraph& GetMeshVertices() {return m_vertices_graph; }
	const TrianglesGraph& GetMeshTriangles() {return m_triangles_graph; }
private:
	MeshKeeper(const MeshKeeper&){};
	MeshKeeper& operator=(const MeshKeeper&){};
	void Clear();//call in ConstructMesh
	void FillVertices(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals);
	void FillTriangles(const std::vector<cv::Point3i>& triangles);
	void FindTriangleNeighboursOfEachVertex();
	void FillMeshVertices();
	void FillMeshTriangles();
private:
	VerticesGraph m_vertices_graph;
	TrianglesGraph m_triangles_graph;
	std::vector<Vertice>	m_vertices;
	std::vector<cv::Point3i>	m_triangles_with_indexes;

	std::vector<std::vector<int>>	m_triangle_neighb_of_vertex;
	std::vector<int>			m_neighbours_numbers;
};

}

#endif