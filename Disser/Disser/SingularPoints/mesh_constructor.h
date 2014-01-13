#pragma once
#ifndef MESH_CONSTRUCTOR_H
#define MESH_CONSTRUCTOR_H

#include <vector>
#include "opencv2/core/core.hpp"
#include "GraphLib\graph_structures.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

//class for constructing and keeping graph
class MeshKeeper
{
public:
	MeshKeeper() { }
	void ConstructMesh(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, const std::vector<cv::Point3i>& triangles);
	const std::vector<MeshVertice>* GetMeshVertices() {return &m_mesh_vertices; }
	const std::vector<MeshTriangle>* GetMeshTriangles() {return &m_mesh_triangles; }
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
	std::vector<Vertice>		m_vertices;
	std::vector<Triangle>	m_triangles;
	std::vector<cv::Point3i>		m_triangles_with_indexes;
	std::vector<MeshVertice> m_mesh_vertices;
	std::vector<MeshTriangle> m_mesh_triangles;

	std::vector<std::vector<int>>	m_triangle_neighb_of_vertex;
	std::vector<int>			m_neighbours_numbers;
};

}

#endif