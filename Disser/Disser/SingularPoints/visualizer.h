#pragma once

#include <vector>
#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include "mesh_types.h"

namespace molecule_descriptor
{

template <typename SphereMap>
void VisualizeSphereMap(const std::vector<Vertice>& vertices, const std::vector<cv::Point3i>& triangles, const SphereMap& prop_map)
{
	Eigen::MatrixXd V(vertices.size(), 3);
	Eigen::MatrixXi F(triangles.size(), 3);

	for (int i = 0; i < vertices.size(); ++i)
	{
		V(i, 0) = vertices[i].Center().x;
		V(i, 1) = vertices[i].Center().y;
		V(i, 2) = vertices[i].Center().z;
	}

	for (int i = 0; i < triangles.size(); ++i)
	{
		F(i, 0) = triangles[i].x;
		F(i, 1) = triangles[i].y;
		F(i, 2) = triangles[i].z;
	}

	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	Eigen::MatrixXd P1(vertices.size(), 3);
	Eigen::MatrixXd P2(vertices.size(), 3);

	for (int i = 0; i < vertices.size(); ++i)
	{
		P1(i, 0) = vertices[i].Center().x;
		P1(i, 1) = vertices[i].Center().y;
		P1(i, 2) = vertices[i].Center().z;

		P2(i, 0) = vertices[i].Center().x + prop_map[i].x;
		P2(i, 1) = vertices[i].Center().y + prop_map[i].y;
		P2(i, 2) = vertices[i].Center().z + prop_map[i].z;
	}
	viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,0,0));
	viewer.launch();
}

template <typename SphereMap>
void VisualizeSphereMapWithSingPts(const std::vector<Vertice>& vertices, const std::vector<cv::Point3i>& triangles, const SphereMap& prop_map, 
								   const std::vector<cv::Point3d>& sing_pts)
{
	Eigen::MatrixXd V(vertices.size(), 3);
	Eigen::MatrixXi F(triangles.size(), 3);

	for (int i = 0; i < vertices.size(); ++i)
	{
		V(i, 0) = vertices[i].Center().x;
		V(i, 1) = vertices[i].Center().y;
		V(i, 2) = vertices[i].Center().z;
	}

	for (int i = 0; i < triangles.size(); ++i)
	{
		F(i, 0) = triangles[i].x;
		F(i, 1) = triangles[i].y;
		F(i, 2) = triangles[i].z;
	}
	std::cout << 1 << std::endl;
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	std::cout << 2 << std::endl;
	Eigen::MatrixXd P1(vertices.size(), 3);
	Eigen::MatrixXd P2(vertices.size(), 3);

	for (int i = 0; i < vertices.size(); ++i)
	{
		P1(i, 0) = vertices[i].Center().x;
		P1(i, 1) = vertices[i].Center().y;
		P1(i, 2) = vertices[i].Center().z;

		P2(i, 0) = vertices[i].Center().x + prop_map[i].x;
		P2(i, 1) = vertices[i].Center().y + prop_map[i].y;
		P2(i, 2) = vertices[i].Center().z + prop_map[i].z;
	}
	//viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,0,0));
	std::cout << 3 << std::endl;
	Eigen::MatrixXd P(sing_pts.size(), 3);

	for (int i = 0; i < sing_pts.size(); ++i)
	{
		P(i, 0) = sing_pts[i].x;
		P(i, 1) = sing_pts[i].y;
		P(i, 2) = sing_pts[i].z;
	}

	//viewer.data().add_points(P,Eigen::RowVector3d(1,0,0));
	std::cout << 4 << std::endl;
	viewer.launch();
	std::cout << 5 << std::endl;
}

}