#include "i_singular_points_finder.h"
#include <memory>
#include <vector>
#include <utility>
#include <iostream>

#include "opencv2/core/core.hpp"

void SingularPointsTest()
{
	auto instance = molecule_descriptor::CreateSingularPointsFinder();
	std::vector<cv::Point3i> triangles(2);
	triangles[0] = cv::Point3i(0, 1, 2);
	triangles[1] = cv::Point3i(3, 4, 5);
	std::vector<cv::Point3d> vertices(6);
	vertices[0] = cv::Point3d(0.0, 0.0, 0.0);
	vertices[1] = cv::Point3d(1.0, 0.0, 0.0);
	vertices[2] = cv::Point3d(1.0, 1.0, 0.0);
	vertices[3] = cv::Point3d(1.0, 1.0, 1.0);
	vertices[4] = cv::Point3d(1.0, 0.0, 1.0);
	vertices[5] = cv::Point3d(0.0, 0.0, 1.0);

	std::vector<cv::Point3d> normals(6);
	normals[0] = cv::Point3d(1.0, 0.0, 0.0);
	normals[1] = cv::Point3d(1.0, 0.0, 0.0);
	normals[2] = cv::Point3d(1.0, 1.0, 0.0);
	normals[3] = cv::Point3d(1.0, 1.0, 1.0);
	normals[4] = cv::Point3d(1.0, 0.0, 1.0);
	normals[5] = cv::Point3d(0.0, 0.0, 1.0);
	std::vector<std::pair<cv::Point3d, double>> charges(10);
	instance->Process(vertices, normals, triangles, charges);
	std::vector<std::pair<cv::Point3d, size_t>> segmented_vertices(10);

	instance->GetSegmentedVertices(segmented_vertices);
	//instance->GetMarkedSingularPoints(segmented_vertices);
	instance->GetVerticesWithTypes(segmented_vertices);
}