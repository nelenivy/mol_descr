#pragma once
#ifndef SURFACE_READER_H
#define SURFACE_READER_H

#include <fstream>
#include <string>
#include <vector>
#include <fstream>
#include <opencv2/core/core.hpp>

namespace molecule_descriptor
{

class SurfaceReader
{
public:
	bool OpenFile(const std::string& file_name);
	void ReadVertices(std::vector<cv::Point3d>& vertices);
	void ReadNormals(std::vector<cv::Point3d>& normals);
	void ReadTriangles(std::vector<cv::Point3i>& triangles);
private:
	template<typename T>
	void ReadSurface(const std::string& points_mark, const std::string& points_start, std::vector<cv::Point3_<T>>& points);
private:
	std::ifstream m_file_in;
	std::string m_string_buf;
};

}

#endif