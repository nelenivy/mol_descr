#include "surface_reader.h"
#include <sstream>

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

bool SurfaceReader::OpenFile(const string& file_name)
{
	if (m_file_in.is_open())
	{
		m_file_in.close();
	}
	m_file_in.clear();
	m_file_in.open(file_name);
	m_string_buf.reserve(500);
	return m_file_in.is_open();
}

void SurfaceReader::ReadVertices(vector<cv::Point3d>& vertices)
{
	const string kVerticesMark = "Triangle";
	const string kVerticesStart = "point";
	ReadSurface(kVerticesMark, kVerticesStart, vertices);
}

void SurfaceReader::ReadNormals(vector<cv::Point3d>& normals)
{
	const string kNormalsMark = "# Normal definition";
	const string kNormalsStart = "vector";
	ReadSurface(kNormalsMark, kNormalsStart, normals);
}

void SurfaceReader::ReadTriangles(vector<Point3i>& triangles)
{
	const string kTrianglesMark = "IndexedFaceSet";
	const string kTrianglesStart = "coordIndex";
	ReadSurface(kTrianglesMark, kTrianglesStart, triangles);
}
//functions for reading point from line
template<typename T>
bool ExtractNumberFromLine(std::stringstream& line_stream, T& num)
{
	line_stream >> num;
	int iterations_num = 0;
	//try to read number while skipping
	while (line_stream.fail() == true	&&
		line_stream.eof()	== false	&& 
		iterations_num < 1000)
	{
		line_stream.clear();
		char c = 0;
		line_stream.get(c);
		line_stream >> num;
		iterations_num++;
	}

	if (line_stream.fail())
	{
		return false;
	}
	else
	{
		return true;
	}
}

template<typename T>
bool ExtractPointFromLine(std::string& curr_line, cv::Point3_<T>& point)
{
	std::stringstream line_stream;
	line_stream << curr_line;

	bool res = ExtractNumberFromLine(line_stream, point.x) && 
	ExtractNumberFromLine(line_stream, point.y) &&
	ExtractNumberFromLine(line_stream, point.z);

	return res;
}
//////////////////////////////////////////////////////////////////////////
template<typename T>
void SurfaceReader::ReadSurface(const string& points_mark, const string& points_start, vector<Point3_<T>>& points)
{
	CV_Assert(m_file_in.is_open() && ! m_file_in.fail());

	while (m_string_buf.find(points_mark) == string::npos && ! m_file_in.eof()) 
	{
		std::getline(m_file_in, m_string_buf);
	}

	while (m_string_buf.find(points_start) == string::npos && !m_file_in.eof())
	{
		std::getline(m_file_in, m_string_buf);
	}

	//read vertices
	size_t vertices_size = 0;      
	cv::Point3_<T> new_point;
	points.clear();

	while(1)
	{			   
		std::getline(m_file_in, m_string_buf);

		if (m_file_in.rdstate() == 0 && ExtractPointFromLine(m_string_buf, new_point))
		{
			points.push_back(new_point);
		}
		else
		{
			break;
		}
	}
	
	m_file_in.clear();
	m_file_in.seekg(0, ios_base::beg);
}

}