#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "mesh_constructor.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

class SingularPointsFinder : public ISingularPointsFinder
{
public:
	SingularPointsFinder() { }
	virtual void Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges);
	virtual void GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points);
	virtual void GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& non_marked_singular_points);
	//TODO: fix compile-time properties num. Maybe put pairs into this class
	virtual size_t GetTypesNum();
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers);
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	virtual void Release() { delete this; }
private:
	void Clear();//call before processing new data
	void CalculateVerticesSurfaceType();
	void SegmentMolecularSurface(const size_t max_segm_size);
	void FindSegmentsCenters();
	void CalculateSingularPointsTypes();

	MeshKeeper m_mesh_keeper;
	std::vector<Vertice> m_singular_points;
};

}