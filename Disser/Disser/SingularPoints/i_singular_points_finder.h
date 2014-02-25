#pragma once 

#include <utility>
#include <vector>
#include <memory>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

//abstract interface to hide implementation details
class ISingularPointsFinder
{
public:
	virtual void Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges) = 0;
	//returns number of types of properties. Is a product of properties ranges
	virtual size_t GetTypesNum() = 0;
	virtual void GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points) = 0;
	virtual void GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& non_marked_singular_points) = 0;
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers) = 0;
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types) = 0;
	static const size_t kHistSize = 9;
	virtual void GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist) = 0;

	virtual void Release() = 0;
};

std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder();

}