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
	virtual void Process(const std::vector<cv::Point3d>& vertices, 
		const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, 
		const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
		const bool calc_prop_as_average) = 0;
	virtual void GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, 
		std::vector<size_t>> & non_marked_singular_points) = 0;
	virtual void GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, 
		std::vector<std::vector<size_t>>>& non_marked_singular_points) = 0;
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers) = 0;
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types) = 0;
	virtual void GetVerticesWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types) {};

	virtual void CalcOnlyProps(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii)
	{

	}
	enum SurfProperty
	{
		GAUSS_CURV = 0, 
		FIRST_PROP = GAUSS_CURV, 
		MEAN_CURV, ELECTR_POTENT, STERIC_POTENT, 
		PROPS_NUM
	};
	virtual void AppendProp(std::vector<double>& prop, SurfProperty prop_type)
	{
	}
	virtual void SetMeanAndSigma(std::vector<std::vector<double>> mean_and_sigma)
	{
	}

	virtual void InitParams(int argc, char** argv) = 0;
	virtual ~ISingularPointsFinder() { }
};

enum SingularPointsAlgorithm
{
	SEGMENTATION, SCALE_SPACE
};
std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder(const SingularPointsAlgorithm alg);

}