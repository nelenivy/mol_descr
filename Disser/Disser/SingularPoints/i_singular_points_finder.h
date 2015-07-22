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
		FIRST_SURF_PROP = FIRST_PROP, 
		MEAN_CURV, ELECTR_POTENT, STERIC_POTENT, 
		SURF_PROPS_NUM,//properties which are used in singular points calculation

		LOG_GAUSS_CURV,
		FIRST_LOG_PROP = LOG_GAUSS_CURV, 
		LOG_MEAN_CURV, LOG_ELECTR_POTENT, LOG_STERIC_POTENT,
		LAST_LOG_PROP = LOG_STERIC_POTENT, 

		PCA_GAUSS_CURV,
		FIRST_PCA_PROP = PCA_GAUSS_CURV, 
		PCA_MEAN_CURV, PCA_ELECTR_POTENT, PCA_STERIC_POTENT,
		LAST_PCA_PROP = PCA_STERIC_POTENT, 

		SEGM_NUM = LAST_PCA_PROP + 1,
		FIRST_ADDITIONAL_PROP = SEGM_NUM, 
		SURF_TYPE, PCA_LOG, PCA_SCALE_SPACE, 		
		PROPS_NUM
	};
	virtual void GetVerticesWithDblProp(std::vector<std::pair<cv::Point3d, double>>& vertices_with_prop, 
		const ISingularPointsFinder::SurfProperty prop_type) {}
	virtual void GetVerticesWithDblPropLevels(std::vector<std::vector<std::pair<cv::Point3d, double>>>& vertices_with_props_lev, 
		const ISingularPointsFinder::SurfProperty prop_type) {}
	virtual void AppendProp(std::vector<double>& prop, SurfProperty prop_type)
	{
	}
	virtual void SetMeanAndSigma(std::vector<std::vector<double>> mean_and_sigma)
	{
	}

	virtual void InitParams(int argc, char** argv) = 0;
	virtual ~ISingularPointsFinder() { }
};

inline const char* SurfPropertyName(const ISingularPointsFinder::SurfProperty prop_type)
{
	if (prop_type == ISingularPointsFinder::GAUSS_CURV)
	{
		return "GAUSS_CURV";
	}
	else if (prop_type == ISingularPointsFinder::MEAN_CURV)
	{
		return "MEAN_CURV";
	}
	else if (prop_type == ISingularPointsFinder::ELECTR_POTENT)
	{
		return "ELECTR_POTENT";
	}
	else if (prop_type == ISingularPointsFinder::STERIC_POTENT)
	{
		return "STERIC_POTENT";
	}
	else 	if (prop_type == ISingularPointsFinder::LOG_GAUSS_CURV)
	{
		return "LOG_GAUSS_CURV";
	}
	else if (prop_type == ISingularPointsFinder::LOG_MEAN_CURV)
	{
		return "LOG_MEAN_CURV";
	}
	else if (prop_type == ISingularPointsFinder::LOG_ELECTR_POTENT)
	{
		return "LOG_ELECTR_POTENT";
	}
	else if (prop_type == ISingularPointsFinder::LOG_STERIC_POTENT)
	{
		return "LOG_STERIC_POTENT";
	}
	else 	if (prop_type == ISingularPointsFinder::PCA_GAUSS_CURV)
	{
		return "PCA_GAUSS_CURV";
	}
	else if (prop_type == ISingularPointsFinder::PCA_MEAN_CURV)
	{
		return "PCA_MEAN_CURV";
	}
	else if (prop_type == ISingularPointsFinder::PCA_ELECTR_POTENT)
	{
		return "PCA_ELECTR_POTENT";
	}
	else if (prop_type == ISingularPointsFinder::PCA_STERIC_POTENT)
	{
		return "PCA_STERIC_POTENT";
	}
	else if (prop_type == ISingularPointsFinder::SEGM_NUM)
	{
		return "SEGM_NUM";
	}
	else if (prop_type == ISingularPointsFinder::SURF_TYPE)
	{
		return "SURF_TYPE";
	}
	else if (prop_type == ISingularPointsFinder::PCA_LOG)
	{
		return "PCA_LOG";
	}
	else if (prop_type == ISingularPointsFinder::PCA_SCALE_SPACE)
	{
		return "PCA_SCALE_SPACE";
	}
	else 
	{
		CV_Assert(0, "Unknown srface type");
	}		
}

enum SingularPointsAlgorithm
{
	SEGMENTATION, SCALE_SPACE
};
std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder(const SingularPointsAlgorithm alg);

}