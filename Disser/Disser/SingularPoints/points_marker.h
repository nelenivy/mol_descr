#pragma once
#include <vector>
#include <utility>

#include "opencv2/core/core.hpp"

#include "CommonUtilities/attributes_container.h"
#include "GraphLib\graph_structures.h"
#include "GraphLib\segmentation_types.h"
#include "singular_point_types.h"
#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

template <typename NodeType>
void CalculatePotential(const std::vector<NodeType>& points, const std::vector<std::pair<cv::Point3d, double>>& charges);


template <typename NodeType>
void CalculatePotential(const std::vector<NodeType>& points, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	for (auto points_iter = points.begin(); points_iter != points.end(); ++points_iter)
	{
		points_iter->attr.Add<ElectricSign>();
		points_iter->attr.Add<ElectricPotential>();
		points_iter->attr.Get<ElectricPotential>() = 0.0;

		for (auto charges_iter = charges.begin(); charges_iter != charges.end(); ++charges_iter)
		{
			const double curr_dist = cv::norm(points_iter->Center() - charges_iter->first);
			points_iter->attr.Get<ElectricPotential>() += charges_iter->second / curr_dist;
		}

		points_iter->attr.Get<ElectricSign>() = Sign(points_iter->attr.Get<ElectricPotential>());
	}
}
}