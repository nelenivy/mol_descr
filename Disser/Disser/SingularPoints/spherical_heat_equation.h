#pragma once

#include <vector>
#include <utility>
#include <limits>
#include "opencv2/core/core.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"
#include "SingularPoints/mesh_types.h"
#include "CommonUtilities/common_functions.h"


namespace molecule_descriptor
{

class SphericalHeatEquation
{
public:

	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> VecDoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Point3d>, VERTEX> SphereMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> HessianMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> CoordBasisMap;

	template <class Graph, class CoordMap, class DistanceMap,class TangentBasisMap>
	void ProcessIteration(const Graph& graph, const CoordMap& coord_map, 
		const SphereMap& prev_iter_func, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map, const double delta_t, SphereMap& curr_iter_func)
	{
		curr_iter_func.SetGraph(graph);
		m_hessian_calculator.Process(graph, coord_map, prev_iter_func, dist_map, dist_thresh, tangent_basis_map, 8);

		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{
			//in normal coordinates x and y of previous solution are equal to zero, and trace is delta
			const double x_sphere = delta_t * (m_hessian_calculator.m_hessian_map_x[*curr_vert](0 ,0) +  m_hessian_calculator.m_hessian_map_x[*curr_vert](1 ,1));
			const double y_sphere = delta_t * (m_hessian_calculator.m_hessian_map_y[*curr_vert](0 ,0) +  m_hessian_calculator.m_hessian_map_y[*curr_vert](1 ,1));
			//std::cout << "(" <<x_sphere << "," << y_sphere << "),";
			const double phi_sphere = sqrt(Sqr(x_sphere) + Sqr(y_sphere));
			const double theta_sphere = atan2(y_sphere, x_sphere);
			cv::Mat_<double> new_coord(3, 1);
			new_coord(0) = cos(theta_sphere) * sin(phi_sphere);
			new_coord(1) = sin(theta_sphere) * sin(phi_sphere);
			new_coord(2) = cos(phi_sphere);
			//std::cout << new_coord;
			const cv::Mat_<double>& curr_coord_sys = m_hessian_calculator.m_coord_sphere_basis_map[*curr_vert];
			//std::cout << curr_coord_sys << curr_coord_sys.inv(DECOMP_SVD);
			cv::Mat_<double> old_coord = curr_coord_sys * new_coord;
			curr_iter_func[*curr_vert].x = old_coord(0);
			curr_iter_func[*curr_vert].y = old_coord(1);
			curr_iter_func[*curr_vert].z = old_coord(1);
			//std::cout << old_coord;
		}
	}

	HessianMatrixCalculatorSpherical m_hessian_calculator;
};

}