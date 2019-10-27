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

template <class Graph, class CoordMap, class FunctionMap, class DistanceMap>
void CalcDifferentialInVertLeastSquares(const Graph& graph, const CoordMap& coord_map, 
						const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
						const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,
						const cv::Mat_<double>& tangent_basis, const double min_grad, 
						const double max_grad, cv::Mat_<double>& grad)//vectors in columns, 3rd vector is normal to surface
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	std::vector<VertexDescriptor> vertices_within_dist;
	GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
	cv::Mat_<double> cent_coord;
	Point_ToMat_Transposed(coord_map[cent_vert], cent_coord);
	//FILL EQUATION FOR GRADIENT CALCULATION
	cv::Mat_<double> vect_in_local_coords, neighb_coord;

	cv::Mat_<double> coeffs(vertices_within_dist.size(), 2);
	cv::Mat_<double> derivatives(vertices_within_dist.size(), 1);

	for (auto neighb_vert = vertices_within_dist.begin(), end_vert = vertices_within_dist.end(); neighb_vert != end_vert; ++neighb_vert)
	{
		const size_t curr_row = neighb_vert - vertices_within_dist.begin();
		if (*neighb_vert == cent_vert)
		{
			coeffs(curr_row, 0) = coeffs(curr_row, 1) = derivatives(curr_row, 0) = 0.0;
			continue;
		}
		Point_ToMat_Transposed(coord_map[*neighb_vert], neighb_coord);
		vect_in_local_coords = tangent_basis * (neighb_coord - cent_coord);
		vect_in_local_coords(2) = 0.0; //vector should belong to tangent plane, so z coordinate should be 0
		const double curve_length = cv::norm(neighb_coord - cent_coord);
		vect_in_local_coords *= curve_length / cv::norm(vect_in_local_coords);
		coeffs(curr_row, 0) = vect_in_local_coords(0);
		coeffs(curr_row, 1) = vect_in_local_coords(1);

		//const double curve_length = sqrt(Sqr(projected_vect(0)) + Sqr(projected_vect(1)));
		
		derivatives(curr_row, 0) = (func_map[*neighb_vert] - func_map[cent_vert]) / curve_length;
	}

	//SOLVE EQUATION
	grad = coeffs.inv(DECOMP_SVD) * derivatives;
	grad(0) = std::min(max_grad, std::max(min_grad, grad(0)));
	grad(1) = std::min(max_grad, std::max(min_grad, grad(1)));
}

template <class Graph, class CoordMap, class FunctionMap, class DistanceMap,class TangentBasisMap, class GradientMap>
void CalcDifferentialMaps(const Graph& graph, const CoordMap& coord_map, 
					 const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
					 const TangentBasisMap& tangent_basis_map, const bool least_squares_grad, 
					 GradientMap& grad_x, GradientMap& grad_y)
{
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		cv::Mat_<double> curr_grad;
		CalcDifferentialInVertLeastSquares(graph, coord_map, func_map, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], -1000.0, 1000.0, curr_grad);
		grad_x[*curr_vert] = curr_grad(0);
		grad_y[*curr_vert] = curr_grad(1);
	}
}

template <class Graph, class CoordMap, class FunctionMap, class DistanceMap>
void CalcHessianByTaylorInVertLeastSquares(const Graph& graph, const CoordMap& coord_map, 
											const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
											const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,
											const cv::Mat_<double>& tangent_basis, //vectors in columns, 3rd vector is normal to surface
											const cv::Mat_<double>& differential, const double min_grad, 
											const double max_grad, cv::Mat_<double>& hessian)
{
	//The realization through the formula H(p-p0, p-p0) = f(p) - (f(p_0) + df(p-p0)) + o(p-p_0)^2
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	std::vector<VertexDescriptor> vertices_within_dist;
	GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
	cv::Mat_<double> cent_coord;
	Point_ToMat_Transposed(coord_map[cent_vert], cent_coord);
	//FILL EQUATION FOR GRADIENT CALCULATION
	cv::Mat_<double> projected_vect, neighb_coord;
	cv::Mat_<double> vect_in_local_coords = cv::Mat_<double>::zeros(1, 2);
	cv::Mat_<double> coeffs = cv::Mat_<double>::zeros(vertices_within_dist.size(), 3);
	cv::Mat_<double> func_diff(vertices_within_dist.size(), 1);

	for (auto neighb_vert = vertices_within_dist.begin(), end_vert = vertices_within_dist.end(); neighb_vert != end_vert; ++neighb_vert)
	{
		const size_t curr_row = neighb_vert - vertices_within_dist.begin();
		if (*neighb_vert == cent_vert)
		{
			coeffs(curr_row, 0) = coeffs(curr_row, 1) = coeffs(curr_row, 2) = func_diff(curr_row, 0) = 0.0;
			continue;
		}
		Point_ToMat_Transposed(coord_map[*neighb_vert], neighb_coord);
		projected_vect = tangent_basis * (neighb_coord - cent_coord);
		projected_vect(2) = 0;
		const double curve_length = cv::norm(neighb_coord - cent_coord);
		projected_vect *= curve_length / cv::norm(projected_vect);
		vect_in_local_coords(0) = projected_vect(0);
		vect_in_local_coords(1) = projected_vect(1);
		const cv::Mat_<double> f_lin = func_map[cent_vert] + vect_in_local_coords * differential;//f_lin(p) = f(p_0)  + df(p-p_0)
		CV_Assert(f_lin.size().area() == 1);
		const double funct_first_ord_approx = f_lin(0);

		coeffs(curr_row, 0) = vect_in_local_coords(0) * vect_in_local_coords(0);
		coeffs(curr_row, 1) = 2.0 * vect_in_local_coords(0) * vect_in_local_coords(1);
		coeffs(curr_row, 2) = vect_in_local_coords(1) * vect_in_local_coords(1);		
		func_diff(curr_row, 0) = func_map[*neighb_vert] - funct_first_ord_approx;
	}

	//SOLVE EQUATION
	cv::Mat_<double> hess_vec = coeffs.inv(DECOMP_SVD) * func_diff;
	hess_vec(0) = std::min(max_grad, std::max(min_grad, hess_vec(0)));
	hess_vec(1) = std::min(max_grad, std::max(min_grad, hess_vec(1)));
	hess_vec(2) = std::min(max_grad, std::max(min_grad, hess_vec(2)));
	hessian.create(2, 2);
	hessian(0, 0) = hess_vec(0);
	hessian(1, 0) = hess_vec(1);
	hessian(0, 1) = hess_vec(1);
	hessian(1, 1) = hess_vec(2);
}


template <class Graph, class CoordMap, class FunctionMap, class DistanceMap>
void CalcCovariantDifferentialInVertLeastSquares(const Graph& graph, const CoordMap& coord_map, 
										const FunctionMap& func_vec_map, const DistanceMap& dist_map, const double dist_thresh, 
										const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,
										const cv::Mat_<double>& tangent_basis, const double min_grad, 
										const double max_grad, cv::Mat_<double>& grad_x, cv::Mat_<double>& grad_y)//vectors in columns, 3rd vector is normal to surface
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	std::vector<VertexDescriptor> vertices_within_dist;
	GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
	cv::Mat_<double> cent_coord;
	Point_ToMat_Transposed(coord_map[cent_vert], cent_coord);
	//FILL EQUATION FOR GRADIENT CALCULATION
	cv::Mat_<double> projected_vect, neighb_coord;

	cv::Mat_<double> coeffs = cv::Mat_<double>::zeros(vertices_within_dist.size(), 2);
	cv::Mat_<double> derivatives_x(vertices_within_dist.size(), 1);
	cv::Mat_<double> derivatives_y(vertices_within_dist.size(), 1);

	for (auto neighb_vert = vertices_within_dist.begin(), end_vert = vertices_within_dist.end(); neighb_vert != end_vert; ++neighb_vert)
	{
		const size_t curr_row = neighb_vert - vertices_within_dist.begin();
		if (*neighb_vert == cent_vert)
		{
			coeffs(curr_row, 0) = coeffs(curr_row, 1) = derivatives_x(curr_row, 0) = derivatives_y(curr_row, 0) = 0.0;
			continue;
		}
		Point_ToMat_Transposed(coord_map[*neighb_vert], neighb_coord);
		projected_vect = tangent_basis * (neighb_coord - cent_coord);
		projected_vect(2) = 0;
		const double curve_length = cv::norm(neighb_coord - cent_coord);
		projected_vect *= curve_length / cv::norm(projected_vect);
		coeffs(curr_row, 0) = projected_vect(0);
		coeffs(curr_row, 1) = projected_vect(1);

		//const double curve_length = sqrt(Sqr(projected_vect(0)) + Sqr(projected_vect(1)));
		const double curve_length = cv::norm(neighb_coord - cent_coord);
		const cv::Mat_<double> vec_diff = func_vec_map[*neighb_vert] - func_vec_map[cent_vert];
		const cv::Mat_<double> proj_vec_diff = tangent_basis * vec_diff;
		derivatives_x(curr_row, 0) = proj_vec_diff(0) / curve_length;
		derivatives_y(curr_row, 0) = proj_vec_diff(1) / curve_length;
	}

	//SOLVE EQUATION
	grad_x = coeffs.inv(DECOMP_SVD) * derivatives_x;
	grad_y = coeffs.inv(DECOMP_SVD) * derivatives_y;
	grad_x(0) = std::min(max_grad, std::max(min_grad, grad_x(0)));
	grad_x(1) = std::min(max_grad, std::max(min_grad, grad_x(1)));

	grad_y(0) = std::min(max_grad, std::max(min_grad, grad_y(0)));
	grad_y(1) = std::min(max_grad, std::max(min_grad, grad_y(1)));
}

template <class Graph, class CoordMap, class FunctionVecMap, class DistanceMap,class TangentBasisMap, class CovDiffMap>
void CalcCovariantDifferentialMaps(const Graph& graph, const CoordMap& coord_map, 
						  const FunctionVecMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
						  const TangentBasisMap& tangent_basis_map, 
						  CovDiffMap& cov_diff)
{
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{		
		cv::Mat_<double> curr_grad_x, curr_grad_y;
		CalcCovariantDifferentialInVertLeastSquares(graph, coord_map, func_map, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert],
			-1000.0, 1000.0, curr_grad_x, curr_grad_y);
		
		cv::Mat_<double> curr_grad(2, 2);
		curr_grad(0,0) = curr_grad_x(0);
		curr_grad(0,1) = (curr_grad_x(1) + curr_grad_y(0)) / 2.0;
		curr_grad(1,0) = (curr_grad_x(1) + curr_grad_y(0)) / 2.0;
		curr_grad(1,1) = curr_grad_y(1);
		cov_diff[*curr_vert] = curr_grad;
	}
}

class HessianMatrixCalculatorSpherical
{
public:
	
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> VecDoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Point3d>, VERTEX> SphereMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> HessianMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> CoordBasisMap;

	template <class Graph, class CoordMap, class DistanceMap,class TangentBasisMap>
	void Process(const Graph& graph, const CoordMap& coord_map, 
		const SphereMap& func_map, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map)
	{
		cv::Mat_<double> sphere_point;
		m_spher_normal_coord_x.SetGraph(graph);
		m_spher_normal_coord_y.SetGraph(graph);
		m_hessian_map_x.SetGraph(graph);
		m_hessian_map_y.SetGraph(graph);
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{
			ConvertToLocalNormalCoords(graph, func_map, dist_map, *curr_vert,dist_thresh);
			cv::Mat_<double> curr_grad;
			//Calculate Hessian for each coordinate (x, y) of local coordinates
			CalcDifferentialInVertLeastSquares(graph, coord_map, m_spher_normal_coord_x, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], 
				-1000.0, 1000.0, curr_grad);
			CalcHessianByTaylorInVertLeastSquares(graph, coord_map, m_spher_normal_coord_x, dist_map, dist_thresh, 
					*curr_vert,	tangent_basis_map[*curr_vert], //vectors in columns, 3rd vector is normal to surface
					curr_grad, -1000.0, 1000.0, m_hessian_map_x[*curr_vert]);

			CalcDifferentialInVertLeastSquares(graph, coord_map, m_spher_normal_coord_y, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], 
				-1000.0, 1000.0, curr_grad);
			CalcHessianByTaylorInVertLeastSquares(graph, coord_map, m_spher_normal_coord_y, dist_map, dist_thresh, 
				*curr_vert,	tangent_basis_map[*curr_vert], //vectors in columns, 3rd vector is normal to surface
				curr_grad, -1000.0, 1000.0, m_hessian_map_y[*curr_vert]);
		}
	}

	template <class Graph, class DistanceMap>
	void ConvertToLocalNormalCoords(const Graph& graph, 
		const SphereMap& func_map, const DistanceMap& dist_map, 
		const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,const double dist_thresh)
	{
		//Find local coordinates
		const cv::Point3d& curr_sphere_point = func_map[cent_vert];
		const double theta = atan(curr_sphere_point.y / (curr_sphere_point.x + 0.0000000001));
		const double phi = acos(curr_sphere_point.z);
		const cv::Point3d r_vect = cv::Point3d(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
		const cv::Point3d theta_vect = cv::Point3d(-sin(theta), cos(theta), 0);
		const cv::Point3d phi_vect = cv::Point3d(cos(theta) * cos(phi), sin(theta) * cos(phi), -sin(phi));
		cv::Mat_<double> new_basis(3, 3), sphere_point;
		Point_ToMat_Transposed(r_vect, new_basis.col(0));
		Point_ToMat_Transposed(theta_vect, new_basis.col(1));
		Point_ToMat_Transposed(r_vect, new_basis.col(2));
		//Fill x and  coords in neighbouring vertices
		std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vertices_within_dist;
		GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
		for (auto neighb_vert = vertices_within_dist.begin(), end_vert = vertices_within_dist.end(); neighb_vert != end_vert; ++neighb_vert)
		{
			if (*neighb_vert == cent_vert)
			{
				m_spher_normal_coord_x[*neighb_vert] = m_spher_normal_coord_y[*neighb_vert] = 0.0;
				continue;
			}
			Point_ToMat_Transposed(func_map[*neighb_vert], sphere_point);
			cv::Mat_<double> sphere_point_in_new_coords = new_basis * sphere_point;
			const double new_phi = acos(sphere_point_in_new_coords(2));
			const double new_theta = atan(sphere_point_in_new_coords(1) / (sphere_point_in_new_coords(0) + 0.0000000001));
			m_spher_normal_coord_x[*neighb_vert] = new_phi * cos(new_theta);
			m_spher_normal_coord_y[*neighb_vert] = new_phi * sin(new_theta);
		}
	}

	HessianMap m_hessian_map_x, m_hessian_map_y;
	DoublePropMap m_spher_normal_coord_x;
	DoublePropMap m_spher_normal_coord_y;
};


class HessianMatrixCalculator
{
public:
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> VecDoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> HessianMap;
	
	template <class Graph, class CoordMap, class FunctionMap, class DistanceMap,class TangentBasisMap, class RatioMap>
	void Process(const Graph& graph, const CoordMap& coord_map, 
		const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map, const bool least_squares_grad,RatioMap& eigenvalues_ratio)
	{
		m_dx.SetGraph(graph);
		m_dy.SetGraph(graph);		
		CalcDifferentialMaps(graph, coord_map, func_map,dist_map, dist_thresh, tangent_basis_map, least_squares_grad, m_dx, m_dy);
		
		ProcessFromGradient(graph, coord_map, func_map,m_dx, m_dy, dist_map, dist_thresh,tangent_basis_map, least_squares_grad, eigenvalues_ratio);
	}
	
	template <class Graph, class CoordMap, class FuncMap, class DistanceMap,class TangentBasisMap>
	void CalcHessianTaylor(const Graph& graph, const CoordMap& coord_map, 
		const FuncMap& func_map, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map)
	{		
		m_grad_2d.SetGraph(graph);
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{	
			m_grad_2d[*curr_vert].create(2, 1);
			m_grad_2d[*curr_vert](0, 0) = m_dx[*curr_vert];
			m_grad_2d[*curr_vert](1, 0) = m_dy[*curr_vert];			
		}
		m_hessian_map.SetGraph(graph);

		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{		
			CalcHessianByTaylorInVertLeastSquares(graph, coord_map, 
				func_map, dist_map, dist_thresh, 
				*curr_vert,
				tangent_basis_map[*curr_vert], //vectors in columns, 3rd vector is normal to surface
				m_grad_2d[*curr_vert], -1000.0, 1000.0, m_hessian_map[*curr_vert]);
		}
	}

	template <class Graph, class CoordMap, class GradMap, class DistanceMap,class TangentBasisMap>
	void CalcHessianCovDiff(const Graph& graph, const CoordMap& coord_map, 
		const GradMap& dx, const GradMap& dy, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map)
	{
		cv::Mat_<double> grad_2d(3, 1);
		m_grad_3d.SetGraph(graph);
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{		
			grad_2d(0, 0) = dx[*curr_vert];
			grad_2d(1, 0) = dy[*curr_vert];
			grad_2d(2, 0) = 0;
			m_grad_3d[*curr_vert] = tangent_basis_map[*curr_vert].t() * grad_2d;			
		}
		m_hessian_map.SetGraph(graph);
		CalcCovariantDifferentialMaps(graph, coord_map, m_grad_3d,dist_map, dist_thresh, tangent_basis_map, m_hessian_map);	
	}

	template <class Graph, class CoordMap, class FuncMap, class GradMap, class DistanceMap,class TangentBasisMap, class RatioMap>
	void ProcessFromGradient(const Graph& graph, const CoordMap& coord_map, const FuncMap& func_map,
		const GradMap& dx, const GradMap& dy, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map, const bool least_squares_grad, RatioMap& eigenvalues_ratio)
	{
		CalcHessianTaylor(graph, coord_map, func_map,dist_map, dist_thresh, tangent_basis_map);
		//CalcHessianCovDiff(graph, coord_map, m_dx, m_dy, dist_map, dist_thresh,tangent_basis_map);

		//find gradient threshold
		m_grad_norm.SetGraph(graph);
		m_hess_det.SetGraph(graph);
		std::vector<double> grad_values(num_vertices(graph));
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{
			m_grad_norm[*curr_vert] = sqrt(Sqr(dx[*curr_vert]) + Sqr(dy[*curr_vert]));
			grad_values[*curr_vert] = m_grad_norm[*curr_vert];
		}
		const size_t max_elem_num = Round(0.8 * grad_values.size());
		std::nth_element(grad_values.begin(), grad_values.begin() + max_elem_num, grad_values.end());
		const size_t min_elem_num = Round(0.3 * grad_values.size());
		std::nth_element(grad_values.begin(), grad_values.begin() + min_elem_num, grad_values.end());
		const double min_grad = *(std::min(grad_values.begin(), grad_values.end()));
		const double grad_thresh = std::min(grad_values[min_elem_num], (grad_values[max_elem_num] - min_grad)  * 0.01 + min_grad);
		std::cout << grad_thresh << " ";
		//determinant and eigenvalues
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{
			cv::Mat_<double> curr_hessian = m_hessian_map[*curr_vert];
			cv::Mat_<double> eigenvalues;
			cv::eigen(curr_hessian, eigenvalues);
			const double min_abs_eigen = std::min(abs(eigenvalues(0)), abs(eigenvalues(1)));
			const double max_abs_eigen = std::max(abs(eigenvalues(0)), abs(eigenvalues(1)));
			m_hess_det[*curr_vert] = curr_hessian(0, 0) * curr_hessian(1, 1) - curr_hessian(1, 0) * curr_hessian(0, 1);//curr_hessian(0, 0) + curr_hessian(1, 1);
			//eigenvalues_ratio[*curr_vert] = /*grad_norm*/determinant(curr_hessian) > 0.0 ? 
			//	Sqr(trace(curr_hessian))[0] * 10.0 - 11.0 * determinant(curr_hessian) : 1.0;*/
			if (min_abs_eigen == 0.0/* || eigenvalues(0) * eigenvalues(1) < 0.0*/)
			{
				eigenvalues_ratio[*curr_vert] = std::numeric_limits<double>::max();
			}
			else if (m_grad_norm[*curr_vert] < grad_thresh)
			{
				eigenvalues_ratio[*curr_vert] = 0;
			}
			else
			{
				eigenvalues_ratio[*curr_vert] = max_abs_eigen / min_abs_eigen;
			}
		}
	}
//private:
	DoublePropMap m_dx, m_dy, m_dxx, m_dxy, m_dyx, m_dyy;
	VecDoublePropMap m_grad_3d, m_grad_2d;
	HessianMap m_hessian_map;
	DoublePropMap m_grad_norm;
	DoublePropMap m_hess_det;
};

}