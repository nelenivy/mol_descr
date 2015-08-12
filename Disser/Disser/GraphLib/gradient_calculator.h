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
void CalcGradientInVertLeastSquares(const Graph& graph, const CoordMap& coord_map, 
						const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
						const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,
						const cv::Mat_<double>& tangent_basis, cv::Mat_<double>& grad)//vectors in columns, 3rd vector is normal to surface
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	std::vector<VertexDescriptor> vertices_within_dist;
	GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
	cv::Mat_<double> cent_coord;
	Point_ToMat_Transposed(coord_map[cent_vert], cent_coord);
	//FILL EQUATION FOR GRADIENT CALCULATION
	cv::Mat_<double> projected_vect, neighb_coord;

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
		projected_vect = tangent_basis*(neighb_coord - cent_coord);
		coeffs(curr_row, 0) = projected_vect(0,0);
		coeffs(curr_row, 1) = projected_vect(0,1);

		const double proj_length = sqrt(Sqr(projected_vect(0, 0)) + Sqr(projected_vect(0, 1)));
		derivatives(curr_row, 0) = (func_map[*neighb_vert] - func_map[cent_vert]) / proj_length;
	}

	//SOLVE EQUATION
	grad = coeffs.inv(DECOMP_SVD) * derivatives;
}

enum Direction {dX = 0, dY = 1};

template <class Graph, class CoordMap, class FunctionMap, class DistanceMap>
double CalcDirectDerivInVert(const Graph& graph, const CoordMap& coord_map, 
							const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
							const typename boost::graph_traits<Graph>::vertex_descriptor cent_vert,
							const cv::Mat_<double>& tangent_basis, const Direction direction)//vectors in columns, 3rd vector is normal to surface
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	std::vector<VertexDescriptor> vertices_within_dist;
	GetVerticesWithinDistPlusAdjacent(cent_vert, graph, dist_map, dist_thresh, vertices_within_dist);
	cv::Mat_<double> cent_coord, neighb_coord, projected_vect;
	Point_ToMat_Transposed(coord_map[cent_vert], cent_coord);
	//FILL EQUATION FOR GRADIENT CALCULATION
	double derivative = 0;
	size_t used_vertices = 0;
	for (auto neighb_vert = vertices_within_dist.begin(), end_vert = vertices_within_dist.end(); neighb_vert != end_vert; ++neighb_vert)
	{
		Point_ToMat_Transposed(coord_map[*neighb_vert], neighb_coord);
		projected_vect = tangent_basis*(neighb_coord - cent_coord);

		const double curr_coord = projected_vect(direction);
		if (curr_coord > 0 && *neighb_vert != cent_vert)
		{
			++used_vertices;
			derivative += (func_map[*neighb_vert] - func_map[cent_vert]) / curr_coord;
		}
	}

	derivative /= used_vertices;
	return derivative;
}

template <class Graph, class CoordMap, class FunctionMap, class DistanceMap,class TangentBasisMap, class GradientMap>
void CalcGradientMaps(const Graph& graph, const CoordMap& coord_map, 
					 const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh, 
					 const TangentBasisMap& tangent_basis_map, const bool least_squares_grad, 
					 GradientMap& grad_x, GradientMap& grad_y)
{
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		if (least_squares_grad)
		{
			cv::Mat_<double> curr_grad;
			CalcGradientInVertLeastSquares(graph, coord_map, func_map, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], curr_grad);
			grad_x[*curr_vert] = curr_grad(0);
			grad_y[*curr_vert] = curr_grad(1);
		}
		else
		{
			grad_x[*curr_vert] = CalcDirectDerivInVert(
				graph, coord_map, func_map, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], dX);
			grad_y[*curr_vert] = CalcDirectDerivInVert(
				graph, coord_map, func_map, dist_map, dist_thresh, *curr_vert, tangent_basis_map[*curr_vert], dY);
		}
	}
}

class HessianMatrixCalculator
{
public:
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoublePropMap;

	template <class Graph, class CoordMap, class FunctionMap, class DistanceMap,class TangentBasisMap, class RatioMap>
	void Process(const Graph& graph, const CoordMap& coord_map, 
		const FunctionMap& func_map, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map, const bool least_squares_grad,RatioMap& eigenvalues_ratio)
	{
		m_dx.SetGraph(graph);
		m_dy.SetGraph(graph);		
		CalcGradientMaps(graph, coord_map, func_map,dist_map, dist_thresh, tangent_basis_map, least_squares_grad, m_dx, m_dy);
		ProcessFromGradient(graph, coord_map, m_dx, m_dy, dist_map, dist_thresh,tangent_basis_map, least_squares_grad, eigenvalues_ratio);
	}
	template <class Graph, class CoordMap, class GradMap, class DistanceMap,class TangentBasisMap, class RatioMap>
	void ProcessFromGradient(const Graph& graph, const CoordMap& coord_map, 
		const GradMap& dx, const GradMap& dy, const DistanceMap& dist_map, const double dist_thresh,
		const TangentBasisMap& tangent_basis_map, const bool least_squares_grad,RatioMap& eigenvalues_ratio)
	{
		m_dxx.SetGraph(graph);
		m_dxy.SetGraph(graph);
		m_dyx.SetGraph(graph);
		m_dyy.SetGraph(graph);
		CalcGradientMaps(graph, coord_map, dx,dist_map, dist_thresh, tangent_basis_map, least_squares_grad, m_dxx, m_dxy);
		CalcGradientMaps(graph, coord_map, dy,dist_map, dist_thresh, tangent_basis_map, least_squares_grad, m_dyx, m_dyy);

		//find gradient threshold
		m_grad_norm.SetGraph(graph);
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
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
		{
			cv::Mat_<double> curr_hessian(2, 2);
			curr_hessian(0,0) = m_dxx[*curr_vert];
			curr_hessian(1,0) = m_dyx[*curr_vert];
			curr_hessian(0,1) = m_dyx[*curr_vert];//m_dxy[*curr_vert];
			curr_hessian(1,1) = m_dyy[*curr_vert];

			cv::Mat_<double> eigenvalues;
			cv::eigen(curr_hessian, eigenvalues);
			const double min_abs_eigen = std::min(abs(eigenvalues(0)), abs(eigenvalues(1)));
			const double max_abs_eigen = std::max(abs(eigenvalues(0)), abs(eigenvalues(1)));

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
	DoublePropMap m_grad_norm;
};

}