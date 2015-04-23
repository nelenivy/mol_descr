#pragma once
#include <vector>
#include <iostream>

#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"
#include "opencv2/core/core.hpp"

#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

//structure for keeping curvature values
struct Curvature
{
	Curvature() :  mean_curv(0), gaussian_curv(0) { }
	double mean_curv;
	double gaussian_curv;
};

struct CubicFittingCurvature
{
	static const int kEquationRowLen = 7;
	static const int kRowsForOneElement = 3;
};

struct ParabolloidFittingCurvature
{
	static const int kEquationRowLen = 3;
	static const int kRowsForOneElement = 1;
};

//class for imitating specialization of template method of template class
template<class Info3DType, class GraphType, class CurvatureType>
struct CurvatureCalculatorHelper;

template<class GraphType>
class CurvatureCalculator
{
public:
	typedef typename boost::graph_traits<GraphType>::vertex_descriptor VertexDescr;
	typedef typename boost::graph_traits<GraphType>::vertex_iterator VertexIter;
	template<class Info3DType, class GraphType, class CurvatureType> friend	struct CurvatureCalculatorHelper;
	//friend struct CurvatureCalculatorHelper<>;//<GraphType, ParabolloidCurvature>;
	//friend struct CurvatureCalculatorHelper<GraphType, CubicCurvature>;

	CurvatureCalculator(const int max_neighbours_num);
	template <typename CoordMapT>
	void CalculateCurvatureCubic(const GraphType& graph, const CoordMapT& coord_map, const VertexDescr mesh_vertice, Curvature& curvature_matr);
	template <typename CoordMapT>
	void CalculateCurvatureParabolloid(const GraphType& graph, const CoordMapT& coord_map, const VertexDescr mesh_vertice, Curvature& curvature_matr);
private:
	template <class CurvatureType, typename CoordMapT>
	void CalculateCurvatureImpl(const GraphType& graph, const CoordMapT& coord_map, const VertexDescr mesh_vertices, Curvature& curvature_matr);
	template <class Info3DType, class CurvatureType>
	void FillEquationMatrices(const Info3DType& curr_vertex, const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect);
	template <class CurvatureType>
	void ReservePlace(const int max_neighbours_num);

	cv::Mat_<double>	m_E, 
					m_rotat_mat,
					m_temp_v,
					m_r1,
					m_r2,
					m_r3,
					m_ans;
	cv::Mat_<double> m_curr_normal;
	cv::Mat_<double> m_curr_coord;
	cv::Mat_<double> m_neighb_normal;
	cv::Mat_<double> m_neighb_coord;
	//matrices for the equation solution
	cv::Mat_<double> m_coefficient_matrix,
				 m_z_vect; 
	int m_max_neighbours_num;
};
//////////////////////////////////////////////////////////////////////////
//class for imitating specialization of template method of template class
//template<class NodeType, class CurvatureType>
//struct CurvatureCalculatorHelper
//{
//	static void FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const NodeType& curr_vertex, 
//		const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect);
//};

template<class Info3DType, class NodeType>
struct CurvatureCalculatorHelper<Info3DType, NodeType, CubicFittingCurvature>
{
	static void FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const Info3DType& curr_vertex, 
		const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect);
};

template<class Info3DType, class NodeType>
struct CurvatureCalculatorHelper<Info3DType, NodeType, ParabolloidFittingCurvature>
{
	static void FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const Info3DType& curr_vertex, 
		const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect);
};
//////////////////////////////////////////////////////////////////////////
template <class NodeType>
CurvatureCalculator<NodeType>::CurvatureCalculator(const int max_neighbours_num)
{
	m_E			= cv::Mat_<double>::eye(3, 3); 
	m_rotat_mat	= cv::Mat_<double>::zeros(3, 3);
	m_temp_v	= cv::Mat_<double>::zeros(1, 3);
	m_r1		= cv::Mat_<double>::zeros(3, 1);
	m_r2		= cv::Mat_<double>::zeros(3, 1);
	m_r3		= cv::Mat_<double>::zeros(3, 1);
	m_max_neighbours_num = max_neighbours_num;
}

//"A Novel Cubic-Order Algorithm for Approximating Principal Direction Vectors"
template <class GraphType>
template <typename CoordMapT>
void CurvatureCalculator<GraphType>::CalculateCurvatureCubic(const GraphType& graph, const CoordMapT& coord_map,
	const VertexDescr mesh_vertices, Curvature& curvature_matr)
{
	CalculateCurvatureImpl<CubicFittingCurvature>(graph, coord_map, mesh_vertices, curvature_matr);
}

//"A Comparison of Gaussian and Mean Curvatures Estimation Methods on Triangular Meshes"
template <class GraphType>
template <typename CoordMapT>
void CurvatureCalculator<GraphType>::CalculateCurvatureParabolloid(const GraphType& graph, const CoordMapT& coord_map,
	const VertexDescr mesh_vertices, Curvature& curvature_matr)
{
	CalculateCurvatureImpl<ParabolloidFittingCurvature>(graph, coord_map, mesh_vertices, curvature_matr);
}

template <class GraphType>
template <class CurvatureType, typename CoordMapT>
void CurvatureCalculator<GraphType>::CalculateCurvatureImpl(const GraphType& graph, const CoordMapT& coord_map,
	const VertexDescr mesh_vertices, Curvature& curvature_values)
{
	typedef boost::graph_traits<GraphType>::adjacency_iterator AdjacencyIter;
	AdjacencyIter curr_neighbour, last;
	boost::tie(curr_neighbour, last) = adjacent_vertices(mesh_vertices, graph);
	if (curr_neighbour == last)
	{//if there are no neighbours we can't run algorithm
		curvature_values.gaussian_curv = 0.0;
		curvature_values.mean_curv = 0.0;
		return;
	}

	ReservePlace<CurvatureType>(m_max_neighbours_num);
	Point_ToMat_Transposed(coord_map[mesh_vertices].Normal(), m_curr_normal);
	Point_ToMat_Transposed(coord_map[mesh_vertices].Center(), m_curr_coord);
	//find coordinate transformation
	m_rotat_mat = m_E - m_curr_normal * m_curr_normal.t();

	m_r1(0) = 1;
	m_r1(1) = 0; 
	m_r1(2) = 0;			        
	m_r1 =  m_rotat_mat.t() * m_r1;
	m_r1 /= norm(m_r1);
	m_curr_normal.copyTo(m_r3);

	m_r2(0)=	m_r3(1) * m_r1(2) - m_r3(2) * m_r1(1);
	m_r2(1)= - (m_r3(0) * m_r1(2) - m_r3(2) * m_r1(0));
	m_r2(2)=	m_r3(0) * m_r1(1) - m_r3(1) * m_r1(0);

	m_r1.copyTo(m_rotat_mat.col(0));
	m_r2.copyTo(m_rotat_mat.col(1));
	m_r3.copyTo(m_rotat_mat.col(2));

	boost::tie(curr_neighbour, last) = adjacent_vertices(mesh_vertices, graph);
	const size_t neighb_num = out_degree(mesh_vertices, graph);
	ReservePlace<CurvatureType>(static_cast<int>(neighb_num));
	typedef decltype(coord_map[mesh_vertices]) Info3DType;
	for (; curr_neighbour != last; ++curr_neighbour)
	{ 
		const size_t neighb_ind = neighb_num - (last - curr_neighbour);
		FillEquationMatrices<Info3DType, CurvatureType>(
			coord_map[mesh_vertices], 
			coord_map[*curr_neighbour], 
			m_coefficient_matrix[CurvatureType::kRowsForOneElement * neighb_ind],
			m_z_vect[CurvatureType::kRowsForOneElement * neighb_ind]);
	}

	m_ans = m_coefficient_matrix.inv(DECOMP_SVD) * m_z_vect;
	curvature_values.gaussian_curv = m_ans(0) * m_ans(2) - m_ans(1) * m_ans(1);
	curvature_values.mean_curv = m_ans(0) + m_ans(2);
}

template <class Info3DType, class GraphType>
void CurvatureCalculatorHelper<Info3DType, GraphType, ParabolloidFittingCurvature>::FillEquationMatrices(CurvatureCalculator<GraphType>& calculator, 
	const Info3DType& curr_vertex, const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect)
{
	Point_ToMat_(curr_vertex.Center(), calculator.m_curr_coord);
	Point_ToMat_(curr_vertex.Normal(), calculator.m_curr_normal);
	Point_ToMat_(neighb_vertex.Center(), calculator.m_neighb_coord);
	Point_ToMat_(neighb_vertex.Normal(), calculator.m_neighb_normal);
	//transform neighbours coordinates and normal in new coordinates
	calculator.m_temp_v = (calculator.m_neighb_coord - calculator.m_curr_coord) * calculator.m_rotat_mat;                      
	const double x = calculator.m_temp_v(0);
	const double y = calculator.m_temp_v(1);
	const double z = calculator.m_temp_v(2);

	coeff_mat_row[0] = pow(x, 2) / 2.0;
	coeff_mat_row[1] = x * y;
	coeff_mat_row[2] = pow(y, 2) / 2.0;
	z_vect[0]	= z;
}

template <class Info3DType, class GraphType>
void CurvatureCalculatorHelper<Info3DType, GraphType, CubicFittingCurvature>::FillEquationMatrices(CurvatureCalculator<GraphType>& calculator,
	const Info3DType& curr_vertex, const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect)
{
	Point_ToMat_(curr_vertex.Center(), calculator.m_curr_coord);
	Point_ToMat_(curr_vertex.Normal(), calculator.m_curr_normal);
	Point_ToMat_(neighb_vertex.Center(), calculator.m_neighb_coord);
	Point_ToMat_(neighb_vertex.Normal(), calculator.m_neighb_normal);
	//transform neighbours coordinates and normal in new coordinates
	calculator.m_temp_v = (calculator.m_neighb_coord - calculator.m_curr_coord) * calculator.m_rotat_mat;                      
	const double x = calculator.m_temp_v(0);
	const double y = calculator.m_temp_v(1);
	const double z = calculator.m_temp_v(2);

	calculator.m_temp_v =  calculator.m_neighb_normal * calculator.m_rotat_mat;
	const double nx = calculator.m_temp_v(0);
	const double ny = calculator.m_temp_v(1);
	const double nz = calculator.m_temp_v(2);	

	coeff_mat_row[0] = pow(x, 2) / 2.0;
	coeff_mat_row[1] = x * y;
	coeff_mat_row[2] = pow(y, 2) / 2.0;
	coeff_mat_row[3] = pow(x, 3);
	coeff_mat_row[4] = pow(x, 2) * y;
	coeff_mat_row[5] = pow(y, 2) * x;
	coeff_mat_row[6] = pow(y, 3);

	coeff_mat_row += CubicFittingCurvature::kEquationRowLen;
	coeff_mat_row[0] = x;
	coeff_mat_row[1] = y;
	coeff_mat_row[2] = 0.0;
	coeff_mat_row[3] = 3.0 * pow(x, 2);
	coeff_mat_row[4] = 2.0 * x * y;
	coeff_mat_row[5] = pow(y, 2);
	coeff_mat_row[6] = 0.0;

	coeff_mat_row += CubicFittingCurvature::kEquationRowLen;
	coeff_mat_row[0] = 0;
	coeff_mat_row[1] = x;
	coeff_mat_row[2] = y;
	coeff_mat_row[3] = 0;
	coeff_mat_row[4] = pow(x, 2);
	coeff_mat_row[5] = 2 * x * y;
	coeff_mat_row[6] = 3 * pow(y, 2);

	z_vect[0]	= z;
	z_vect[1]	= - nx / nz;
	z_vect[2] = - ny / nz;
}

template <class GraphType>
template <class Info3DType, class CurvatureType>
void CurvatureCalculator<GraphType>::FillEquationMatrices(const Info3DType& curr_vertex, 
	const Info3DType& neighb_vertex, double* coeff_mat_row, double* z_vect)
{
	CurvatureCalculatorHelper<Info3DType,GraphType, CurvatureType>::FillEquationMatrices(*this, curr_vertex, neighb_vertex, coeff_mat_row, z_vect);
}

template <class NodeType>
template <class CurvatureType>
void CurvatureCalculator<NodeType>::ReservePlace(const int max_neighbours_num)
{
	m_coefficient_matrix = cv::Mat_<double>::zeros(CurvatureType::kRowsForOneElement * max_neighbours_num, CurvatureType::kEquationRowLen);
	m_z_vect = cv::Mat_<double>::zeros(CurvatureType::kRowsForOneElement * max_neighbours_num, 1); 
}

}