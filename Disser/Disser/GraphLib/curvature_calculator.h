#pragma once
#include <vector>
#include <iostream>

#include "opencv2/core/core.hpp"

#include "graph_structures.h"
#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

struct CubicCurvature
{
	static const int kEquationRowLen = 7;
	static const int kRowsForOneElement = 3;
};

struct ParabolloidCurvature
{
	static const int kEquationRowLen = 3;
	static const int kRowsForOneElement = 1;
};

//class for imitating specialization of template method of template class
template<class NodeType, class CurvatureType>
struct CurvatureCalculatorHelper;

template<class NodeType>
class CurvatureCalculator
{
public:
	friend struct CurvatureCalculatorHelper<NodeType, ParabolloidCurvature>;
	friend struct CurvatureCalculatorHelper<NodeType, CubicCurvature>;

	CurvatureCalculator(const int max_neighbours_num);
	void CalculateCurvatureCubic(const GraphNode<NodeType>& mesh_vertices, Vec2d& curvature_matr);
	void CalculateCurvatureParabolloid(const GraphNode<NodeType>& mesh_vertices, Vec2d& curvature_matr);
private:
	template <class CurvatureType>
	void CalculateCurvatureImpl(const GraphNode<NodeType>& mesh_vertices, Vec2d& curvature_matr);
	template <class CurvatureType>
	void FillEquationMatrices(const NodeType& curr_vertex, const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect);
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

template<class NodeType>
struct CurvatureCalculatorHelper<NodeType, CubicCurvature>
{
	static void FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const NodeType& curr_vertex, 
		const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect);
};

template<class NodeType>
struct CurvatureCalculatorHelper<NodeType, ParabolloidCurvature>
{
	static void FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const NodeType& curr_vertex, 
		const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect);
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

template <class NodeType>
void CurvatureCalculator<NodeType>::CalculateCurvatureCubic(const GraphNode<NodeType>& mesh_vertices, cv::Vec2d& curvature_matr)
{
	CalculateCurvatureImpl<CubicCurvature>(mesh_vertices, curvature_matr);
}

template <class NodeType>
void CurvatureCalculator<NodeType>::CalculateCurvatureParabolloid(const GraphNode<NodeType>& mesh_vertices, cv::Vec2d& curvature_matr)
{
	CalculateCurvatureImpl<ParabolloidCurvature>(mesh_vertices, curvature_matr);
}

template <class NodeType>
template <class CurvatureType>
void CurvatureCalculator<NodeType>::CalculateCurvatureImpl(const GraphNode<NodeType>& mesh_vertex, cv::Vec2d& curvature_values)
{
	if (mesh_vertex.neighbours.size() == 0)
	{//if there are no neighbours we can't run algorithm
		curvature_values[0] = 0.0;
		curvature_values[1] = 0.0;
		return;
	}

	ReservePlace<CurvatureType>(m_max_neighbours_num);
	Point_ToMat_Transposed(mesh_vertex.element->Normal(), m_curr_normal);
	Point_ToMat_Transposed(mesh_vertex.element->Center(), m_curr_coord);
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

	//write normal and coordinates to the new coordinate system
	const auto& curr_neighbours = mesh_vertex.neighbours;
	ReservePlace<CurvatureType>(static_cast<int>(curr_neighbours.size()));

	for (int neighb_ind = 0; neighb_ind < curr_neighbours.size(); neighb_ind++)
	{ 
		FillEquationMatrices<CurvatureType>(
			*(mesh_vertex.element), 
			*(mesh_vertex.neighbours[neighb_ind]->element), 
			m_coefficient_matrix[CurvatureType::kRowsForOneElement * neighb_ind],
			m_z_vect[CurvatureType::kRowsForOneElement * neighb_ind]);
	}

	m_ans = m_coefficient_matrix.inv(DECOMP_SVD) * m_z_vect;
	curvature_values[0] = m_ans(0) * m_ans(2) - m_ans(1) * m_ans(1);
	curvature_values[1] = m_ans(0) + m_ans(2);
}

template <class NodeType>
void CurvatureCalculatorHelper<NodeType, ParabolloidCurvature>::FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const NodeType& curr_vertex, 
	const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect)
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

template <class NodeType>
void CurvatureCalculatorHelper<NodeType, CubicCurvature>::FillEquationMatrices(CurvatureCalculator<NodeType>& calculator, const NodeType& curr_vertex, 
	const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect)
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

	coeff_mat_row += CubicCurvature::kEquationRowLen;
	coeff_mat_row[0] = x;
	coeff_mat_row[1] = y;
	coeff_mat_row[2] = 0.0;
	coeff_mat_row[3] = 3.0 * pow(x, 2);
	coeff_mat_row[4] = 2.0 * x * y;
	coeff_mat_row[5] = pow(y, 2);
	coeff_mat_row[6] = 0.0;

	coeff_mat_row += CubicCurvature::kEquationRowLen;
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

template <class NodeType>
template <class CurvatureType>
void CurvatureCalculator<NodeType>::FillEquationMatrices(const NodeType& curr_vertex, 
	const NodeType& neighb_vertex, double* coeff_mat_row, double* z_vect)
{
	CurvatureCalculatorHelper<NodeType, CurvatureType>::FillEquationMatrices(*this, curr_vertex, neighb_vertex, coeff_mat_row, z_vect);
}

template <class NodeType>
template <class CurvatureType>
void CurvatureCalculator<NodeType>::ReservePlace(const int max_neighbours_num)
{
	m_coefficient_matrix = cv::Mat_<double>::zeros(CurvatureType::kRowsForOneElement * max_neighbours_num, CurvatureType::kEquationRowLen);
	m_z_vect = cv::Mat_<double>::zeros(CurvatureType::kRowsForOneElement * max_neighbours_num, 1); 
}

}