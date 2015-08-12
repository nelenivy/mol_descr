#pragma once

#include <utility>
#include "opencv2/core/core.hpp"
#include "CommonUtilities/common_functions.h"

namespace molecule_descriptor
{

class VectorOrientedCoordFinder
{
public:
	VectorOrientedCoordFinder()
	{
		m_E			= cv::Mat_<double>::eye(3, 3); 
		m_rotat_mat	= cv::Mat_<double>::zeros(3, 3);
		m_r1		= cv::Mat_<double>::zeros(3, 1);
		m_r2		= cv::Mat_<double>::zeros(3, 1);
		m_r3		= cv::Mat_<double>::zeros(3, 1);
	}
//normal is 3 dimensional vector
//new_basis is 3x3 matrix, where in columns are written basis vectors of new orthonormal coordinate system in old coordinates
//in new coordinate system normal is 3rd basis vectors
	void Process(const cv::Mat_<double>& normal, cv::Mat_<double>& new_basis)
	{
		//this matrix transforms vectors of unit basis into vectors, othogonal to vector normal
		m_rotat_mat = m_E - normal * normal.t();

		m_r1(0) = 1;
		m_r1(1) = 0; 
		m_r1(2) = 0;			        
		m_r1 =  m_rotat_mat.t() * m_r1;
		m_r1 /= norm(m_r1);
		normal.copyTo(m_r3);

		m_r2(0)=	m_r3(1) * m_r1(2) - m_r3(2) * m_r1(1);
		m_r2(1)= - (m_r3(0) * m_r1(2) - m_r3(2) * m_r1(0));
		m_r2(2)=	m_r3(0) * m_r1(1) - m_r3(1) * m_r1(0);

		new_basis.create(3, 3);
		m_r1.copyTo(new_basis.col(0));
		m_r2.copyTo(new_basis.col(1));
		m_r3.copyTo(new_basis.col(2));
	}
	void ProcessConnected(const cv::Mat_<double>& normal, cv::Mat_<double>& new_basis)
	{
		//IN ORDER TO FIND TANGENT COORDINATE SYSTEM WE ROTATE WORLD COORDINATE SYSTEM BY THE VECTOR IN OXY PLANE

		//1. WE FIND BASIS IN WHICH AXIS OF ROTATION IS ONE OF COORDINATES
		//3rd basis vector
		m_r3(0) = 0;
		m_r3(1) = 0; 
		m_r3(2) = 1;	
		//Axis of rotation - vector which is perpendicular to 3rd basis vector and normal
		m_r2(0)=	m_r3(1) * normal(2) - m_r3(2) * normal(1);
		m_r2(1)= - (m_r3(0) * normal(2) - m_r3(2) * normal(0));
		m_r2(2)=	m_r3(0) * normal(1) - m_r3(1) * normal(0);
		m_r2 /= norm(m_r2);
		//vector which is perpendicular to r2 and r3
		m_r1(0)=	m_r3(1) * m_r2(2) - m_r3(2) * m_r2(1);
		m_r1(1)= - (m_r3(0) * m_r2(2) - m_r3(2) * m_r2(0));
		m_r1(2)=	m_r3(0) * m_r2(1) - m_r3(1) * m_r2(0);
		m_r1 /= norm(m_r1);
		//in this basis matrix has the next view^
		//    ( cos(a) 0   -sin(a))
		//R = ( 0      1    0     )
		//    ( sin(a) 0    cos(a))
		//R * r3 = normal => (-sin(a), 0, cos(a)) = normal => -sin(a) = (normal. r1), cos(a) = (normal. r3)
		const double sin_alpha = - normal.dot(m_r1) / norm(normal);
		const double cos_alpha = normal.dot(m_r3) / norm(normal);
// 		std::cout << sin_alpha << " " << cos_alpha << " " << Sqr(cos_alpha) + Sqr(sin_alpha) << std::endl;
		//this matrix rotates basis around vector r2 in basis r1, r2, r3
		m_rotat_mat.setTo(0);
		m_rotat_mat(0, 0) = m_rotat_mat(2, 2) = cos_alpha;
		m_rotat_mat(0, 2) = -sin_alpha;
		m_rotat_mat(2, 0) = sin_alpha;
		m_rotat_mat(1, 1) = 1.0;
		//matrix of transforming coordinates to r1, r2, r3
		m_r1.copyTo(m_E.col(0));
		m_r2.copyTo(m_E.col(1));
		m_r3.copyTo(m_E.col(2));
		//matrix of rotation in initial coordinates
		m_rotat_mat = m_E.inv(DECOMP_SVD) * m_rotat_mat * m_E;
		//find images of basis vectors by this rotation
		m_r1(0) = 1;
		m_r1(1) = 0; 
		m_r1(2) = 0;	
		m_r1 = m_rotat_mat * m_r1;

		m_r2(0) = 0;
		m_r2(1) = 1; 
		m_r2(2) = 0;	
		m_r2 = m_rotat_mat * m_r2;

		m_r3(0) = 0;
		m_r3(1) = 0; 
		m_r3(2) = 1;	
		m_r3 = m_rotat_mat * m_r3;
		//std::cout << m_r3 << normal << "\n";
		//final matrix of coordinates transform
		new_basis.create(3, 3);
		m_r1.copyTo(new_basis.col(0));
		m_r2.copyTo(new_basis.col(1));
		m_r3.copyTo(new_basis.col(2));
		new_basis = new_basis.inv(DECOMP_SVD);
		//std::cout << new_basis * normal << " " << sin_alpha << " " << cos_alpha << "\n";
	}
private:
	cv::Mat_<double> m_E, 
		m_rotat_mat,
		m_r1,
		m_r2,
		m_r3;
};

template <class Graph, class CoordMap, class NormalMap, class BasisMap>
void CalcTangentCoordSystemMap(const Graph& graph, const CoordMap& coord_map, const NormalMap& norm_map, BasisMap& basis_map)
{
	VectorOrientedCoordFinder coord_finder;
	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		if (adjacent_vertices(*curr_vert, graph).first == adjacent_vertices(*curr_vert, graph).second)
		{//if there are no neighbours we can't run algorithm
			basis_map[*curr_vert] = cv::Mat_<double>::eye(3, 3); 
			return;
		}

		cv::Mat_<double> curr_normal, curr_coord;
		Point_ToMat_Transposed(norm_map[*curr_vert], curr_normal);
		Point_ToMat_Transposed(coord_map[*curr_vert], curr_coord);
		//find coordinate transformation
		cv::Mat_<double> basis;
		coord_finder./*Process*/ProcessConnected(curr_normal, basis);
		basis_map[*curr_vert] = basis;
	}
}

}