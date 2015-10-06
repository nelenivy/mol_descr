#pragma once

#include <utility>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/tags.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/bindings/ublas.hpp>
//#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "opencv2/core/core.hpp"
#include "CommonUtilities/common_functions.h"
#include "GraphLib/interpolator.h"

namespace molecule_descriptor
{

//vector, which is perpendicular to normal and r3
inline void VectorProduct(const cv::Mat_<double>& vec_1, const cv::Mat_<double>& vec_2, cv::Mat_<double>& vec_prod)
{
	if (vec_prod.total() != 3)
	{
		vec_prod.create(3, 1);
	}
	vec_prod(0)=	vec_1(1) * vec_2(2) - vec_1(2) * vec_2(1);
	vec_prod(1)= - (vec_1(0) * vec_2(2) - vec_1(2) * vec_2(0));
	vec_prod(2)=	vec_1(0) * vec_2(1) - vec_1(1) * vec_2(0);
}

//3rd basis vector in both basises 
inline void BasisesSum(const cv::Mat_<double>& basis_1, const cv::Mat_<double>& basis_2, cv::Mat_<double>& basises_sum)
{
	basises_sum.create(3, 3);
	basises_sum = (basis_1 + basis_2) * 0.5;

	const double kMinNorm = 0.00001;
	bool collinear = false;
	for (int col_num = 0; col_num < 2; ++col_num)
	{
		if (cv::norm(basises_sum.col(col_num)) < kMinNorm)
		{
			collinear = true;
		}
	}
	
	if (collinear )
	{
		basis_1.col(0).copyTo(basises_sum.col(1));
		basis_1.col(1).copyTo(basises_sum.col(0));
		basis_1.col(2).copyTo(basises_sum.col(2));
	}
	else
	{
		for (int col_num = 0; col_num < 3; ++col_num)
		{
			basises_sum.col(col_num) /= cv::norm(basises_sum.col(col_num));
		}
	}
}


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
	
	//function computes the matrix which performs rotation around the vector which is perpendicular to source and dest vectors
	void ShortestRotationMatrix(const cv::Mat_<double>& source, const cv::Mat_<double>& dest, const bool normalized, cv::Mat_<double>& rotation_mat)
	{
		//1. WE FIND BASIS IN WHICH AXIS OF ROTATION IS ONE OF COORDINATES
		//3rd basis vector
		if (!normalized)
		{
			m_r3 = source * (1.0 / cv::norm(source));
			m_dest = dest * (1.0 / cv::norm(dest));
		}
		const cv::Mat_<double>& r3_old = normalized ? source : m_r3;
		const cv::Mat_<double>& r3_new =  normalized ? dest : m_dest;
		//Axis of rotation - vector which is perpendicular to 3rd basis vector and normal
		VectorProduct(r3_old, r3_new, m_r2);
		m_r2 /= norm(m_r2);
		//vector which is perpendicular to r2 and r3
		VectorProduct(r3_old, m_r2, m_r1);		
		m_r1 /= norm(m_r1);
		//in this basis matrix has the next view^
		//    ( cos(a) 0   -sin(a))
		//R = ( 0      1    0     )
		//    ( sin(a) 0    cos(a))
		//R * r3 = normal => (-sin(a), 0, cos(a)) = normal => -sin(a) = (normal. r1), cos(a) = (normal. r3)
		const double sin_alpha = - r3_new.dot(m_r1);
		const double cos_alpha = r3_new.dot(r3_old);
		//std::cout << sin_alpha << " " << cos_alpha << " " << Sqr(cos_alpha) + Sqr(sin_alpha) << std::endl;
		//this matrix rotates basis around vector r2 in basis r1, r2, r3
		rotation_mat.create(3, 3);
		rotation_mat.setTo(0);
		rotation_mat(0, 0) = rotation_mat(2, 2) = cos_alpha;
		rotation_mat(0, 2) = -sin_alpha;
		rotation_mat(2, 0) = sin_alpha;
		rotation_mat(1, 1) = 1.0;

		//matrix of transforming coordinates to r1, r2, r3
		m_r1.copyTo(m_E.col(0));
		m_r2.copyTo(m_E.col(1));
		r3_old.copyTo(m_E.col(2));
		//matrix of rotation in initial coordinates
		//std::cout << m_E.inv(DECOMP_SVD) * m_r1 << m_E.inv(DECOMP_SVD) * m_r2 << m_E.inv(DECOMP_SVD) * m_r3 << "\n";
		//m_rotat_mat = m_E.inv(DECOMP_SVD) * m_rotat_mat * m_E;
		rotation_mat = m_E * rotation_mat * m_E.inv(DECOMP_SVD);
	}
	void ProcessShortestRotation(const cv::Mat_<double>& old_basis, const cv::Mat_<double>& normal, cv::Mat_<double>& new_basis)
	{
		const bool kNormalized = false;
		ShortestRotationMatrix(old_basis.col(2), normal, kNormalized, m_rotat_mat);		
		//find images of basis vectors by this rotation
		new_basis = m_rotat_mat * old_basis;
		//m_r1(0) = 1;
		//m_r1(1) = 0; 
		//m_r1(2) = 0;	
		//m_r1 = m_rotat_mat * m_r1;

		//m_r2(0) = 0;
		//m_r2(1) = 1; 
		//m_r2(2) = 0;	
		//m_r2 = m_rotat_mat * m_r2;

		//m_r3(0) = 0;
		//m_r3(1) = 0; 
		//m_r3(2) = 1;	
		//m_r3 = m_rotat_mat * m_r3;
		//std::cout << cv::norm(m_r1) << " " << cv::norm(m_r2) << " " << cv::norm(m_r3) << "\n";
		////final matrix of coordinates transform
		//new_basis.create(3, 3);
		//m_r1.copyTo(new_basis.col(0));
		//m_r2.copyTo(new_basis.col(1));
		//m_r3.copyTo(new_basis.col(2));
	}
private:
	cv::Mat_<double> m_E, 
		m_rotat_mat,
		m_r1,
		m_r2,
		m_r3,
		m_dest;
};

template <class Graph, class NormalMap, class CoordMap, class BasisMap>
void CalcConsistentBasises(const Graph& graph, const CoordMap& coord_map, const NormalMap& norm_map, const typename boost::graph_traits<Graph>::vertex_descriptor seed, 
						   BasisMap& basis_map)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	ContPropMap<Graph, std::vector<uint8_t>, VERTEX> visited_vert(graph);
	visited_vert.SetTo(0);
	ContPropMap<Graph, std::vector<uint8_t>, VERTEX> finished_vert(graph);
	finished_vert.SetTo(0);
	ContPropMap<Graph, std::vector<double>, VERTEX> visits_count(graph);
	visits_count.SetTo(0);
	ContPropMap<Graph, std::vector<cv::Mat_<double>>, VERTEX> basises_sum(graph);
	VectorOrientedCoordFinder vert_coord_finder;
	std::deque<vertex_descriptor> vertices_queue;
	vertices_queue.push_back(seed);
	Mat_<double> curr_basis;
	while (!vertices_queue.empty())
	{
		vertex_descriptor curr_vertice = vertices_queue.front();
		vertices_queue.pop_front();
		visited_vert[curr_vertice] = 1;
		finished_vert[curr_vertice] = 1;

		if (visits_count[curr_vertice] > 0.0)
		{
			//std::cout << visits_count[curr_vertice] << " ";

			basises_sum[curr_vertice].copyTo(basis_map[curr_vertice]);
			for (int col_num = 0; col_num < 3; ++col_num)
			{
				basis_map[curr_vertice].col(col_num) *= (1.0 / /*visits_count[curr_vertice]);*/cv::norm(basis_map[curr_vertice].col(col_num)));
			}
		}

		const Mat_<double>& curr_basis = basis_map[curr_vertice];
		/*basis_map[curr_vertice].copyTo(curr_basis);
		for (int col_num = 0; col_num < 3; ++col_num)
		{
		curr_basis.col(col_num) *= 1.0 / cv::norm(curr_basis.col(col_num));
		}*/
		
		cv::Mat_<double> curr_normal;
		Point_ToMat_Transposed(norm_map[curr_vertice], curr_normal);

		for (auto curr_neighb = adjacent_vertices(curr_vertice, graph).first,
			end_neighb = adjacent_vertices(curr_vertice, graph).second; curr_neighb != end_neighb; ++curr_neighb)
		{
			if (finished_vert[*curr_neighb] == 0)
			{
				cv::Mat_<double> neighb_basis, neighb_normal;
				Point_ToMat_Transposed(norm_map[*curr_neighb], neighb_normal);
				neighb_normal /= cv::norm(neighb_normal);
				vert_coord_finder.ProcessShortestRotation(curr_basis, neighb_normal, neighb_basis);
				//FIND WEIGHTS
				cv::Point3d connect = coord_map[*curr_neighb] - coord_map[curr_vertice];
				cv::Mat_<double> vect_connect;
				Point_ToMat_Transposed(connect, vect_connect);
				vect_connect /= cv::norm(vect_connect);
				cv::Mat_<double> vect_prod;
				VectorProduct(curr_normal, neighb_normal, vect_prod);
				vect_prod /= cv::norm(vect_prod);

				const double prod = std::min(1.0, std::max(0.0, 1.0 - vect_prod.dot(vect_connect)));
				const double weight = pow(prod, 10.0);

				neighb_basis *= weight;

				if (visited_vert[*curr_neighb] == 0)
				{
					basises_sum[*curr_neighb] = neighb_basis;
					vertices_queue.push_back(*curr_neighb);
				}
				else
				{
					basises_sum[*curr_neighb] += neighb_basis;
				}
				visited_vert[*curr_neighb] = 1;
				visits_count[*curr_neighb] += 1;//weight;
			}
		}
	}	
}

template <class Graph, class BasisMap>
void MakeBasisesConsistent(const Graph& graph, const typename boost::graph_traits<Graph>::vertex_descriptor seed, BasisMap& basis_map)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	ContPropMap<Graph, std::vector<uint8_t>, VERTEX> visited_vert(graph);
	ContPropMap<Graph, std::vector<uint8_t>, VERTEX> visited_first_time_vert(graph);
	ContPropMap<Graph, std::vector<uint8_t>, VERTEX> errors_vert(graph);
	ContPropMap<Graph, std::vector<std::vector<cv::Mat_<double>>>, VERTEX> uncoin_basises(graph);
	std::fill(visited_vert.begin(), visited_vert.end(), 0);
	std::fill(visited_first_time_vert.begin(), visited_first_time_vert.end(), 0);
	std::deque<vertex_descriptor> vertices_queue;
	vertices_queue.push_back(seed);

	while (!vertices_queue.empty())
	{
		vertex_descriptor curr_vertice = vertices_queue.front();
		vertices_queue.pop_front();
		const Mat_<double> curr_basis = (basis_map[curr_vertice]).inv(DECOMP_SVD);
		visited_vert[curr_vertice] = 1;

		for (auto curr_neighb = adjacent_vertices(curr_vertice, graph).first,
			end_neighb = adjacent_vertices(curr_vertice, graph).second; curr_neighb != end_neighb; ++curr_neighb)
		{
			if (visited_vert[*curr_neighb] == 0)
			{
				Mat_<double> neighb_basis = (basis_map[*curr_neighb]).inv(DECOMP_SVD);
				const int8_t normals_dot_sign = Sign(curr_basis.col(2).dot(neighb_basis.col(2)));
				bool errors = false;
				for (int col_num = 0; col_num < 2; ++col_num)
				{
					const int8_t curr_dot_sign = Sign(curr_basis.col(col_num).dot(neighb_basis.col(col_num)));

					if (curr_dot_sign != normals_dot_sign)
					{
						neighb_basis.col(col_num) *= -1;
						errors= true;
						/*if (visited_first_time_vert[*curr_neighb] == 1)
						{
							std::cout << curr_basis.col(col_num) << " " << neighb_basis.col(col_num) << "\n";
						}*/
					}
				}
				basis_map[*curr_neighb] = neighb_basis.inv(DECOMP_SVD);

				if (errors)
				{
					uncoin_basises[*curr_neighb].push_back(curr_basis);
				}
				if (visited_first_time_vert[*curr_neighb] == 1 && errors)
				{
					errors_vert[*curr_neighb] = 1;
				}
				else if (visited_first_time_vert[*curr_neighb] == 0)
				{
					visited_first_time_vert[*curr_neighb] = 1;
					vertices_queue.push_back(*curr_neighb);
				}
			}
		}
	}	

	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		if (uncoin_basises[*curr_vert].size() > 1)
		{
			std::cout << basis_map[*curr_vert].inv(DECOMP_SVD);
			
			for (size_t ind = 0; ind < uncoin_basises[*curr_vert].size(); ++ind)
			{
				std::cout << "\n" << uncoin_basises[*curr_vert][ind];
			}
			std::cout << "\n \n";

		}
	}
	double errors_num = 0;
	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		errors_num += (errors_vert[*curr_vert] == 1);
	}

	std::cout << errors_num / num_vertices(graph) << "\n";
}

template <class Graph, class CoordMap, class NormalMap, class BasisMap>
void CalcTangentCoordSystemMap(const Graph& graph, const CoordMap& coord_map, const NormalMap& norm_map, 
							   const double norm_thresh,BasisMap& basis_map)
{
	const int kItersNum = 600;
	std::vector<ContPropMap<Graph, std::vector<cv::Mat_<double>>, VERTEX>> basises(
		kItersNum, ContPropMap<Graph, std::vector<cv::Mat_<double>>, VERTEX>(graph));
	const int step = (num_vertices(graph) - 1) / kItersNum;

	for (int iter = 0; iter < kItersNum; ++iter)
	{
		auto first_vert = boost::vertices(graph).first + iter * step;
		auto end_vert = boost::vertices(graph).second;
		for (; first_vert != end_vert; ++first_vert)
		{
			if (adjacent_vertices(*first_vert, graph).first == adjacent_vertices(*first_vert, graph).second)
			{//if there are no neighbours we can't run algorithm
				basis_map[*first_vert] = cv::Mat_<double>::eye(3, 3); 
			}
			else
			{
				break;
			}
		}
		if (first_vert == end_vert)
		{
			return;
		}

		if (iter == 0)
		{
			cv::Mat_<double> curr_normal;
			Point_ToMat_Transposed(norm_map[*first_vert], curr_normal);
			curr_normal /= cv::norm(curr_normal);
			cv::Mat_<double> basis;
			VectorOrientedCoordFinder coord_finder;
			coord_finder.Process(curr_normal, basis);
			basises[iter][*first_vert] = basis;
		}
		else
		{
			basises[iter - 1][*first_vert].copyTo(basises[iter][*first_vert]);
			for (int col_num = 0; col_num < 3; ++col_num)
			{
				basises[iter][*first_vert].col(col_num) *= (1.0 / cv::norm(basises[iter][*first_vert].col(col_num)));
			}
		}

		CalcConsistentBasises(graph, coord_map, norm_map, *first_vert, basises[iter]);
	}

	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		basis_map[*curr_vert] = basises[0][*curr_vert];
		
		for (int iter = 1; iter < kItersNum; ++iter)
		{
			basis_map[*curr_vert] += basises[iter][*curr_vert];
		}
		basis_map[*curr_vert] *= (1.0 / kItersNum);

		for (int col_num = 0; col_num < 3; ++col_num)
		{
			basis_map[*curr_vert].col(col_num) *= (1.0 / cv::norm(basis_map[*curr_vert].col(col_num)));
		}

		/*for (int col_num = 0; col_num < 3; ++col_num)
		{
		std::cout << (basis_map[*curr_vert].col(col_num)).dot(basis_map[*curr_vert].col((col_num + 1) % 3)) << " ";
		}
		cv::Mat_<double> curr_normal;
		Point_ToMat_Transposed(norm_map[*curr_vert], curr_normal);
		std::cout << cv::norm(basis_map[*curr_vert].col(2) - curr_normal) << "\n";*/
		basis_map[*curr_vert] = basis_map[*curr_vert].inv(DECOMP_SVD);
	}
}



template <class Graph, class CoordMap, class NormalMap, class BasisMap>
void CalcTangentCoordSystemMap_SystemSolving(const Graph& graph, const CoordMap& coord_map, const NormalMap& norm_map, BasisMap& basis_map)
{
	namespace ublas = boost::numeric::ublas;
	namespace umf = boost::numeric::bindings::umfpack;
	//determine a number of non-zero elements and equations number
	int equations_num = 1;
	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		equations_num += out_degree(*curr_vert, graph);
	}
	equations_num *= 3;
	const int non_zero_elements = equations_num * 4;
	const int variables_num = 3 * num_vertices(graph);
	//
	//ublas::compressed_matrix<double, ublas::column_major, 0,
	//	ublas::unbounded_array<int>, ublas::unbounded_array<double> > coeff_matrix(equations_num, variables_num, non_zero_elements); 
	//ublas::vector<double> zero_column(equations_num), variables(variables_num);
	cv::Mat_<double> coeff_matrix(equations_num, variables_num), zero_column(equations_num, 1), variables(variables_num, 1); 	
	std::fill(zero_column.begin(), zero_column.end(), 0.0);
	coeff_matrix.setTo(0);
	//FILL MATRICES
	VectorOrientedCoordFinder coord_finder;
	cv::Mat_<double> curr_normal, neighb_normal, rotation_mat;
	const auto index_map = boost::get(boost::vertex_index, graph);
	size_t curr_equation = 0;
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		Point_ToMat_Transposed(norm_map[*curr_vert], curr_normal);
		curr_normal /= cv::norm(curr_normal);
		const size_t curr_index = index_map[*curr_vert];

		if (curr_index == 0)
		{
			coord_finder.Process(curr_normal, rotation_mat);

			for (size_t curr_coord = 0; curr_coord < 3; ++curr_coord)
			{
				coeff_matrix(curr_coord, curr_coord) = 1;
				zero_column(curr_coord) = rotation_mat(curr_coord, 0);
			}
			++curr_equation;
		}
		for (auto curr_neighb = adjacent_vertices(*curr_vert, graph).first, end_neighb = adjacent_vertices(*curr_vert, graph).second;
			curr_neighb != end_neighb; ++curr_neighb, ++curr_equation)
		{
			Point_ToMat_Transposed(norm_map[*curr_neighb], neighb_normal);
			neighb_normal /= cv::norm(neighb_normal);
			coord_finder.ShortestRotationMatrix(curr_normal, neighb_normal, true, rotation_mat);
			const size_t neighb_index = index_map[*curr_neighb];

			for (size_t curr_row = 0; curr_row < 3; ++curr_row)
			{
				for (size_t curr_coord = 0; curr_coord < 3; ++curr_coord)
				{
					coeff_matrix(3 * curr_equation + curr_row, 3* curr_index + curr_coord) = rotation_mat(curr_row, curr_coord);
				}

				coeff_matrix(3 * curr_equation + curr_row, 3* neighb_index + curr_row) = -1.0;
			}
		}
	}

	/*umf::symbolic_type<double> symbolic_p;
	umf::numeric_type<double> numeric_p;

	umf::symbolic (coeff_matrix, symbolic_p); 
	umf::numeric (coeff_matrix, symbolic_p, numeric_p); 
	umf::solve (coeff_matrix, variables, zero_column, numeric_p);   */
	//FILL BASIS FROM SOLUTION
	variables = coeff_matrix.inv(DECOMP_QR) * zero_column;
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		cv::Mat_<double> curr_basis(3, 3);
		const size_t curr_index = index_map[*curr_vert];
		Point_ToMat_Transposed(norm_map[*curr_vert], curr_normal);
		curr_normal /= cv::norm(curr_normal);
		curr_normal.copyTo(curr_basis.col(2));

		for (size_t curr_coord = 0; curr_coord < 3; ++curr_coord)
		{
			curr_basis(curr_coord, 1) = variables(3 * curr_index + curr_coord);
		}
		curr_basis.col(1) /= cv::norm(curr_basis.col(1));
		VectorProduct(curr_basis.col(1), curr_basis.col(2), curr_basis.col(0));
		curr_basis.col(0) /= cv::norm(curr_basis.col(0));
		basis_map[*curr_vert] = curr_basis.inv(DECOMP_SVD);
	}
} 

template <class Graph, class BasisMap, class MaskMap>
void MarkUncoincidedBasises(const Graph& graph, const BasisMap& basis_map, MaskMap& mask_map)
{
	for (auto curr_vert = boost::vertices(graph).first, end_vert = boost::vertices(graph).second; curr_vert != end_vert; ++curr_vert)
	{
		double uncoins_vectors = 0;

		for (auto curr_neighb = adjacent_vertices(*curr_vert, graph).first, end_neighb = adjacent_vertices(*curr_vert, graph).second;
			curr_neighb != end_neighb; ++curr_neighb)
		{
			bool uncoin = false;

			for (int col_num = 0; col_num < 3; ++col_num)
			{
				if ((basis_map[*curr_vert].inv(DECOMP_SVD)).col(col_num).dot((basis_map[*curr_neighb].inv(DECOMP_SVD)).col(col_num)) < 0.0)
				{
					mask_map[*curr_vert] = 1;
					break;
				}
			}	
		}
		/*if (uncoins_vectors / out_degree(*curr_vert, graph) > 0.33)
		{
			++alarm;
			bad_vertices_map[*curr_vert] = 1;
		}
		else
		{
			bad_vertices_map[*curr_vert] = 0;
		}*/
	}
}
}