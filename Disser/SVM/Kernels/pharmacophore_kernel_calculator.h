/************************************************************************/
/* this is implementation of the                                                                      */
/************************************************************************/
#pragma once
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <type_traits>
#include <iterator>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "shark/Models/Kernels/AbstractKernelFunction.h"


namespace molecule_descriptor
{

using std::vector;
using std::pair;
using std::tuple;
using std::array;

//SingularPoint<> is a pair of 3d point and its Property
//for Property is used kernel PropsKernel in method ProcessPairs
template<typename PropType>
class PharmacophoreKernelBase
{
public:
	typedef SingularPoint<PropType> singular_point;
	typedef vector<singular_point> singular_points_cont;
	template <typename IteratorType>
	void Init(IteratorType points_1_begin, IteratorType points_1_end,
		IteratorType points_2_begin, IteratorType points_2_end);
	template <class DistanceKernel, class PropsKernel>
	double ProcessPairs(DistanceKernel distance_kernel, PropsKernel props_kernel);
	double ProcessTriples();
	double GetPairsKernelVal() { 
		return m_curr_pairs_kernel_val; }
	double GetTriplesKernelVal() { 
		return m_curr_triples_kernel_val; }
private:
	/*double CalcKernelAll();
	double CalcKernelMax();*/
	singular_points_cont m_pts_1;
	singular_points_cont m_pts_2;
	vector<tuple<int, int, float>> m_precalc_vector_for_triples;//index1, index2, kernel val
	vector<float> m_kernels_for_pairs;
	vector<float> m_precalc_matrix_for_triples;
	double m_curr_pairs_kernel_val;
	double m_curr_triples_kernel_val;
};


template<typename PropType>
template <typename IteratorType>
void PharmacophoreKernelBase<PropType>::Init(IteratorType points_1_begin, IteratorType points_1_end,
	IteratorType points_2_begin, IteratorType points_2_end)
{
	typedef std::remove_const<typename std::iterator_traits<IteratorType>::value_type>::type value_type;
	static_assert(std::is_same<value_type, singular_point>::value, "wrong_input");
	m_pts_1.clear();
	m_pts_2.clear();
	std::copy(points_1_begin, points_1_end, std::back_inserter(m_pts_1));
	std::copy(points_2_begin, points_2_end, std::back_inserter(m_pts_2));

	const size_t pairs_num = m_pts_1.size() * m_pts_2.size();
	m_precalc_vector_for_triples.resize(pairs_num);
	//m_kernels_for_pairs.resize(pairs_num * pairs_num);
	//m_precalc_matrix_for_triples.resize(pairs_num * pairs_num);
}


template<typename PropType>
template <class DistanceKernel, class PropsKernel>
double PharmacophoreKernelBase<PropType>::ProcessPairs(DistanceKernel distance_kernel, PropsKernel props_kernel)
{
	int curr_pair = 0;
	//calculate properties kernels
	vector<int> indexes_1(m_precalc_vector_for_triples.size());
	vector<int> indexes_2(m_precalc_vector_for_triples.size());
	vector<float> kernels(m_precalc_vector_for_triples.size());

	for (int ind_1 = 0; ind_1 < m_pts_1.size(); ++ind_1)
	{
		for (int ind_2 = 0; ind_2 < m_pts_2.size(); ++ind_2)
		{
			float res = static_cast<float>(props_kernel(m_pts_1[ind_1].Property(), m_pts_2[ind_2].Property()));

			if (res > 0.0f)
			{
				//m_precalc_vector_for_triples[curr_pair] = std::make_tuple(ind_1, ind_2, res);
				indexes_1[curr_pair] = ind_1;
				indexes_2[curr_pair] = ind_2;
				kernels[curr_pair] = res;
				curr_pair++;
			}
		}
	}

	const int pair_num = curr_pair;
	m_precalc_vector_for_triples.resize(pair_num);
	indexes_1.resize(pair_num);
	indexes_2.resize(pair_num);
	kernels.resize(pair_num);
	m_precalc_matrix_for_triples.resize(pair_num * pair_num);
	//std::fill(m_kernels_for_pairs.begin(), m_kernels_for_pairs.end(), 0);
	m_curr_pairs_kernel_val = 0;
	const int pairs_1_num = m_pts_1.size() * m_pts_1.size();
	const int pairs_2_num = m_pts_2.size() * m_pts_2.size();

	//fill matrices
	const int kMaxThreads = 7;
	double res[kMaxThreads];
	cv::Point3d temp_point[kMaxThreads];
	std::fill(res, res + kMaxThreads, 0);
	omp_set_num_threads(kMaxThreads);
#pragma omp parallel for
	for (int y = 0; y < pair_num - 1; y++)
	{
		const int ind_1_y = indexes_1[y];
		const int ind_2_y = indexes_2[y];
		const float kernel_y = kernels[y];
		const int curr_thread = omp_get_thread_num();

		int* ind_1_x_ptr = &indexes_1[y + 1];
		int* ind_2_x_ptr = &indexes_2[y + 1];
		float* kernel_x_ptr = &kernels[y + 1];

		for (int x = y + 1; x < pair_num; x++, ++ind_1_x_ptr, ++ind_2_x_ptr, ++kernel_x_ptr )
		{
			//const auto& pair_1 = m_precalc_vector_for_triples[y];
			//const auto& pair_2 = m_precalc_vector_for_triples[x];
			////pairs from the same array
			//from_1_arr = std::make_pair(std::get<0>(pair_1), std::get<0>(pair_2));
			//from_2_arr = std::make_pair(std::get<1>(pair_1), std::get<1>(pair_2));
			const int ind_1_x = *ind_1_x_ptr;
			const int ind_2_x = *ind_2_x_ptr;
			const float kernel_x = *kernel_x_ptr;

			if (ind_1_x == ind_1_y || ind_2_x == ind_2_y)
			{//if points are the same
				m_precalc_matrix_for_triples[y * pair_num + x] = 
				m_precalc_matrix_for_triples[x * pair_num + y] = 0.0f;
				continue;
			}
			//calculate distance kernel
			temp_point[curr_thread] = m_pts_1[ind_1_x].Coord() - 
				m_pts_1[ind_1_y].Coord();
			const double dist_from_1_arr = cv::norm(temp_point[curr_thread]
				);
			temp_point[curr_thread] = m_pts_2[ind_2_x].Coord() - 
				m_pts_2[ind_2_y].Coord();
			const double dist_from_2_arr = cv::norm(temp_point[curr_thread]);

			const float dist_kern_val = static_cast<float>(distance_kernel(dist_from_1_arr, dist_from_2_arr));
			//assign properties
			/*const float prop_kern_pair_1 = std::get<2>(pair_1);
			const float prop_kern_pair_2 = std::get<2>(pair_2);*/
			m_precalc_matrix_for_triples[y * pair_num + x] = kernel_y * dist_kern_val;
			m_precalc_matrix_for_triples[x * pair_num + y] = kernel_x * dist_kern_val;
			//assign pairs kernels
			/*const int pair_1_index = from_1_arr.first * m_pts_1.size() + from_1_arr.second;
			const int pair_2_index = from_2_arr.first * m_pts_2.size() + from_2_arr.second;
			const int pair_index = pairs_2_num * pair_1_index + pair_2_index;*/
			res[curr_thread] += (kernel_y) * dist_kern_val * (kernel_x);
		}

		
	}
	//calculate kernel
	/*if (max_similar_to)
	{
		return CalcKernelMax();
	}
	else*/
	/*{
		return CalcKernelAll();
	}*/

	m_curr_pairs_kernel_val = std::accumulate(res, res + kMaxThreads, 0.0);
	return m_curr_pairs_kernel_val;
}

template<typename PropType>
double PharmacophoreKernelBase<PropType>::ProcessTriples()
{
	const int pairs_num = static_cast<size_t>(m_precalc_vector_for_triples.size());	
	const int kMaxThreads = 16;
	double res[kMaxThreads];
	std::fill(res, res + kMaxThreads, 0);
	omp_set_num_threads(kMaxThreads);

#pragma omp parallel for
	for (int i = 0;	i < pairs_num; i++)
	{
		const int curr_thread = omp_get_thread_num();
		const int i_row = i * pairs_num;

		for (int j = 0; j < pairs_num; j++)
		{
			
			const float precalc_i_j = m_precalc_matrix_for_triples[i_row + j];
			const float* precalc_j_row = &m_precalc_matrix_for_triples[j * pairs_num];
			if (precalc_i_j == 0.0)
			{
				continue;
			}

			for (int k = 0; k < pairs_num; k++, precalc_j_row++)
			{
				const int k_row = k * pairs_num;

				res[curr_thread] +=	precalc_i_j *
									*precalc_j_row * 
									m_precalc_matrix_for_triples[k_row + i];
			}
		}
	}

	m_curr_triples_kernel_val = std::accumulate(res, res + kMaxThreads, 0.0);
	return m_curr_triples_kernel_val;
}


}