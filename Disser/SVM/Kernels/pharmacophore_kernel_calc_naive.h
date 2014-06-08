/************************************************************************/
/* this is implementation of the                                                                      */
/************************************************************************/
#pragma once
#include "stdafx.h"
#include <vector>
#include <type_traits>
#include <iterator>
#include <numeric>
#include <algorithm>

#include <omp.h>
#include "shark/Models/Kernels/AbstractKernelFunction.h"
#include "element_with_index_and_id.h"


namespace molecule_descriptor
{

using  shark::RealVector;
//SingularPoint<> is a pair of 3d point and its PropType
//PropType is used for kernel PropKernel
template<typename SingPtsT, class KernelT>//KernelT must have default constructor
class PharmSequenceKernelNaive : public shark::AbstractKernelFunction<
	ElemWithIndexAndID<//it's to not calculate gramm matrix each time
	std::vector<SingPtsT> > >
{
public:
	typedef std::vector<SingPtsT> SingularPointsSequence;

	PharmSequenceKernelNaive()
		: m_kernel(KernelT())	{ }
	PharmSequenceKernelNaive(KernelT kernel)
		: m_kernel(kernel)
	{
		CalcParamsVector();
	}
	virtual shark::RealVector parameterVector() const
	{
		return m_params_vect;
	}	 
	/// set the vector of hyper-parameters
	virtual void setParameterVector(RealVector const& new_parameters)
	{
		SHARK_ASSERT(new_parameters.size() == m_params_vect.size());

		if (std::equal(m_params_vect.begin(), m_params_vect.end(), new_parameters.begin()))
		{
			return;
		}

		init(new_parameters) >> shark::blas::parameters(m_kernel);
		CalcParamsVector();
		m_mat_cache.Clear();
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
	virtual void eval(ConstBatchInputReference batchX1,	ConstBatchInputReference batchX2,
		shark::RealMatrix& result,	shark::State& state) const  
	{
		EvalImplBatch(batchX1, batchX2, result);
	}
	virtual double eval(ConstInputReference batchX1, ConstInputReference batchX2) const  
	{
		return EvalImplSingle(batchX1, batchX2);
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_kernel.numberOfParameters());
		init(m_params_vect) << shark::blas::parameters(m_kernel);
	}
	void EvalImplBatch(ConstBatchInputReference batchX1,	ConstBatchInputReference batchX2,
		shark::RealMatrix& result) const;
	double EvalImplSingle(ConstInputReference elem_1,ConstInputReference elem_2) const;
private:
	mutable KernelT m_kernel;
	shark::RealVector m_params_vect;
	mutable GrammMatrixCacheMatBased<SingularPointsSequence> m_mat_cache;
};

template <typename IteratorType, typename SingPtsKernelT>
double CalcPharmKernelNaive(IteratorType points_1_begin, IteratorType points_1_end,
							IteratorType points_2_begin, IteratorType points_2_end, SingPtsKernelT& kernel)
{
	typedef std::remove_const<typename std::iterator_traits<IteratorType>::value_type>::type SingPtsT;
	const int64_t size_1 = points_1_end - points_1_begin;
	SHARK_ASSERT(size_1 >= 0);
	const size_t kThreadsNum = 7;
	std::array<double, kThreadsNum> res_arr = {};
	omp_set_num_threads(kThreadsNum);
#pragma omp parallel for
	for (int64_t ind_1 = 0; ind_1 < size_1; ++ind_1)
	{
		const int curr_thread = omp_get_thread_num();
		IteratorType iter_1 = points_1_begin + ind_1;

		for (IteratorType iter_2 = points_2_begin; iter_2 != points_2_end; ++iter_2)
		{
			res_arr[curr_thread] += kernel(*iter_1, *iter_2);
		}
	}
	return std::accumulate(res_arr.begin(), res_arr.end(), 0.0);
}

template<typename SingPtsT, class KernelT>
double PharmSequenceKernelNaive<SingPtsT, KernelT>::EvalImplSingle(
	ConstInputReference elem_1,	ConstInputReference elem_2) const
{
	//try to read cached values
	m_mat_cache.InitIfNeed(elem_1.DataID(), elem_1.SizeOfSet());
	m_mat_cache.InitIfNeed(elem_2.DataID(), elem_2.SizeOfSet());
	double res = 0;
	bool success = m_mat_cache.TryGetCached(elem_1, elem_2, res);
	if (!success) //no cache
	{//calculate and write to cache
		auto sing_points_1 = elem_1.ElemConst();
		auto sing_points_2 = elem_2.ElemConst();
		res = CalcPharmKernelNaive(boost::begin(sing_points_1), boost::end(sing_points_1),
			boost::begin(sing_points_2), boost::end(sing_points_2), m_kernel);
		m_mat_cache.TryWriteCache(elem_1, elem_2, res);
	}
	//std::cout << elem_1.Index() << " " << elem_2.Index() << " " << res << "\n";
	return res;
}

template<typename SingPtsT, class KernelT>
void PharmSequenceKernelNaive<SingPtsT, KernelT>::EvalImplBatch(
	ConstBatchInputReference pts_1, ConstBatchInputReference pts_2, 
	shark::RealMatrix& kernel_matrix) const
{
	const size_t mol_num_1 = boost::size(pts_1);
	const size_t mol_num_2 = boost::size(pts_2);
	kernel_matrix.resize(mol_num_1, mol_num_2);

	for (int ind_1 = 0; ind_1 < mol_num_1; ind_1++)
	{
		for (int ind_2 = 0; ind_2 < mol_num_2; ind_2++)
		{
			kernel_matrix(ind_1, ind_2) = EvalImplSingle(pts_1[ind_1], pts_2[ind_2]);
		}
	}
}




}