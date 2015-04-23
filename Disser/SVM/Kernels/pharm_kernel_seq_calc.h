#pragma once
#include "stdafx.h"
#include <utility>
#include <string>
#include <vector>
#include <algorithm>

#include "opencv2/core/core.hpp"

#include "shark/Models/Kernels/AbstractKernelFunction.h"

#include "../Common/singular_point.h"
#include "basic_kernels.h"
#include "pharmacophore_kernel_calculator.h"
#include "gramm_matrix_cache.h"

namespace molecule_descriptor
{

using  shark::RealVector;
//SingularPoint<> is a pair of 3d point and its PropType
//PropType is used for kernel PropKernel
template<typename PropType, class DistKernel, class PropKernel>//PropKernel must have default constructor
class PharmSequenceKernel : public shark::AbstractKernelFunction<
	ElemWithIndexAndID<//it's to not calculate gramm matrix each time
	std::vector<SingularPoint<PropType> > 
						> >
{
public:
	typedef ElemWithIndexAndID<//it's to not calculate gramm matrix each time
		std::vector<SingularPoint<PropType> > 
	> value_type;
	typedef std::vector<SingularPoint<PropType>> SingularPointsSequence;

	PharmSequenceKernel()
		: m_props_kernel(PropKernel()),
		m_dist_kernel(DistKernel())
	{

	}

	PharmSequenceKernel(DistKernel dist_kernel, PropKernel prop_kernel)
		: m_props_kernel(prop_kernel),
		m_dist_kernel(dist_kernel)
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

		init(new_parameters) >> shark::blas::parameters(m_dist_kernel), shark::blas::parameters(m_props_kernel);
		CalcParamsVector();
		m_triples_mat_cache.Clear();
		m_pairs_mat_cache.Clear();

	}	
	/*virtual void eval(BatchInputType const& batchX1, BatchInputType const& batchX2, shark::RealMatrix& result, shark::State& state) const
	{
		eval(batchX1, batchX2, result);
	}*/
		
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_props_kernel.numberOfParameters() + m_dist_kernel.numberOfParameters());
		init(m_params_vect) << shark::blas::parameters(m_dist_kernel), shark::blas::parameters(m_props_kernel);
	}
	void EvalImplBatch(ConstBatchInputReference batchX1,	ConstBatchInputReference batchX2,
		shark::RealMatrix& result,	shark::State& state, const bool process_pairs) const;
	double EvalImplSingle(ConstInputReference elem_1,	ConstInputReference elem_2, const bool process_pairs) const;
private:
	mutable PharmacophoreKernelBase<PropType> m_pharm_kernel_for_elem;
	mutable PropKernel m_props_kernel;
	mutable DistKernel m_dist_kernel;
	shark::RealVector m_params_vect;
	mutable GrammMatrixCacheMatBased<SingularPointsSequence> m_pairs_mat_cache;
	mutable GrammMatrixCacheMatBased<SingularPointsSequence> m_triples_mat_cache;

};

template<typename PropType, class DistKernel, class PropKernel>
class PharmSequenceKernelPairs : public PharmSequenceKernel<PropType, DistKernel, PropKernel>
{
public:
	PharmSequenceKernelPairs()
		: PharmSequenceKernel<PropType, DistKernel, PropKernel>()
	{

	}
	PharmSequenceKernelPairs(DistKernel dist_kernel, PropKernel prop_kernel)
		: PharmSequenceKernel<PropType, DistKernel, PropKernel>(dist_kernel, prop_kernel)
	{

	}
	virtual std::string name()  
	{
		return "PharmSequenceKernelPairs";
	}

	virtual void eval(ConstBatchInputReference batchX1,	ConstBatchInputReference batchX2,
		shark::RealMatrix& result,	shark::State& state) const  
	{
		const bool process_pairs = true;
		EvalImplBatch(batchX1, batchX2, result, state, process_pairs);
	}

	virtual double eval(ConstInputReference batchX1,	ConstInputReference batchX2) const  
	{
		const bool process_pairs = true;
		return EvalImplSingle(batchX1, batchX2, process_pairs);
	}
};

template<typename PropType, class DistKernel, class PropKernel>
class PharmSequenceKernelTriples : public PharmSequenceKernel<PropType, DistKernel, PropKernel>
{
public:
	PharmSequenceKernelTriples()
		: PharmSequenceKernel<PropType, DistKernel, PropKernel>()
	{
	}
	PharmSequenceKernelTriples(DistKernel dist_kernel, PropKernel prop_kernel)
		: PharmSequenceKernel<PropType, DistKernel, PropKernel>(dist_kernel, prop_kernel)
	{

	}

	virtual std::string name() 
	{
		return "PharmSequenceKernelTriples";
	}

	virtual void eval(ConstBatchInputReference batchX1,	ConstBatchInputReference batchX2,
		shark::RealMatrix& result,	shark::State& state) const  
	{
		const bool process_pairs = false;
		EvalImplBatch(batchX1, batchX2, result, state, process_pairs);
	}
	virtual double eval(ConstInputReference batchX1,	ConstInputReference batchX2) const   
	{
		const bool process_pairs = false;
		return EvalImplSingle(batchX1, batchX2, process_pairs);
	}
};

template<typename PropType, class DistKernel, class PropKernel>
double PharmSequenceKernel<PropType, DistKernel, PropKernel>::EvalImplSingle(
	ConstInputReference elem_1,	ConstInputReference elem_2, const bool process_pairs) const
{
	//TODO:
	//SPLIT TO 2 FUNCTIONS
	//try to read cached values
	m_pairs_mat_cache.InitIfNeed(elem_1.DataID(), elem_1.SizeOfSet());
	m_pairs_mat_cache.InitIfNeed(elem_2.DataID(), elem_2.SizeOfSet());
	double pairs_res = 0;
	bool success_pairs = m_pairs_mat_cache.TryGetCached(elem_1, elem_2, pairs_res);
	m_triples_mat_cache.InitIfNeed(elem_1.DataID(), elem_1.SizeOfSet() );
	m_triples_mat_cache.InitIfNeed(elem_2.DataID(), elem_2.SizeOfSet());
	double triples_res = 0;
	bool success_triples = false;
	if (!process_pairs) {
		success_triples = m_triples_mat_cache.TryGetCached(elem_1, elem_2, triples_res);
	}

	if ((process_pairs && !success_pairs) ||//processing pairs and no cache for pairs
		(!process_pairs && !success_triples)) //processing triples and no cache for triples
	{//calculate and write to cache
		auto sing_points_1 = elem_1.ElemConst();
		auto sing_points_2 = elem_2.ElemConst();

		m_pharm_kernel_for_elem.Init(boost::begin(sing_points_1), boost::end(sing_points_1),
			boost::begin(sing_points_2), boost::end(sing_points_2));

		pairs_res = m_pharm_kernel_for_elem.ProcessPairs(m_dist_kernel, m_props_kernel);
		m_pairs_mat_cache.TryWriteCache(elem_1, elem_2, pairs_res);

		if (!process_pairs)
		{
			triples_res = m_pharm_kernel_for_elem.ProcessTriples();
			m_triples_mat_cache.TryWriteCache(elem_1, elem_2, triples_res);
		}
	}

	if (process_pairs)
	{
		//std::cout << elem_1.Index() << " " << elem_2.Index() << " " << pairs_res << "\n";

		return pairs_res;
	}
	else
	{
		return triples_res;
	}
}

template<typename PropType, class DistKernel, class PropKernel>
void PharmSequenceKernel<PropType, DistKernel, PropKernel>::EvalImplBatch(
	ConstBatchInputReference pts_1, ConstBatchInputReference pts_2, 
	shark::RealMatrix& kernel_matrix, shark::State&, const bool process_pairs) const
{
	const size_t mol_num_1 = boost::size(pts_1);
	const size_t mol_num_2 = boost::size(pts_2);
	kernel_matrix.resize(mol_num_1, mol_num_2);

	for (int ind_1 = 0; ind_1 < mol_num_1; ind_1++)
	{
		for (int ind_2 = 0; ind_2 < mol_num_2; ind_2++)
		{
			kernel_matrix(ind_1, ind_2) = EvalImplSingle(pts_1[ind_1], pts_2[ind_2], process_pairs);
		}
	}
}

}