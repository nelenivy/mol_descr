#pragma once
#include "stdafx.h"
#include "shark\Core\IParameterizable.h"
#include "../Common/singular_point.h"
#include "cross_validation_svm_trainer.h"
#include "Kernels/pharmacophore_kernel_calc_naive.h"
#include "zero_one_loss_modified.h"
#include "types_serialization.h"

namespace molecule_descriptor
{

template<typename PropType, typename PropKernel, typename DistKernel>
struct SingularPointsPairKernel : public shark::IParameterizable
{
	typedef SingularPointsPair<PropType> elem_type; 
	SingularPointsPairKernel(){	}
	SingularPointsPairKernel(const PropKernel prop_kernel, const DistKernel dist_kernel)
		: m_prop_kernel(prop_kernel),
		m_dist_kernel(dist_kernel)
	{
		CalcParamsVector();
	}
	inline double operator()(const elem_type& x, const elem_type& y)
	{
		const double res_1 = m_prop_kernel(x.elem1, y.elem1) * m_prop_kernel(x.elem2, y.elem2);
		const double res_2 = m_prop_kernel(x.elem1, y.elem2) * m_prop_kernel(x.elem2, y.elem1);
		const double res = std::max(res_1, res_2);
		double dist_res = 0.0;
		res > 0.0 && (dist_res = m_dist_kernel(x.dist, y.dist));
		return res * dist_res;
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

		init(new_parameters) >> shark::blas::parameters(m_prop_kernel), shark::blas::parameters(m_dist_kernel);
		CalcParamsVector();
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_prop_kernel.numberOfParameters() + m_dist_kernel.numberOfParameters());
		init(m_params_vect) << shark::blas::parameters(m_prop_kernel), shark::blas::parameters(m_dist_kernel);
	}
private:
	PropKernel m_prop_kernel;
	DistKernel m_dist_kernel;
	shark::RealVector m_params_vect;
};

template <typename PropT, typename PropKernel, typename DistKernel>
CrossValidationSvmTrainer<std::vector<SingularPointsPair<PropT>>, unsigned int, 
	shark::KernelClassifier<typename ElemTypeForTrainer<std::vector<SingularPointsPair<PropT>>, true>::type>,true> 
	GetNaiveKernelSVMTrainerForClassification(const PropKernel& prop_kernel, const DistKernel& dist_kernel)
{
	typedef SingularPointsPair<PropT> SingPtsPair;
	typedef std::vector<SingPtsPair> sing_pts_seq;
	typedef CrossValidationSvmTrainer<std::vector<SingularPointsPair<PropT>>, unsigned int, 
		shark::KernelClassifier<ElemWithIndexAndID<std::vector<SingularPointsPair<PropT>>>>,true> TrainerType;
	typedef typename TrainerType::ElemType ElemType;
	typedef SingularPointsPairKernel<PropT, PropKernel, DistKernel> PairKernel;
	typedef PharmSequenceKernelNaive<SingPtsPair, PairKernel> TrSetKernel;
	PairKernel pair_kernel(prop_kernel, dist_kernel);
	TrSetKernel set_kernel(pair_kernel);
	TrainerType cv_trainer;
	cv_trainer.SetTrainer(shark::CSvmTrainer<ElemType>(nullptr, 0.001, true, false));
	cv_trainer.SetKernel(set_kernel);
	cv_trainer.SetLoss(shark::ZeroOneLossNonAverage<unsigned int, unsigned int>());
	cv_trainer.SetModel(shark::KernelClassifier<ElemType>());
	return cv_trainer;
}

}

