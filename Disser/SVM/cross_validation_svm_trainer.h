#pragma once
#include "stdafx.h"
#include <fstream>
#include <utility>
#include "shark/Models/Kernels/AbstractKernelFunction.h"
#include "shark/Models/Kernels/LinearKernel.h"
#include "shark/Models/Kernels/NormalizedKernel.h"

#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/Algorithms/Trainers/EpsilonSvmTrainer.h"
#include "shark/ObjectiveFunctions/Loss/AbstractLoss.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include "shark/ObjectiveFunctions/Loss/AbsoluteLoss.h"
#include <shark/Data/CVDatasetTools.h>
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>
#include "modified_grid_search.h"
#include "modified_cv_folds.h"
#include "zero_one_loss_modified.h"
#include "cross_vaidation_modified.h"

namespace molecule_descriptor
{

template <typename InputTypeT, typename LabelT, bool kCacheKernel = false>
struct FoldsForSvmCreator;

template <typename InputTypeT, typename LabelT>
struct FoldsForSvmCreator<InputTypeT, LabelT, false>
{
	static shark::CVFolds<shark::LabeledData<InputTypeT, LabelT>> Do(
		shark::LabeledData<InputTypeT, LabelT>& train_set, const int folds_num)
	{
		 return shark::createCVSameSize(train_set, folds_num);
	}
};
//for classification
template <typename InputTypeT>
struct FoldsForSvmCreator<InputTypeT, unsigned int, true>
{
	typedef ElemWithIndexAndID<InputTypeT> InputIndT;
	static shark::CVFolds<shark::LabeledData<InputIndT, unsigned int>> Do(
		shark::LabeledData<InputTypeT, unsigned int>& train_set, const int folds_num)
	{
		std::vector<ElemWithIndexAndID<InputTypeT>> data_wth_ind;
		CreateSeqWithIndexAndIdFromSeq(train_set.inputs().elements().begin(), train_set.inputs().elements().end(), 0, data_wth_ind);
		shark::Data<ElemWithIndexAndID<InputTypeT>> input_data = shark::createDataFromRange(data_wth_ind);
		auto labeled_data_with_ind = //shark::LabeledData<ElemWithIndexAndID<InputTypeT>, unsigned int>(input_data, train_set.labels());
			PrepareDataForSVM(data_wth_ind.begin(), data_wth_ind.end(), train_set.labels().elements().begin(), train_set.labels().elements().end());
		size_t conflicts_num = 0;
		DataConflictsDeleter deleter;
		CV_Assert(deleter.DeleteConflict(labeled_data_with_ind, conflicts_num, GetLabelsCmp(labeled_data_with_ind)));
		return shark::createCVSameSizeBalancedIndexedElems(labeled_data_with_ind, folds_num);
	}
};
//
template <typename InputTypeT>
struct FoldsForSvmCreator<InputTypeT, unsigned int, false>
{
	static shark::CVFolds<shark::LabeledData<InputTypeT, unsigned int>> Do(
		shark::LabeledData<InputTypeT, unsigned int>& train_set, const int folds_num)
	{
		std::vector<ElemWithIndexAndID<InputTypeT>> data_wth_ind;
		CreateSeqWithIndexAndIdFromSeq(train_set.inputs().elements().begin(), train_set.inputs().elements().end(), 0, data_wth_ind);
		auto labeled_data_with_ind = 
			PrepareDataForSVM(data_wth_ind.begin(), data_wth_ind.end(), train_set.labels().elements().begin(), train_set.labels().elements().end());
		return shark::createCVSameSizeBalancedNoIndexedElems(labeled_data_with_ind, folds_num);
	}
};

template <typename InputTypeT, bool kUseInd>
struct ElemTypeForTrainer
{
	typedef typename std::conditional<kUseInd, ElemWithIndexAndID<InputTypeT>, InputTypeT>::type type;
};
template <typename InputTypeT, typename LabelT, class Model, bool kCacheKernel = false>
class CrossValidationSvmTrainer
{
public:
	typedef typename ElemTypeForTrainer<InputTypeT, kCacheKernel>::type ElemType;

	CrossValidationSvmTrainer()
		: m_kernel(nullptr),
		m_svm_trainer(nullptr),
		m_grid(nullptr)
	{}

	void SetOutput(const std::string& out_name) {
		m_out.reset(new std::ofstream(out_name));
		//SHARK_ASSERT(m_out->is_open());
		//m_grid.SetOutput(*m_out);
	}
	template<typename LossType>
	void SetLoss(const LossType loss)
	{
		m_loss.reset(new LossType(loss));
	}
	template<typename KernelType>
	void SetKernel(KernelType kernel) {
		m_kernel.reset(new KernelType(kernel));
	}
	void SetGrid(const std::vector<double>& min, const std::vector<double>& max, const std::vector<size_t>& sections) {
		m_grid.reset(new GridSearchType());
		m_grid->configure(min, max, sections);
	}
	template<typename TrainerType>
	void SetTrainer(const TrainerType& trainer)
	{
		m_svm_trainer.reset(new TrainerType(trainer));		
	}
	template<typename ModelType>
	void SetModel(const ModelType& model)
	{
		m_model.reset(new ModelType(model));
	}
	
//returns training error
	double train(shark::LabeledData<InputTypeT, LabelT>& train_set, const size_t folds_num)
	{
		CV_Assert(m_kernel != nullptr);
		CV_Assert(m_model != nullptr);
		CV_Assert(m_svm_trainer != nullptr);
		CV_Assert(m_loss != nullptr);
		m_folds = FoldsForSvmCreator<InputTypeT, LabelT, kCacheKernel>::Do(train_set, folds_num);
		auto kernel = GetNormalizedKernel<ElemType, false/*kCacheKernel*/>::Get(m_kernel);	

		shark::AbstractSvmTrainer<ElemType, LabelT, Model>* svm_trainter_casted = 
			dynamic_cast<shark::AbstractSvmTrainer<ElemType, LabelT, Model>*>(m_svm_trainer.get());
		if (svm_trainter_casted)
		{
			svm_trainter_casted->setKernel(kernel.get());
		}

		shark::IParameterizable* meta_obj = dynamic_cast<shark::IParameterizable*>(m_svm_trainer.get());
		SHARK_ASSERT(meta_obj);
		shark::CrossValidationErrorElemAv<Model> cv_error(
			m_folds, meta_obj, m_model.get(), m_svm_trainer.get(), m_loss.get());
		m_grid->step(cv_error);
		//evaluate trained
		meta_obj->setParameterVector(m_grid->solution().point);
		return m_grid->solution().value;		
	}

	const Model& TrainedModel(shark::LabeledData<InputTypeT, LabelT>& train_set) 
	{
		m_svm_trainer->train(*m_model, train_set);
		return *m_model;
	}
	const shark::RealVector OptimalParams() const
	{
		return m_grid->solution().point;
	}
	CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>& operator=(CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>&& trainer)
	{
		m_kernel = std::move(trainer.m_kernel);
		m_svm_trainer = std::move(trainer.m_svm_trainer);
		m_model = std::move(trainer.m_model);
		m_loss = std::move(trainer.m_loss);
		m_out = std::move(trainer.m_out);
		//m_eval_result = std::move(trainer.m_eval_result);

		//m_grid = std::move(trainer.m_grid);
		return *this;
	}
	CrossValidationSvmTrainer(CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>&& trainer)
	{
		*this = std::forward<CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>>(trainer);		
	}
	std::shared_ptr<shark::AbstractKernelFunction<ElemType>> m_kernel;

private:
	CrossValidationSvmTrainer(const CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>&);
	CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel> operator=(const CrossValidationSvmTrainer<InputTypeT, LabelT, Model, kCacheKernel>&);

	double TrainImpl(shark::LabeledData<InputTypeT, LabelT>& train_set)
	{	
		if (kCacheKernel)
		m_svm_trainer->train(*m_model, train_set);
		m_eval_result = (*m_model)(train_set.inputs());
		return m_loss->eval(train_set.labels(), m_eval_result);
	}
	typedef typename std::conditional<kCacheKernel, GridSearchFromEnd, shark::GridSearch>::type GridSearchType;

	std::unique_ptr<std::ofstream> m_out;
	typedef shark::AbstractSvmTrainer<ElemType, LabelT, Model> SVMType;
	std::unique_ptr<shark::AbstractTrainer<Model,LabelT>> m_svm_trainer;
	shark::CVFolds<shark::LabeledData<ElemType, LabelT>> m_folds;
	std::unique_ptr<shark::AbstractLoss<LabelT, typename Model::OutputType>> m_loss;
	std::unique_ptr<Model> m_model;
	std::unique_ptr<GridSearchType> m_grid;
	//std::unique_ptr<shark::Data<shark::blas::vector<double>>> m_eval_result;
};

template <typename ElemType, bool kNormalizedKernel>
class GetNormalizedKernel;

template <typename ElemType>
class GetNormalizedKernel<ElemType, true>
{
public:
	static std::shared_ptr<shark::NormalizedKernel<ElemType>> Get(std::shared_ptr<shark::AbstractKernelFunction<ElemType>>& kernel)
	{
		return std::shared_ptr<shark::NormalizedKernel<ElemType>>(new shark::NormalizedKernel<ElemType>(kernel.get()));
	}
};
template <typename ElemType>
class GetNormalizedKernel<ElemType, false>
{
public:
	static std::shared_ptr<shark::AbstractKernelFunction<ElemType>> Get(std::shared_ptr<shark::AbstractKernelFunction<ElemType>>& kernel)
	{
		return kernel;
	}
};


template <typename InputT>
CrossValidationSvmTrainer<InputT, unsigned int, shark::KernelClassifier<InputT>, false> GetSVMTrainerForClassification()
{
	CrossValidationSvmTrainer<InputT, unsigned int, shark::KernelClassifier<InputT>, false> cv_trainer;
	cv_trainer.SetTrainer(shark::CSvmTrainer<InputT>(nullptr, 0.001, true, false));
	cv_trainer.SetKernel(shark::LinearKernel<InputT>());
	cv_trainer.SetLoss(shark::ZeroOneLossNonAverage<unsigned int,unsigned int>());
	cv_trainer.SetModel(shark::KernelClassifier<InputT>());
	return cv_trainer;
}

template <typename InputT>
CrossValidationSvmTrainer<InputT, unsigned int, shark::LinearClassifier<InputT>, false> GetSVMTrainerForLinearClassification()
{
	CrossValidationSvmTrainer<InputT, unsigned int, shark::LinearClassifier<InputT>, false> cv_trainer;
	cv_trainer.SetTrainer(shark::LinearCSvmTrainer<InputT>(0.001, true));
	cv_trainer.SetKernel(shark::LinearKernel<InputT>());
	cv_trainer.SetLoss(shark::ZeroOneLossNonAverage<unsigned int,unsigned int>());
	cv_trainer.SetModel(shark::LinearClassifier<InputT>());
	return cv_trainer;
}

template <typename InputT>
CrossValidationSvmTrainer<InputT, shark::RealVector, shark::KernelExpansion<InputT>, false> GetSVMTrainerForRegression(const double epsilon)
{
	CrossValidationSvmTrainer<InputT, shark::RealVector, shark::KernelExpansion<InputT>, false> cv_trainer;
	cv_trainer.SetTrainer(shark::EpsilonSvmTrainer<InputT>(nullptr, 0.001, epsilon, true));
	cv_trainer.SetKernel(shark::LinearKernel<InputT>());
	cv_trainer.SetLoss(shark::AbsoluteLoss<shark::RealVector>());
	cv_trainer.SetModel(shark::KernelExpansion<InputT>());
	return cv_trainer;
}
}