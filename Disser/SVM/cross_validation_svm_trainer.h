#pragma once
#include <fstream>
#include "shark/Models/Kernels/AbstractKernelFunction.h"
#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include <shark/Data/CVDatasetTools.h>
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>
#include "modified_grid_search.h"

namespace molecule_descriptor
{

template <typename InputTypeT>
class CrossValidationSvmTrainer
{
public:
	CrossValidationSvmTrainer()
		: m_kernel(nullptr),
		m_svm_trainer(nullptr, 0.001, true, false)
	{}

	void SetOutput(const std::string& out_name) {
		m_out.reset(new std::ofstream(out_name));
		//SHARK_ASSERT(m_out->is_open());
		//m_grid.SetOutput(*m_out);
	}
	template<typename KernelType>
	void SetKernel(KernelType kernel) {
		m_kernel.reset(new KernelType(kernel));
		m_svm_trainer.setKernel(m_kernel.get());
	}
	void SetGrid(const std::vector<double>& min, const std::vector<double>& max, const std::vector<size_t>& sections) {
		m_grid.configure(min, max, sections);
	}
	//returns training error
	double train(shark::KernelExpansion<InputTypeT>& model, shark::LabeledData<InputTypeT, unsigned int>& train_set)
	{
		shark::CVFolds<shark::LabeledData<InputTypeT, unsigned int>> folds = 
			shark::createCVSameSizeBalanced(train_set, train_set.numberOfElements());
		shark::CrossValidationError<shark::KernelExpansion<InputTypeT>, unsigned int> cv_error(
			folds, &m_svm_trainer, &model, &m_svm_trainer, &m_loss);
		m_grid.step(cv_error);
		//evaluate trained
		m_svm_trainer.setParameterVector(m_grid.solution().point);
		m_svm_trainer.train(model, train_set);
		m_eval_result = model(train_set.inputs());
		return m_loss.eval(train_set.labels(), m_eval_result);
	}
private:
	CrossValidationSvmTrainer(const CrossValidationSvmTrainer<InputTypeT>&);
	CrossValidationSvmTrainer<InputTypeT> operator=(const CrossValidationSvmTrainer<InputTypeT>&);

	std::unique_ptr<shark::AbstractKernelFunction<InputTypeT>> m_kernel;
	std::unique_ptr<std::ofstream> m_out;
	shark::CSvmTrainer<InputTypeT> m_svm_trainer;
	shark::ZeroOneLoss<unsigned int, shark::RealVector> m_loss;
	/*GridSearchFromEnd*/shark::GridSearch m_grid;
	shark::Data<shark::blas::vector<double>> m_eval_result;
};

}