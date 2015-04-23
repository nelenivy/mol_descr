#include "stdafx.h"

#include "descriptor_svm_trainer.h"
#include <vector>
#include <fstream>
#include "boost/archive/polymorphic_text_oarchive.hpp"

#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/Models/Kernels/LinearKernel.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include <shark/Data/CVDatasetTools.h>
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>
//#include "shark/Data/Csv.h"
#include "element_with_index_and_id.h"
#include "modified_grid_search.h"
#include "prepare_data_for_svm.h"
#include "zero_one_loss_fuzzy.h"
#include "modified_cv_folds.h"

namespace molecule_descriptor
{	
	void DescriptorSVMTrainerManager::SetData(const cv::Mat_<size_t>& data, const vector<unsigned int>& labels)
	{
		//convert data to shark format
		std::vector<RealVector> temp(data.rows, RealVector(data.cols));
		for (int y = 0; y < data.rows; ++y)
		{
			for (int x = 0; x < data.cols; ++x)
			{
				temp[y][x] = data(y, x);
			}
		}

		std::vector<sing_pts_seq_wth_ind> data_wth_ind;

		CreateSeqWithIndexAndIdFromSeq(temp.begin(), temp.end(), m_curr_dataset_id, data_wth_ind);
		++m_curr_dataset_id;
		
		m_labeled_data_with_ind = 
			PrepareDataForSVM(data_wth_ind.begin(), data_wth_ind.end(), labels.begin(), labels.end());
		m_labeled_data = CopyFromIndexedToNonIndexed(m_labeled_data_with_ind);
		//shark::Data<sing_pts_seq_wth_ind> input_data = shark::createDataFromRange(temp/*_wth_ind*/);
		//shark::Data<unsigned int> input_labels = shark::createDataFromRange(labels);
		//const shark::LabeledData<sing_pts_seq_wth_ind, unsigned int> labeled_data(input_data, input_labels);
		//m_labeled_data = labeled_data;
	}

	/*template <typename PropType, class DistKernel, class PropKernel>
	void DescriptorSVMTrainerManager<PropType, DistKernel, PropKernel>::SetKernels(DistKernel dist_kernel, PropKernel prop_kernel)
	{
	PharmSequenceKernelPairs<PropType, DistKernel, PropKernel> pairs_kernel(dist_kernel, prop_kernel);
	m_pairs_kernel = pairs_kernel;
	}*/

	void DescriptorSVMTrainerManager::Train(const std::string& file_name)
	{
		//trainer
		const double C = 0.001;
		shark::LinearKernel<sing_pts_seq/*_wth_ind*/> linear_kernel;
		shark::CSvmTrainer<sing_pts_seq/*_wth_ind*/> svm_trainer(&linear_kernel, C, true, false);
		// cross-validation error
		const unsigned int N= 5;//m_labeled_data.elements().size();  // number of folds
		shark::ZeroOneLoss/*Fuzzy*/<unsigned int, unsigned int> loss/*(-1.0, 1.0)*/;
		shark::CVFolds<shark::LabeledData<sing_pts_seq/*_wth_ind*/, unsigned int>> folds = 
			shark::createCVSameSizeBalancedNoIndexedElems(m_labeled_data_with_ind, N);
		shark::CrossValidationError<shark::KernelClassifier<sing_pts_seq/*_wth_ind*/>, unsigned int> cv_error(
			folds, &svm_trainer, &m_kernel_expansion, &svm_trainer, &loss);
		// find best parameters
		std::vector<double> min(1);
		std::vector<double> max(1);
		std::vector<size_t> sections(1);
		min[0] = -8; max[0] = 2; sections[0] = 11;  // regularization parameter C
		GridSearchFromEnd grid;
		grid.configure(min, max, sections);
		std::ofstream out(file_name + "_res_simple.txt");
		grid.SetOutput(out);
		grid.step(cv_error);
		// train model on the full dataset
		svm_trainer.setParameterVector(grid.solution().point);
		std::cout << grid.solution().value << "\n";
		std::cout << grid.solution().point << "\n";
		svm_trainer.train(m_kernel_expansion, m_labeled_data/*_with_ind*/);

		EvaluateTrained();
	}



	void DescriptorSVMTrainerManager::EvaluateTrained()
	{
		/*m_eval_result = m_kernel_expansion(m_labeled_data.inputs());
		shark::ZeroOneLossFuzzy<unsigned int, shark::RealVector> loss(-1.0, 1.0);
		m_training_error = loss.eval(m_labeled_data.labels(), m_eval_result);*/
	}

	void DescriptorSVMTrainerManager::Write(const std::string& file_name)
	{
		std::ofstream out(file_name + "_res_normal_simple.txt");
		for (auto iter= m_eval_result.elements().begin(); iter != m_eval_result.elements().end() ; ++iter)
		{
			out << *iter << " ";
		}
		out.close();

		/*std::ofstream out(file_name + "_res_simple.txt");
		boost::archive::polymorphic_text_oarchive oa(out);
		m_eval_result.write(oa);
		out.close();

		
		out.open(file_name + "_labels_in_simple.txt");
		boost::archive::polymorphic_text_oarchive oa3(out);
		m_labeled_data.labels().write(oa);
		out.close();

		out.open(file_name + "_svm_simple.txt");
		boost::archive::polymorphic_text_oarchive oa1(out);
		m_kernel_expansion.write(oa1);
		out.close();
		out.open(file_name + "_res_set_simple.txt");
		out << m_training_error;*/
	}
}