#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "opencv2/core/core.hpp"
#include "shark/Models/Kernels/KernelExpansion.h"

#include "cross_validation_svm_trainer.h"
//#include "shark/Data/Csv.h"
#include "element_with_index_and_id.h"
#include "prepare_data_for_svm.h"

namespace molecule_descriptor
{

struct DescriptorsSetWithError
{
	DescriptorsSetWithError()
		: error(1.0)
	{}
	DescriptorsSetWithError(size_t capacity, double error_)
		: error(error_)
	{
		descriptors.reserve(capacity);
	}
	std::vector<size_t> descriptors;
	double error;//from 0.0 to 1.0
};

class MGUATrainer
{
public:
	typedef shark::RealVector data_type;
	template<typename InputTypeT>
	void SetData(const cv::Mat_<InputTypeT>& data,
		const std::vector<unsigned int>& labels);//Convert input to shark format
	void SetParameters(size_t max_iterations, size_t descr_to_select_num);
	template <typename TrainerType>
	void Train(TrainerType& trainer, const std::string& file_name);
	void Write(const std::string& file_name);
private:
	size_t m_max_iterations;
	size_t m_descr_to_select_num;
	size_t m_descriptors_num;
	size_t m_curr_iteration;
	shark::LabeledData<data_type, unsigned int> m_labeled_data;
	shark::LabeledData<data_type, unsigned int> m_sparsed_descriptors_labeled_data;
	std::vector<DescriptorsSetWithError> m_base_descr_set;
	std::vector<DescriptorsSetWithError> m_curr_best_set;
	DataConflictsDeleter m_conflicts_deleter;
};

template <typename InputTypeT>
void MGUATrainer::SetData(const cv::Mat_<InputTypeT>& data,
							 const std::vector<unsigned int>& labels)
{
	//convert data to shark format
	std::vector<data_type> temp(data.rows, data_type(data.cols));

	for (int y = 0; y < data.rows; ++y)
	{
		for (int x = 0; x < data.cols; ++x)
		{
			temp[y][x] = data(y, x);
		}
	}
	m_descriptors_num = data.cols;
	m_labeled_data = PrepareDataForSVM(temp.begin(), temp.end(), labels.begin(), labels.end());	
	SHARK_ASSERT(DataConflictsDeleter().DeleteConflict(m_labeled_data));

	//init labeled data for sparsed descriptors
	temp.resize(m_labeled_data.numberOfElements(), data_type(data.cols));
	std::vector<unsigned int> temp1(m_labeled_data.numberOfElements());
	shark::Data<data_type> input_data_1 = shark::createDataFromRange(temp/*_wth_ind*/);
	shark::Data<unsigned int> input_labels_1 = shark::createDataFromRange(temp1);
	const shark::LabeledData<data_type, unsigned int> labeled_data_1(input_data_1, input_labels_1);
	m_sparsed_descriptors_labeled_data = labeled_data_1;

	/*for (int y = 0; y < m_labeled_data.numberOfElements(); ++y)
	{
		for (int x = 1; x < m_labeled_data.inputs().element(y).size(); ++x)
		{
			m_sparsed_descriptors_labeled_data.inputs().element(y)[x] = 0;
			std::cout << m_labeled_data.inputs().element(y)[x] << " ";
		}
		std::cout << "\n";

	}*/


}

struct DescriptorsSetWithErrorLess {
	bool operator()(const DescriptorsSetWithError& a, const DescriptorsSetWithError& b)	{
		return a.error < b.error;
	}
};

//it's assumed that memory for output was already allocated
//return - if there are any different elements
template <typename T>
bool SparseDescriptors(const shark::LabeledData<T, unsigned int>& data_in, 
						DataConflictsDeleter& conflicts_deleter,
					   std::vector<size_t> descr_to_select, shark::LabeledData<T, unsigned int>& data_out)
{
	if (data_in.numberOfElements() == 0)
	{
		return false;
	}
	if (data_out.numberOfElements() != data_in.numberOfElements()
		|| data_out.inputs().element(0).size() != descr_to_select.size())
	{//init if need
		std::vector<T> temp(data_out.numberOfElements(), T(descr_to_select.size()));
		shark::Data<T> new_inputs = shark::createDataFromRange(temp);
		const shark::LabeledData<T, unsigned int> labeled_data(new_inputs, data_out.labels());
		data_out = labeled_data;
	}

	const size_t elems_num = data_in.numberOfElements();
	size_t different_elems = 0;
	for (size_t ind = 0; ind < elems_num; ++ind)
	{
		data_out.labels().element(ind) = data_in.labels().element(ind);

		for (size_t curr_descr_num = 0; curr_descr_num < descr_to_select.size(); ++curr_descr_num)
		{
			const size_t curr_descr_to_select = descr_to_select[curr_descr_num];
			data_out.inputs().element(ind)[curr_descr_num] = data_in.inputs().element(ind)[curr_descr_to_select];

			if (ind > 0 && 
				data_out.inputs().element(ind - 1)[curr_descr_num] != data_out.inputs().element(ind)[curr_descr_num]) {
				++different_elems;
			}
			std::cout << data_out.inputs().element(ind)[curr_descr_num] << " ";
		}
		std::cout << "\n";

		auto& curr_data_out_vect = data_out.inputs().element(ind);
		std::fill(curr_data_out_vect.begin() + descr_to_select.size(), curr_data_out_vect.end(), 0);
	}

	if (!conflicts_deleter.DeleteConflict(data_out)) {
		return false;
	}
	std::cout << "=============================\n";
	return (different_elems > 0);
}

const double kBigError = 2.0;

void WriteSetToFile(const std::vector<DescriptorsSetWithError>& base_descr_set, const std::string& file_name);

template <typename TrainerType>
void MGUATrainer::Train(TrainerType& trainer, const std::string& file_prefix)
{
	shark::KernelExpansion<data_type> kernel_expansion(true);
	DescriptorsSetWithError curr_descr_set(m_max_iterations, kBigError);

	for (int curr_iteration = 0; curr_iteration < m_max_iterations; ++curr_iteration)
	{
		std::cout << curr_iteration << "\n";
		curr_descr_set.descriptors.resize(curr_iteration + 1);
		//assign big error values to new best set
		for (auto iter = m_curr_best_set.begin(); iter != m_curr_best_set.end(); ++iter)
		{
			iter->error = kBigError;
		}
		std::pair<size_t, double> max_error_and_pos_among_best(0, kBigError);

		size_t prev_set_size = m_descr_to_select_num;
		if (curr_iteration == 0){
			prev_set_size = 1;//if we are in 0th iteration there are no saved descriptors
		}

		for (size_t curr_descr_num = 0; curr_descr_num < m_descriptors_num; ++curr_descr_num)
		{
			for (size_t curr_best_set_num = 0; curr_best_set_num < prev_set_size; ++curr_best_set_num)
			{
				const vector<size_t>& prev_iter_best_set = m_base_descr_set[curr_best_set_num].descriptors;
				//find if current descriptor is already in the set
				//in that case we don't process it
				if (prev_iter_best_set.end() != std::find(prev_iter_best_set.begin(), prev_iter_best_set.end(), curr_descr_num)) {
					continue;
				}
				//add current descriptor to set
				curr_descr_set.descriptors[0] = curr_descr_num;
				std::copy(prev_iter_best_set.begin(), prev_iter_best_set.end(), curr_descr_set.descriptors.begin() + 1);

				if (!SparseDescriptors(m_labeled_data, m_conflicts_deleter, curr_descr_set.descriptors, m_sparsed_descriptors_labeled_data)) {
					continue;
				}
				//find error
				curr_descr_set.error = trainer.train(kernel_expansion, m_sparsed_descriptors_labeled_data);
				//find if it can be placed in the best set
				if (curr_descr_set.error < max_error_and_pos_among_best.second)
				{
					m_curr_best_set[max_error_and_pos_among_best.first] = curr_descr_set;
					auto iter_max = std::max_element(m_curr_best_set.begin(), m_curr_best_set.end(), DescriptorsSetWithErrorLess());
					max_error_and_pos_among_best.first = iter_max - m_curr_best_set.begin();
					max_error_and_pos_among_best.second = iter_max->error;
				}
			}
		}

		m_base_descr_set = m_curr_best_set;
		WriteSetToFile(m_base_descr_set, file_prefix + "_descr_set.txt");
	}
}

}