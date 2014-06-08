#pragma once
#include "stdafx.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

#include "opencv2/core/core.hpp"
#include "shark/Models/Kernels/KernelExpansion.h"

#include "cross_validation_svm_trainer.h"
#include "shark/Data/Libsvm.h"
#include "element_with_index_and_id.h"
#include "prepare_data_for_svm.h"

namespace molecule_descriptor
{

struct DescriptorsSetWithError
{
	DescriptorsSetWithError()
		: recall(0.0),
		precision(0.0)
	{}
	DescriptorsSetWithError(size_t capacity, double recall_, double precision_)
		: recall(recall_), precision(precision_)
	{
		descriptors.reserve(capacity);
	}
	double FMeasure() const
	{
		return 2 * recall * precision / (recall + precision + 0.000001) + 0.35 * precision;
	}
	std::vector<size_t> descriptors;
	double recall;
	double precision;//error without penalty
};

template <typename LabelT>
class MGUATrainer
{
public:
	typedef shark::RealVector data_type;
	template<typename InputTypeT>
	void SetData(const cv::Mat_<InputTypeT>& data,
		const std::vector<LabelT>& labels, const bool normalize);//Convert input to shark format
	void SetParameters(size_t max_iterations, size_t descr_to_select_num);
	template <typename TrainerType>
	void Train(std::vector<TrainerType>& trainer_set, const std::string& file_prefix, const bool use_only_good_descriptors);
	template <typename TrainerType>
	void TrainOneSet(TrainerType& trainer, const std::string& file_prefix);
	const std::vector<bool>& GetGoodDescriptors()
	{
		return m_good_descr;
	}
	void Write(const std::string& file_name);
private:
	size_t m_max_iterations;
	size_t m_descr_to_select_num;
	size_t m_descriptors_num;
	shark::LabeledData<data_type, LabelT> m_labeled_data;
	shark::LabeledData<data_type, LabelT> m_sparsed_descriptors_labeled_data;
	std::vector<DescriptorsSetWithError> m_base_descr_set;
	std::vector<DescriptorsSetWithError> m_curr_best_set;
	std::vector<bool> m_good_descr;
	DataConflictsDeleter m_conflicts_deleter;
};

template <typename InputType, typename LabelT> 
struct ExportLibSvm
{
	static void Do(shark::LabeledData<InputType, LabelT>& dataset, const std::string &fn, 
		bool dense=false, bool oneMinusOne = true, bool sortLabels = false, bool append = false){}
};
template <typename InputType> 
struct ExportLibSvm<InputType, unsigned int>
{
	static void Do(shark::LabeledData<InputType, unsigned int>& dataset, const std::string &fn, 
		bool dense=false, bool oneMinusOne = true, bool sortLabels = false, bool append = false)
	{
		shark::export_libsvm(dataset, fn, dense, oneMinusOne, sortLabels, append);
	}
};
template <typename LabelT>
template <typename InputTypeT>
void MGUATrainer<LabelT>::SetData(const cv::Mat_<InputTypeT>& data,
							 const std::vector<LabelT>& labels, const bool normalize)
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

	if (normalize)
	{
		for (int x = 0; x < data.cols; ++x)
		{
			double norm = 0;
			for (int y = 0; y < data.rows; ++y)
			{
				norm += data(y, x) * data(y, x);
			}
			norm = sqrt(norm);

			if (norm > 0.0)
			{
				for (int y = 0; y < data.rows; ++y)
				{
					temp[y][x] /= norm;
				}
			}
		}
	}
	m_descriptors_num = data.cols;
	shark::Data<data_type> input_data = shark::createDataFromRange(temp/*_wth_ind*/);
	shark::Data<LabelT> input_labels = shark::createDataFromRange(labels);
	m_labeled_data = shark::LabeledData<data_type, LabelT>(input_data, input_labels);
	ExportLibSvm<data_type, LabelT>::Do(m_labeled_data, "data.txt");
	/*size_t conflicts_num = 0;
	SHARK_ASSERT(DataConflictsDeleter().DeleteConflict(m_labeled_data, conflicts_num, GetLabelsCmp(m_labeled_data)));
	std::cout << conflicts_num << "\n";*/
	//init labeled data for sparsed descriptors
	temp.resize(m_labeled_data.numberOfElements(), data_type(data.cols));
	std::vector<LabelT> temp1(m_labeled_data.numberOfElements(), labels[0]);
	shark::Data<data_type> input_data_1 = shark::createDataFromRange(temp/*_wth_ind*/);
	shark::Data<LabelT> input_labels_1 = shark::createDataFromRange(temp1);
	const shark::LabeledData<data_type, LabelT> labeled_data_1(input_data_1, input_labels_1);
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
		return a.FMeasure() < b.FMeasure();
	}
};

struct DescriptorsSetWithErrorGreater {
	bool operator()(const DescriptorsSetWithError& a, const DescriptorsSetWithError& b)	{
		return a.FMeasure() > b.FMeasure();
	}
};
//it's assumed that memory for output was already allocated
//return - if there are any different elements
template <typename T, typename LabelT>
bool SparseDescriptors(const shark::LabeledData<T, LabelT>& data_in, const bool normalize,
						DataConflictsDeleter& conflicts_deleter,
					   std::vector<size_t> descr_to_select, shark::LabeledData<T, LabelT>& data_out)
{
	if (data_in.numberOfElements() == 0)
	{
		return false;
	}
	if (data_out.numberOfElements() != data_in.numberOfElements()
		|| data_out.inputs().element(0).size() != descr_to_select.size())
	{//init if need
		std::vector<T> temp(data_in.numberOfElements(), T(descr_to_select.size()));
		shark::Data<T> new_inputs = shark::createDataFromRange(temp);
		std::vector<LabelT> temp_labels(data_in.numberOfElements(), data_in.labels().element(0));
		shark::Data<LabelT> new_labels = shark::createDataFromRange(temp_labels);
		shark::LabeledData<T, LabelT> labeled_data(new_inputs, new_labels);
		swap(data_out,labeled_data);
	}

	const size_t elems_num = data_in.numberOfElements();
	size_t different_elems = 0;
	for (size_t ind = 0; ind < elems_num; ++ind)
	{
		//std::cout << " " << ind << " ";

		data_out.labels().element(ind) = data_in.labels().element(ind);
		double norm = 0;
		for (size_t curr_descr_num = 0; curr_descr_num < descr_to_select.size(); ++curr_descr_num)
		{
			const size_t curr_descr_to_select = descr_to_select[curr_descr_num];
			data_out.inputs().element(ind)[curr_descr_num] = data_in.inputs().element(ind)[curr_descr_to_select];
			norm += data_out.inputs().element(ind)[curr_descr_num] * data_out.inputs().element(ind)[curr_descr_num];
			if (ind > 0 && 
				data_out.inputs().element(ind - 1)[curr_descr_num] != data_out.inputs().element(ind)[curr_descr_num]) {
				++different_elems;
			}
			//std::cout << data_out.inputs().element(ind)[curr_descr_num] << " ";
		}
		norm = sqrt(norm);

		if (normalize && norm > 0.0)
		{
			for (size_t curr_descr_num = 0; curr_descr_num < descr_to_select.size(); ++curr_descr_num)
			{
				data_out.inputs().element(ind)[curr_descr_num] /= norm;
			}
		}
		//std::cout << "\n";

		auto& curr_data_out_vect = data_out.inputs().element(ind);
		std::fill(curr_data_out_vect.begin() + descr_to_select.size(), curr_data_out_vect.end(), 0);
	}

	//std::cout << "=============================\n";
	return (different_elems > 0);
}

const double kBigError = std::numeric_limits<float>::max();

void WriteSetToFile(const std::vector<DescriptorsSetWithError>& base_descr_set, const std::string& file_name);

template <typename LabelT>
template <typename TrainerType>
void MGUATrainer<LabelT>::Train(std::vector<TrainerType>& trainer_set, const std::string& file_prefix, const bool use_only_good_descriptors)
{
	const size_t threads_num = trainer_set.size();
	const double kGoodDescrThresh = 0.5;
	m_good_descr.resize(m_descriptors_num);
	std::fill(m_good_descr.begin(), m_good_descr.end(), false);
	std::vector<DescriptorsSetWithError> curr_descr_set_thread(threads_num, DescriptorsSetWithError(m_max_iterations, 0.0, 0.0));
	std::vector<std::vector<DescriptorsSetWithError>> curr_best_set_thread(threads_num, m_curr_best_set);
	std::vector<std::pair<size_t, double>> least_fmeasure_and_pos_among_best_thread(threads_num, std::pair<size_t, double>(0, 0));
	std::vector<shark::LabeledData<data_type, LabelT>> labeled_data_sparsed_threads(threads_num);
	std::vector<DataConflictsDeleter> conflicts_deleter_threads(threads_num);
	omp_set_num_threads(threads_num);

	for (int curr_iteration = 0; curr_iteration < m_max_iterations; ++curr_iteration)
	{
		std::cout << "Iteration num " << curr_iteration << "\n";
		for (size_t ind = 0; ind < curr_descr_set_thread.size(); ++ind)
		{
			curr_descr_set_thread[ind].descriptors.resize(curr_iteration + 1);
		}
		//assign big error values to new best set
		for (auto thread_iter = curr_best_set_thread.begin(); thread_iter != curr_best_set_thread.end(); ++thread_iter)
		{
			for (auto iter = thread_iter->begin(); iter != thread_iter->end(); ++iter)
			{
				iter->precision = 0;
				iter->recall = 0;
			}
		}
		std::fill(least_fmeasure_and_pos_among_best_thread.begin(), least_fmeasure_and_pos_among_best_thread.end(), std::pair<size_t, double>(0, 0));

		size_t prev_set_size = m_descr_to_select_num;
		if (curr_iteration == 0){
			prev_set_size = 1;//if we are in 0th iteration there are no saved descriptors
		}

		size_t positive_descr = 0;
#pragma omp parallel for
		for (int curr_descr_num = 0; curr_descr_num < m_descriptors_num; ++curr_descr_num)
		{
			std::cout << "curr_descr_num " << curr_descr_num << "\n";
			const int curr_thread_num = omp_get_thread_num();
			std::pair<size_t, double>& least_fmeasure_and_pos_among_best = least_fmeasure_and_pos_among_best_thread[curr_thread_num];
			DescriptorsSetWithError& curr_descr_set = curr_descr_set_thread[curr_thread_num];
			std::vector<DescriptorsSetWithError>& curr_best_set = curr_best_set_thread[curr_thread_num];
			shark::LabeledData<data_type, LabelT>& labeled_data_sparsed = labeled_data_sparsed_threads[curr_thread_num];
			DataConflictsDeleter& conflicts_deleter = conflicts_deleter_threads[curr_thread_num];
			auto& trainer = trainer_set[curr_thread_num];

			if (use_only_good_descriptors && curr_iteration > 0 && !m_good_descr[curr_descr_num])
			{
				continue;
			}
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

				if (!SparseDescriptors(m_labeled_data, false, conflicts_deleter, curr_descr_set.descriptors, labeled_data_sparsed)) {
					continue;
				}
				size_t conflicts_num = 0;
				if (!conflicts_deleter.DeleteConflict(labeled_data_sparsed, conflicts_num, 
					GetLabelsCmp(labeled_data_sparsed))) 
				{
					continue;
				}
				const int kNonConflictsThresh = 10;
				if (m_labeled_data.numberOfElements() - conflicts_num < kNonConflictsThresh)
				{
					if (curr_iteration == 0)
					{
						m_good_descr[curr_descr_num] = true;
					}
					std::cout << "conflicts_num "<< conflicts_num << "\n";
					continue;
				}
				//find error
				curr_descr_set.precision = 1.0 - trainer.train(labeled_data_sparsed, 5);
				curr_descr_set.recall = curr_descr_set.precision * static_cast<double>(m_labeled_data.numberOfElements() - conflicts_num) / m_labeled_data.numberOfElements();
				
				std::cout << "precision " << curr_descr_set.precision << " recall= " << curr_descr_set.recall
					<< " f_measure= " << curr_descr_set.FMeasure() 
					<< " best f_measure= "  << least_fmeasure_and_pos_among_best.second << "\n";

				if (curr_descr_set.precision >= 1.0 - kGoodDescrThresh)
				{
					++positive_descr;
					if (curr_iteration == 0)
					{
						m_good_descr[curr_descr_num] = true;
					}
				}
				//find if it can be placed in the best set
				if (curr_descr_set.FMeasure() > least_fmeasure_and_pos_among_best.second)
				{
					curr_best_set[least_fmeasure_and_pos_among_best.first] = curr_descr_set;
					auto iter_min = std::min_element(curr_best_set.begin(), curr_best_set.end(), DescriptorsSetWithErrorLess());
					least_fmeasure_and_pos_among_best.first = iter_min - curr_best_set.begin();
					least_fmeasure_and_pos_among_best.second = iter_min->FMeasure();
				}
			}
		}

		m_curr_best_set.clear();
		for (size_t ind = 0; ind < curr_best_set_thread.size(); ++ind)
		{
			m_curr_best_set.insert(m_curr_best_set.end(), curr_best_set_thread[ind].begin(), curr_best_set_thread[ind].end());
		}

		std::sort(m_curr_best_set.begin(), m_curr_best_set.end(), DescriptorsSetWithErrorGreater());
		m_base_descr_set.assign(m_curr_best_set.begin(), m_curr_best_set.begin() + m_descr_to_select_num);
		WriteSetToFile(m_base_descr_set, file_prefix + "_descr_set.txt");
		const double positive_rate = positive_descr / static_cast<double>(m_descriptors_num * prev_set_size);
		std::cout << "positive_rate = " << positive_rate << "\n";
		std::ofstream f(file_prefix + "_descr_set.txt", std::ofstream::app);
		f << "\n" << positive_rate << "\n";
	}
}

template <typename LabelT>
template <typename TrainerType>
void MGUATrainer<LabelT>::TrainOneSet(TrainerType& trainer, const std::string& file_prefix)
{
	/*{{174, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{878, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{1199, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{27, 323, 324, 122, 1197, 38, 37, 903, 134, 180}, 
	{9, 188, 27, 37, 321, 165, 185, 903, 134, 180},
	{170, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{110, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{167, 188, 27, 37, 321, 165, 185, 903, 134, 180},
	{33, 323, 324, 122, 1197, 38, 37, 903, 134, 180},
	{171, 162, 18, 168, 326, 164, 324, 903, 134, 180}};*/
	const size_t descriptors_arr[][12] = 
	{ {5785, 5646, 8831, 5803, 9871, 4869, 6219, 9741, 7934, 1303, 8936, 6257}
		//{272, 14, 128, 176, 872, 135}//,
	/*{1295, 157, 1009, 1048, 118, 38},
	{1038, 14, 128, 176, 872, 135},
	{901, 37, 271, 1048, 118, 38}, 
	{1157, 180, 1009, 1048, 118, 38},
	{185, 6, 116, 1048, 118, 38},
	{180, 157, 1009, 1048, 118, 38 },
	{415, 157, 1009, 1048, 118, 38},
	{1053, 14, 128, 176, 872, 135},
	{1295, 258, 128, 176, 872, 135}*/};
	const int kSetNum = sizeof(descriptors_arr) / sizeof(descriptors_arr[0]);
	const int kDescrNum = sizeof(descriptors_arr[0]) / sizeof(descriptors_arr[0][0]);


	for (size_t curr_set = 0; curr_set < kSetNum; ++curr_set)
	{
		std::vector<size_t> descriptors(descriptors_arr[curr_set], descriptors_arr[curr_set] + kDescrNum);
				
		if (!SparseDescriptors(m_labeled_data, false, m_conflicts_deleter, descriptors, m_sparsed_descriptors_labeled_data)) {
			return;
		}
		{std::ofstream t("m4.txt");
		for (size_t y = 0; y < m_sparsed_descriptors_labeled_data.numberOfElements(); ++y)
		{
			const auto& curr_mol = m_sparsed_descriptors_labeled_data.inputs().element(y);
			for (int ind = 0; ind < curr_mol.size(); ++ind)
			{
				t << curr_mol[ind] << " ";
			}
			t << "\n";
			for (size_t x = 0; x < m_sparsed_descriptors_labeled_data.numberOfElements(); ++x)
			{
				t << y << " " << x << " " <<shark::LinearKernel<data_type>()(m_sparsed_descriptors_labeled_data.inputs().element(y), m_sparsed_descriptors_labeled_data.inputs().element(x)) << " ";
				t << "\n";
			}
		}}
		size_t conflicts_num = 0;
		if (!m_conflicts_deleter.DeleteConflict(m_sparsed_descriptors_labeled_data, conflicts_num, 
			GetLabelsCmp(m_sparsed_descriptors_labeled_data))) 
		{
			return;
		}
		const int kNonConflictsThresh = 10;
		if (m_labeled_data.numberOfElements() - conflicts_num < kNonConflictsThresh)
		{
		
			std::cout << "conflicts_num "<< conflicts_num << "\n";
			return;
		}
		//find error
		const double precision = 1.0 - trainer.train(m_sparsed_descriptors_labeled_data, 5/*m_sparsed_descriptors_labeled_data.numberOfElements() / 2*/);
		const double recall = precision * static_cast<double>(m_labeled_data.numberOfElements() - conflicts_num) / m_labeled_data.numberOfElements();

		std::cout << "precision " << precision << " recall= " << recall
			<< " f_measure= " << 2 * precision * recall * (precision + recall);	
	}
}

template <typename LabelT>
void MGUATrainer<LabelT>::SetParameters(size_t max_iterations, size_t descr_to_select_num)
{
	m_max_iterations = max_iterations;
	m_descr_to_select_num = descr_to_select_num;

	//change descriptors set containers
	m_base_descr_set.resize(m_descr_to_select_num);
	m_curr_best_set.resize(m_descr_to_select_num);

	std::fill(m_base_descr_set.begin(), m_base_descr_set.end(), DescriptorsSetWithError(m_max_iterations, kBigError, kBigError));
	std::fill(m_curr_best_set.begin(), m_curr_best_set.end(), DescriptorsSetWithError(m_max_iterations, kBigError, kBigError));
}

}