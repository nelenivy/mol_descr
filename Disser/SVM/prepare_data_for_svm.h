#pragma once

#include <iterator>
#include <vector>
#include <type_traits>
#include <algorithm>

//#include <boost/range/iterator_range.hpp>
#include <shark/Data/CVDatasetTools.h>

namespace molecule_descriptor
{

	//delete conflicts from data. Conflicts are elements,which have equal data, but different labels
//conflicts are filled by existing elements
class DataConflictsDeleter
{
public:
	template <typename T>
	bool DeleteConflict(shark::LabeledData<T, unsigned int>& data_in);
private:
	std::vector<size_t> m_conflict_positions;
	std::vector<unsigned char> m_is_conflict;
	std::vector<size_t> m_positive_class_elems;
	std::vector<size_t> m_negative_class_elems;
};

enum Labels {kNegativeClass = 0, kPositiveClass = 1};

template <typename DataIterT, typename LabelIterT>
shark::LabeledData<typename std::iterator_traits<DataIterT>::value_type, unsigned int>
	PrepareDataForSVM(const DataIterT data_begin, const DataIterT data_end, 
	const LabelIterT labels_begin, const LabelIterT labels_end)
{
	//check
	typedef typename std::iterator_traits<LabelIterT>::value_type LabelType;
	static_assert(std::is_same<LabelType, unsigned int>::value, "label type must be unsigned int");
	typedef typename std::iterator_traits<DataIterT>::value_type DataType;
	const int data_size = data_end - data_begin;
	SHARK_ASSERT(data_size == labels_end - labels_begin);

	int positives = 0, negatives = 0;
	for (auto iter = labels_begin; iter != labels_end; ++iter)	{
		SHARK_ASSERT(*iter == kNegativeClass || *iter == kPositiveClass);
		*iter == 1 ? positives += 1 : negatives += 1;
	}
	SHARK_ASSERT(positives > 0 && negatives > 0);
	std::vector<DataType> input_vector;
	input_vector.reserve(2 * data_size);
	input_vector.assign(data_begin, data_end);
	std::vector<unsigned int> labels_vector;
	labels_vector.reserve(2 * data_size);
	labels_vector.assign(labels_begin, labels_end);
	//if there are more objects of one class, add them
	if (positives != negatives)
	{
		const int elems_to_add = std::abs(positives - negatives);
		const Labels elems_to_process = positives > negatives ? kNegativeClass : kPositiveClass;
		const int less_class_size = positives > negatives ? negatives : positives;
		std::vector<size_t> less_class_positions;
		less_class_positions.reserve(less_class_size);

		for (auto iter = labels_vector.begin(); iter != labels_vector.end(); ++iter) {
			if (*iter == elems_to_process) {
				less_class_positions.push_back(iter - labels_vector.begin());
			}
		}

		for (int curr_elem = 0; curr_elem < elems_to_add; ++curr_elem) 
		{
			const int ind_to_add = less_class_positions[curr_elem % less_class_size];
			labels_vector.push_back(labels_vector[ind_to_add]);
			input_vector.push_back(input_vector[ind_to_add]);
		}
	}

	shark::Data<DataType> input_data = shark::createDataFromRange(input_vector);
	shark::Data<unsigned int> input_labels = shark::createDataFromRange(labels_vector);
	shark::LabeledData<DataType, unsigned int> labeled_data(input_data, input_labels);
	return labeled_data;
}

template <typename VectorType>
bool AreVectorsSame(const VectorType& vector_1, const VectorType& vector_2)
{
	if (vector_1.size() != vector_2.size()) {
		return false;
	}

	bool are_same = true;
	const size_t size = vector_1.size();

	for (int ind = 0; ind < size; ++ind) {
		if (vector_1[ind] != vector_2[ind]) {
			are_same = false;
		}
	}

	return are_same;
}

template <typename T>
bool DataConflictsDeleter::DeleteConflict(shark::LabeledData<T, unsigned int>& data_in)
{
	if (data_in.numberOfElements() == 0)
	{
		return true;
	}

	const size_t elems_num = data_in.numberOfElements();
	m_conflict_positions.reserve(elems_num);
	m_conflict_positions.clear();
	m_positive_class_elems.reserve(elems_num);
	m_positive_class_elems.clear();
	m_negative_class_elems.reserve(elems_num);
	m_negative_class_elems.clear();
	m_is_conflict.resize(elems_num);
	std::fill(m_is_conflict.begin(), m_is_conflict.end(), 0);
	size_t conflicts_num = 0;

	//mark conflicts
	for (size_t ind_1 = 0; ind_1 < elems_num; ++ind_1)
	{
		if (m_is_conflict[ind_1]) {
			++conflicts_num;
			continue;
		}

		size_t curr_conflicts = 0;

		for (size_t ind_2 = 0; ind_2 < elems_num; ++ind_2)
		{
			if (data_in.labels().element(ind_1) != data_in.labels().element(ind_2) &&
				AreVectorsSame(data_in.inputs().element(ind_1), data_in.inputs().element(ind_2))) {
				m_is_conflict[ind_1] = 1;
				m_is_conflict[ind_2] = 1;
				++conflicts_num;
				break;
			}
		}

		if (!m_is_conflict[ind_1]) {
			if (data_in.labels().element(ind_1) == kPositiveClass) {
				m_positive_class_elems.push_back(ind_1);
			} else {
				m_negative_class_elems.push_back(ind_1);
			}
		}
	}

	if (conflicts_num >= elems_num || m_positive_class_elems.size() == 0 || m_negative_class_elems.size() == 0)
	{
		return false;
	}

	if (conflicts_num == 0) {
		return true;
	}
	//fill conflicts
	size_t conflict_ind = 0;
	for (size_t ind = 0; ind < elems_num; ++ind)
	{
		if (!m_is_conflict[ind]) {
			continue;
		}
		
		const std::vector<size_t>& class_to_borrow_elems = conflict_ind % 2 ? m_positive_class_elems : m_negative_class_elems;
		const size_t index_in_vect = (conflict_ind / 2) % class_to_borrow_elems.size();
		const size_t index_to_borrow = class_to_borrow_elems[index_in_vect];
		data_in.labels().element(ind) = data_in.labels().element(index_to_borrow);
		data_in.inputs().element(ind).assign(data_in.inputs().element(index_to_borrow));//assume that elements are vectors
		
		++conflict_ind;
	}

	return true;
}
}