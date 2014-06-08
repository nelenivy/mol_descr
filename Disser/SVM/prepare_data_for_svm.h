#pragma once
#include "stdafx.h"
#include <iterator>
#include <vector>
#include <type_traits>
#include <algorithm>

#include <shark/Data/CVDatasetTools.h>

namespace molecule_descriptor
{
//////////////////////////////////////////////////////////////////////////
//labels comparator
//////////////////////////////////////////////////////////////////////////
template<typename T>
struct LabelsCmp;

//delete conflicts from data. Conflicts are elements,which have equal data, but different labels
//conflicts are filled by existing elements
class DataConflictsDeleter
{
public:
	template <typename InputT, typename LabelT, typename LabelsComparator>
	bool DeleteConflict(shark::LabeledData<InputT, LabelT>& data_in, size_t& conflicts_num, LabelsComparator cmp = LabelsComparator());
private:
	std::vector<unsigned char> m_is_conflict;	
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
	std::cout << "positives " << positives <<  " negatives " << negatives << "\n";
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
//////////////////////////////////////////////////////////////////////////
template <typename VectorType>
bool AreVectorsSame(const VectorType& vector_1, const VectorType& vector_2)
{
	if (vector_1.size() != vector_2.size()) {
		return false;
	}

	bool are_same = true;
	const size_t size = vector_1.size();

	for (size_t ind = 0; ind < size; ++ind) {
		if (vector_1[ind] != vector_2[ind]) {
			are_same = false;
			break;
		}
	}

	return are_same;
}
//////////////////////////////////////////////////////////////////////////
template<>
struct LabelsCmp<unsigned int>
{
	typedef unsigned int value_type;

	bool operator()(const value_type elem1, const value_type elem2)
	{
		return elem1 == elem2;
	}
};
template <typename T>
T Sqr(const T num)
{
	return num * num;
}
template<>
struct LabelsCmp<shark::RealVector>
{
	typedef shark::RealVector value_type;
	explicit LabelsCmp(const double thresh)
		: m_thresh(thresh)
	{	}
	bool operator()(const value_type& elem1, const value_type& elem2)
	{
		CV_Assert(elem1.size() == elem2.size());
		double diff = 0;
		for (size_t ind = 0; ind < elem1.size(); ++ind)
		{
			diff += Sqr(elem1[ind] - elem2[ind]);
		} 
		diff = sqrt(diff);

		return diff < m_thresh;
	}
	double m_thresh;
};

template <typename InputT, typename LabelT>
LabelsCmp<LabelT> GetLabelsCmp(const shark::LabeledData<InputT, LabelT>& data_in)
{
	return LabelsCmp<LabelT>();
}
template <typename InputT>
LabelsCmp<shark::RealVector> GetLabelsCmp(const shark::LabeledData<InputT, shark::RealVector>& data_in)
{
	const double thresh = 1e-10;
	return LabelsCmp<shark::RealVector>(thresh);
}

//////////////////////////////////////////////////////////////////////////
template <typename InputT, typename LabelT, typename LabelsComparator>
bool DataConflictsDeleter::DeleteConflict(shark::LabeledData<InputT, LabelT>& data_in, size_t& conflicts_num, LabelsComparator cmp)
{
	if (data_in.numberOfElements() == 0)
	{
		return true;
	}

	const size_t elems_num = data_in.numberOfElements();
	
	m_is_conflict.resize(elems_num);
	std::fill(m_is_conflict.begin(), m_is_conflict.end(), 0);
	conflicts_num = 0;

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
			if (!cmp(data_in.labels().element(ind_1), data_in.labels().element(ind_2)) &&
				AreVectorsSame(data_in.inputs().element(ind_1), data_in.inputs().element(ind_2))) {
				m_is_conflict[ind_1] = 1;
				m_is_conflict[ind_2] = 1;
				++conflicts_num;
				break;
			}
		}		
	}

	if (conflicts_num == elems_num)
	{
		return false;
	}
	if (conflicts_num == 0) 
	{
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
	//found if there is only one class without conflicts
	bool are_same = true;
	{
		bool first = true;
		LabelT prev_label(data_in.labels().element(0));
		for (size_t ind = 0; ind < m_is_conflict.size() - 1; ++ind)
		{
			if (m_is_conflict[ind])
			{
				continue;
			}
			if (!first && !cmp(data_in.labels().element(ind), prev_label))
			{
				are_same = false;
				break;
			}
			first = false;
			prev_label = data_in.labels().element(ind);			
		}
	}

	if (are_same)
	{
		return false;
	}
	
	//delete conflicts
	std::vector<InputT> inputs_without_confl(elems_num - conflicts_num, data_in.inputs().element(0));
	std::vector<LabelT> labels_without_confl(elems_num - conflicts_num, data_in.labels().element(0));
	size_t non_conflict_ind = 0;
	for (size_t ind = 0; ind < elems_num; ++ind)
	{
		if (m_is_conflict[ind]) 
		{
			continue;
		}
		inputs_without_confl[non_conflict_ind] = data_in.inputs().element(ind);
		labels_without_confl[non_conflict_ind] = data_in.labels().element(ind);		
		++non_conflict_ind;
	}

	shark::Data<InputT> input_data = shark::createDataFromRange(inputs_without_confl);
	shark::Data<LabelT> input_labels = shark::createDataFromRange(labels_without_confl);
	shark::LabeledData<InputT, LabelT> new_labeled_data(input_data, input_labels);
	swap(new_labeled_data, data_in);
	return true;
}
}