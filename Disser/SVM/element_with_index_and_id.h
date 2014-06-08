#pragma once
#include "stdafx.h"
#include <tuple>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <type_traits>

namespace molecule_descriptor
{

template <class ElemType, typename DataIDType = int>
class ElemWithIndexAndID : public ElemType
{
public:
	ElemWithIndexAndID() {

	}
	ElemWithIndexAndID(const ElemType& elem, const int ind, const DataIDType data_id, const size_t set_size)
		: ElemType(elem), m_elem_with_ind(ind, data_id, set_size)
	{

	}
	ElemType& Elem()
	{
		return *this;
	}

	const ElemType& ElemConst() const
	{
		return *this;
	}

	int Index() const
	{
		return std::get<0>(m_elem_with_ind);
	}

	DataIDType DataID() const//id of the set which element belongs to
	{
		return std::get<1>(m_elem_with_ind);
	}

	size_t SizeOfSet() const//id of the set which element belongs to
	{
		return std::get<2>(m_elem_with_ind);
	}
private:
	std::tuple<size_t, DataIDType, size_t> m_elem_with_ind;
};

template <typename InIterType, typename ElemType>
void CreateSeqWithIndexAndIdFromSeq(InIterType in_begin, InIterType in_end, int dataset_id, 
									std::vector<ElemWithIndexAndID<ElemType>>& output)
{
	static_assert(std::is_same<
		std::remove_const<
		std::iterator_traits<InIterType>::value_type
		>::type, 
		ElemType>::value, "incompatible element and iterator types");
	const size_t set_size = static_cast<size_t>(in_end - in_begin);
	output.reserve(set_size);

	for (int ind = 0; in_begin != in_end; ++ind, ++in_begin)
	{
		output.push_back(ElemWithIndexAndID<ElemType>(*in_begin, ind, dataset_id, set_size));
	}
}

template <typename ElemT, typename LabelT, typename DataIdT>
shark::LabeledData<ElemT, LabelT> CopyFromIndexedToNonIndexed(
	shark::LabeledData<ElemWithIndexAndID<ElemT, DataIdT>, LabelT>& labeled_data)
{
	typedef ElemWithIndexAndID<ElemT, DataIdT> ElemWithInd;

	shark::UnlabeledData<ElemT> unlabeled_data(labeled_data.numberOfBatches());
	for (size_t batch_ind = 0; batch_ind < labeled_data.numberOfBatches(); ++batch_ind)
	{
		auto& curr_batch_in = labeled_data.inputs().batch(batch_ind);
		auto& curr_batch_out = unlabeled_data.inputs().batch(batch_ind);
		curr_batch_out = shark::Batch<ElemT>::createBatch(shark::get(curr_batch_in, 0).Elem(), curr_batch_in.size());
		//shark::Batch<ElemT>::resize(curr_batch_out, curr_batch_in.size(), 0);

		for (size_t elem_ind = 0; elem_ind < curr_batch_in.size(); ++elem_ind)
			shark::get(curr_batch_out, elem_ind)=shark::get(curr_batch_in, elem_ind).Elem();
			//curr_batch_out.push_back(curr_batch_in.inputs().element(ind).Elem());
	}

	return shark::LabeledData<ElemT, LabelT>(unlabeled_data, labeled_data.labels());
}
}