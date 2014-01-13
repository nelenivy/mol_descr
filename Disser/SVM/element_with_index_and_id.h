#pragma once
#include <tuple>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <type_traits>

namespace molecule_descriptor
{

template <class ElemType, typename DataIDType = int>
class ElemWithIndexAndID
{
public:
	ElemWithIndexAndID() {

	}
	ElemWithIndexAndID(const ElemType& elem, const int ind, const DataIDType data_id, const size_t set_size)
		: m_elem_with_ind(elem, ind, data_id, set_size)
	{

	}
	ElemType& Elem()
	{
		return std::get<0>(m_elem_with_ind);
	}

	const ElemType& ElemConst() const
	{
		return std::get<0>(m_elem_with_ind);
	}

	int Index() const
	{
		return std::get<1>(m_elem_with_ind);
	}

	DataIDType DataID() const//id of the set which element belongs to
	{
		return std::get<2>(m_elem_with_ind);
	}

	size_t SizeOfSet() const//id of the set which element belongs to
	{
		return std::get<3>(m_elem_with_ind);
	}
private:
	std::tuple<ElemType, size_t, DataIDType, size_t> m_elem_with_ind;
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
}