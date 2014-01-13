#pragma once
//kernel type for not recalculating kernel each time
#include <tuple>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <limits>
#include "..\element_with_index_and_id.h"

namespace molecule_descriptor
{

struct LexicographicLess
{
	typedef std::pair<int, int> elem_type;
	bool operator()(const elem_type& elem1, const elem_type& elem2) const
	{
		if (elem1.first < elem2.first)
		{
			return true;
		}
		else if (elem1.first > elem2.first)
		{
			return false;
		}
		else
		{
			if (elem1.second < elem2.second)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
};

template <class ElemType, typename DataIDType = int, typename KernelValType = double>
class GrammMatrixCacheMapBased
{
public:
	static_assert(std::is_floating_point<KernelValType>::value, "Value type of kernel must be floating point type");
	GrammMatrixCacheMapBased()
		: m_inited(false){
	}
	typedef ElemWithIndexAndID<ElemType, DataIDType> elem_type;
	void InitIfNeed(const DataIDType data_id);
	bool IsInited()	const {
		return m_inited;
	}
	bool TryGetCached(const elem_type& elem1, const elem_type& elem2, KernelValType& result) const;
	bool TryWriteCache(const elem_type& elem1, const elem_type& elem2, const KernelValType& result);
	void Clear();
private:
	bool IsFromCurrSet(const elem_type& elem) const;

	std::map<std::pair<int, int>, KernelValType, LexicographicLess> m_gramm_mat;
	DataIDType m_set_id;
	bool m_inited;
};

template <class ElemType, typename DataIDType, typename KernelValType>
void GrammMatrixCacheMapBased<ElemType, DataIDType, KernelValType>::InitIfNeed(const DataIDType data_id)
{
	if (m_inited && data_id == m_set_id)
	{
		return;
	}
	else
	{
		m_set_id = data_id;
		m_gramm_mat.clear();
		m_inited = true;
		return;
	}
}

template <class ElemType, typename DataIDType, typename KernelValType>
void GrammMatrixCacheMapBased<ElemType, DataIDType, KernelValType>::Clear()
{
	m_gramm_mat.clear();
	m_inited = false;
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMapBased<ElemType, DataIDType, KernelValType>::IsFromCurrSet(const elem_type& elem) const
{
	if (!m_inited)
	{
		return false;
	}

	if (elem.DataID() != m_set_id)
	{
		return false;
	}

	return true;
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMapBased<ElemType, DataIDType, KernelValType>::TryGetCached(const elem_type& elem1, const elem_type& elem2, KernelValType& result) const
{
	if (!IsFromCurrSet(elem1) || !IsFromCurrSet(elem2))
	{
		return false;
	}

	auto index = std::make_pair(elem1.Index(), elem2.Index());
	auto res = m_gramm_mat.find(index);

	if (res != m_gramm_mat.end())
	{
		result = (*res).second;
		//std::cout << "srabotalo" << "\n";

		return true;
	}
	else
	{
		//std::cout << "vse norm" << elem1.Index() << " " << elem2.Index() << "\n";
		return false;
	}	
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMapBased<ElemType, DataIDType, KernelValType>::TryWriteCache(const elem_type& elem1, const elem_type& elem2, const KernelValType& result)
{
	if (!IsFromCurrSet(elem1) || !IsFromCurrSet(elem2))
	{
		return false;
	}

	auto index1 = std::make_pair(elem1.Index(), elem2.Index());
	auto index2 = std::make_pair(elem2.Index(), elem1.Index());
	auto res = m_gramm_mat.find(index1);

	if (res != m_gramm_mat.end())
	{
		return false;
	}
	else
	{
		m_gramm_mat[index1] = m_gramm_mat[index2] = result;
		return true;
	}	
}

template <class ElemType, typename DataIDType = int, typename KernelValType = double>
class GrammMatrixCacheMatBased
{
public:
	static_assert(std::is_floating_point<KernelValType>::value, "Value type of kernel must be floating point type");
	GrammMatrixCacheMatBased()
		: m_inited(false){
	}
	typedef ElemWithIndexAndID<ElemType, DataIDType> elem_type;
	void InitIfNeed(const DataIDType data_id, const size_t size);
	bool IsInited()	const {
		return m_inited;
	}
	bool TryGetCached(const elem_type& elem1, const elem_type& elem2, KernelValType& result) const;
	bool TryWriteCache(const elem_type& elem1, const elem_type& elem2, const KernelValType& result);
	void Clear();
private:
	bool IsFromCurrSet(const elem_type& elem) const;

	std::vector<KernelValType> m_gramm_mat;
	std::vector<char> m_gramm_mat_is_filled;
	size_t m_set_size;
	DataIDType m_set_id;
	bool m_inited;
};

template <class ElemType, typename DataIDType, typename KernelValType>
void GrammMatrixCacheMatBased<ElemType, DataIDType, KernelValType>::InitIfNeed(const DataIDType data_id, const size_t size)
{
	if (m_inited && data_id == m_set_id)
	{
		return;
	}
	else
	{
		m_set_id = data_id;
		m_set_size = size;
		m_gramm_mat.swap(std::vector<KernelValType>(size * size, std::numeric_limits<KernelValType>::min()));
		m_gramm_mat_is_filled.swap(std::vector<char>(size * size, 0));
		m_inited = true;
		return;
	}
}

template <class ElemType, typename DataIDType, typename KernelValType>
void GrammMatrixCacheMatBased<ElemType, DataIDType, KernelValType>::Clear()
{
	m_gramm_mat.clear();
	m_gramm_mat_is_filled.clear();
	m_inited = false;
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMatBased<ElemType, DataIDType, KernelValType>::IsFromCurrSet(const elem_type& elem) const
{
	if (!m_inited)
	{
		return false;
	}

	if (elem.DataID() != m_set_id)
	{
		return false;
	}

	return true;
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMatBased<ElemType, DataIDType, KernelValType>::TryGetCached(const elem_type& elem1, const elem_type& elem2, KernelValType& result) const
{
	if (!IsFromCurrSet(elem1) || !IsFromCurrSet(elem2))
	{
		return false;
	}


	const size_t index = elem1.Index() * m_set_size + elem2.Index();

	if (m_gramm_mat_is_filled[index] > 0)
	{
		result = m_gramm_mat[index];
		//std::cout << "srabotalo" << "\n";

		return true;
	}
	else
	{
		//std::cout << "vse norm" << elem1.Index() << " " << elem2.Index() << "\n";
		return false;
	}	
}

template <class ElemType, typename DataIDType, typename KernelValType>
bool GrammMatrixCacheMatBased<ElemType, DataIDType, KernelValType>::TryWriteCache(const elem_type& elem1, const elem_type& elem2, const KernelValType& result)
{
	if (!IsFromCurrSet(elem1) || !IsFromCurrSet(elem2))
	{
		return false;
	}

	const size_t index_1 = elem1.Index() * m_set_size + elem2.Index();
	const size_t index_2 = elem2.Index() * m_set_size + elem1.Index();

	if (m_gramm_mat_is_filled[index_1] > 0)
	{
		return false;
	}
	else
	{
		m_gramm_mat_is_filled[index_1] = m_gramm_mat_is_filled[index_2] = 1;
		m_gramm_mat[index_1] = m_gramm_mat[index_2] = result;
		return true;
	}	
}

}