#pragma once

#include <vector>
#include <tuple>

namespace molecule_descriptor
{
enum ThrowExcOrNot { eThrow = 0, eNotThrow};

//priority heap, which is based on array
template <class IteratorType, typename KeyType, class Comparator = std::less<KeyType>>
class PriorityQueueOnRange
{
public:
	typedef typename iterator_traits<IteratorType>::value_type value_type;
	typedef typename iterator_traits<IteratorType>::pointer pointer;
	typedef typename iterator_traits<IteratorType>::reference reference;

	class HeapElement
	{
	public:
		HeapElement(const std::tuple<KeyType, pointer, size_t>& tuple) : m_tuple(tuple) { }
		KeyType& Key() { return std::get<0>(m_tuple); }
		pointer Pointer() { return std::get<1>(m_tuple); }
		size_t Index() { return std::get<2>(m_tuple); }
	private:
		std::tuple<KeyType, pointer, size_t> m_tuple;
	};

	PriorityQueueOnRange(IteratorType begin, IteratorType end, const KeyType default_key);
	HeapElement Top() { return m_heap[0]; }
	void Pop();
	size_t Size() { return m_heap_size; }
	void ChangeKey(const size_t index, const KeyType new_key, const ThrowExcOrNot no_throw);/**< changes key of an element with index "index", if Comparator(new_key, old_key)*/
private:
	void SwapElements(const size_t heap_index_1, const size_t heap_index_2);
	size_t Heapify(const size_t heap_index);/**< returns the number of the smallest element among current and its children*/	

	std::vector<HeapElement> m_heap;
	std::vector<size_t> m_positions;//m_positions[i] stores position in m_heap of ith element of initial range
	size_t m_heap_size;
};

//создать очередь с элементами в порядке возрастания
template <class IteratorType, typename KeyType>
PriorityQueueOnRange<IteratorType, KeyType, std::less<KeyType>> MakeRisingProrityQueue(IteratorType begin, IteratorType end, const KeyType default_key)
{
	return PriorityQueueOnRange<IteratorType, KeyType, std::less<KeyType>>(begin, end, default_key);
}

//создать очередь с элементами в порядке убывания
template <class IteratorType, typename KeyType>
PriorityQueueOnRange<IteratorType, KeyType, std::less<KeyType>> MakeDecayProrityQueue(IteratorType begin, IteratorType end, const KeyType default_key)
{
	return PriorityQueueOnRange<IteratorType, KeyType, std::greater<KeyType>>(begin, end, default_key);
}

template <class IteratorType, typename KeyType, class Comparator>
PriorityQueueOnRange<IteratorType, KeyType, Comparator>::PriorityQueueOnRange(IteratorType begin, IteratorType end, const KeyType default_key)
{
	m_heap_size = end - begin;
	CV_Assert(m_heap_size > 0);
	m_heap.reserve(m_heap_size);
	m_positions.reserve(m_heap_size);
	int curr_ind = 0;

	for (IteratorType iter = begin; iter != end; ++iter)
	{
		m_heap.push_back(HeapElement(make_tuple(default_key, &*iter, curr_ind)));
		m_positions.push_back(curr_ind);
		curr_ind++;
	}
}

template <class IteratorType, typename KeyType, class Comparator>
void PriorityQueueOnRange<IteratorType, KeyType, Comparator>::ChangeKey(const size_t index, const KeyType new_key, const ThrowExcOrNot no_throw)
{
	if (index >= m_heap.size())
	{
		if (no_throw == eNotThrow)
		{
			return;
		}
		else
		{
			CV_Assert(index < m_heap.size());
		}
	}

	size_t curr_index = m_positions[index];

	if (curr_index >= m_heap_size)
	{
		if (no_throw == eNotThrow)
		{
			return;
		}
		else
		{
			CV_Assert(curr_index < m_heap_size);
		}
	}

	Comparator comparator;

	if(comparator(new_key, m_heap[curr_index].Key()))
	{
		m_heap[curr_index].Key() = new_key;
		size_t parent_index = curr_index / 2;

		while (comparator(m_heap[curr_index].Key(), m_heap[parent_index].Key()) &&
			curr_index != parent_index)
		{
			SwapElements(curr_index, parent_index);
			curr_index = parent_index;
			parent_index = curr_index / 2;
		}		
	}
}

template <class IteratorType, typename KeyType, class Comparator>
void PriorityQueueOnRange<IteratorType, KeyType, Comparator>::Pop()
{
	SwapElements(0, m_heap_size - 1);
	m_heap_size--;
	size_t curr_heap_ind = 0;
	size_t next_heap_ind = 0;
	//find the the location of the element
	do 
	{
		curr_heap_ind = next_heap_ind;
		next_heap_ind = Heapify(curr_heap_ind);
	} while (curr_heap_ind != next_heap_ind);
}

template <class IteratorType, typename KeyType, class Comparator>
size_t PriorityQueueOnRange<IteratorType, KeyType, Comparator>::Heapify(const size_t heap_index)
{
	const size_t left_child_ind = heap_index * 2;
	const size_t right_child_ind = heap_index * 2 + 1;
	size_t smallest_ind = m_heap_size;
	Comparator comparator;

	if (left_child_ind < m_heap_size && 
		comparator(m_heap[left_child_ind].Key(), m_heap[heap_index].Key()))
	{
		smallest_ind = left_child_ind;
	}
	else
	{
		smallest_ind  = heap_index;
	}

	if (right_child_ind < m_heap_size && 
		comparator(m_heap[right_child_ind].Key(), m_heap[smallest_ind].Key()))
	{
		smallest_ind = right_child_ind;
	}

	if (smallest_ind != heap_index)
	{
		SwapElements(smallest_ind, heap_index);
	}

	return smallest_ind;
}

template <class IteratorType, typename KeyType, class Comparator>
void PriorityQueueOnRange<IteratorType, KeyType, Comparator>::SwapElements(const size_t heap_index_1, const size_t heap_index_2)
{
	CV_Assert(heap_index_1 < m_heap_size);
	CV_Assert(heap_index_2 < m_heap_size);

	m_positions[m_heap[heap_index_1].Index()] = heap_index_2;
	m_positions[m_heap[heap_index_2].Index()] = heap_index_1;
	std::swap(m_heap[heap_index_1], m_heap[heap_index_2]);
}

}