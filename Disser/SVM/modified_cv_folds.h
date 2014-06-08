#pragma once
#include "stdafx.h"
#include <algorithm>
#include <vector>

#include <shark/Data/Dataset.h>
#include <shark/Rng/DiscreteUniform.h>
//#include <shark/SharkDefs.h>
#include <shark/Data/DataView.h>
#include <shark/Data/CVDatasetTools.h>

#include "element_with_index_and_id.h"

namespace shark {

namespace detail {

template<class ElemType, typename DataIDType, class L>
void FindElementDublicates(const LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set, 
						   const std::vector< std::vector<std::size_t> >& members,
						   std::vector<std::vector<std::vector<size_t> > >& elems)
{
	const size_t set_size = set.numberOfElements();
	std::vector<bool> visited(set_size, false);
	elems = std::vector<std::vector<std::vector<size_t> > >(members.size());
	const auto& set_inputs = set.inputs();

	for (size_t class_num = 0; class_num < members.size(); ++class_num)
	{
		const std::vector<std::size_t>& curr_members = members[class_num];
		for (size_t elem_num = 0; elem_num < curr_members.size(); ++elem_num)
		{
			const size_t curr_elem_num = curr_members[elem_num];
			if (visited[curr_elem_num])
			{
				continue;
			}
			
			const auto& curr_elem = set_inputs.element(curr_elem_num);
			const int curr_ind = curr_elem.Index();
			elems[class_num].push_back(std::vector<size_t>());

			for (size_t set_ind = 0; set_ind < set_size; ++set_ind)
			{
				if (set_inputs.element(set_ind).Index() == curr_ind)
				{
					visited[set_ind] = true;
					elems[class_num].back().push_back(set_ind);
				}
			}
		}
	}

	//sort elems according to their size
	struct Comparator
	{
		bool operator()(const std::vector<size_t>& elem_1, const std::vector<size_t>& elem_2)
		{
			return elem_1.size() < elem_2.size();
		}
	};

	for (size_t class_num = 0; class_num < members.size(); ++class_num)
	{
		std::sort(elems[class_num].begin(), elems[class_num].end(), Comparator());
	}
}
///\brief Version of createCVSameSizeBalanced which works regardless of the label type
///
/// This function for every class requires one vector to store the indices of.. todo
// todo: help: correct me:
//This functions needs for every class a vector which stores the indices of the vectors  of "set" which are part of this class.
	template<class ElemType, typename DataIDType, class L>
	std::vector<std::size_t> createCVSameSizeBalancedIndexedElemsHelp(
	LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set,
	std::size_t numberOfPartitions,
	std::vector< std::vector<std::size_t> > members,
	std::size_t batchSize
	) {
		typedef molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType> I;

		std::size_t numInputs = set.numberOfElements();
		std::size_t numClasses = members.size();

		//shuffle elements in members
		DiscreteUniform< Rng::rng_type > uni(shark::Rng::GetCurrThreadRng()) ;
		for (std::size_t c = 0; c != numClasses; c++) {
			std::random_shuffle(members[c].begin(), members[c].end(), uni);
		}

		//calculate number of elements per validation subset in the new to construct container
		std::size_t nn = numInputs / numberOfPartitions;
		std::size_t leftOver = numInputs % nn;
		std::vector<std::size_t> validationSize(numberOfPartitions,nn);
		for (std::size_t partition = 0; partition != leftOver; partition++) {
			validationSize[partition]++;
		}

		//calculate the size of the batches for every validation part
		std::vector<std::size_t> partitionStart;
		std::vector<std::size_t> batchSizes;
		std::size_t numBatches = batchPartitioning(validationSize,partitionStart,batchSizes,batchSize);


		LabeledData<I,L> newSet(numBatches);//set of empty batches
		DataView<LabeledData<I,L> > setView(set);//fast access to single elements of the original set
		std::vector<std::size_t> validationSetStart = partitionStart;//current index for the batch of every fold
		//partition classes into the validation subsets of newSet
		std::vector<std::vector<std::size_t> > batch_elements(numberOfPartitions);
		std::vector<std::vector<std::vector<size_t>>> eq_members;
		FindElementDublicates(set, members, eq_members);

		for (std::size_t c = 0, fold = 0; c != numClasses; c++) {
			for (std::size_t i = 0; i != eq_members[c].size(); i++) {
				const std::vector<std::size_t> oldPos = eq_members[c][i];
				std::size_t batchNumber = validationSetStart[fold];

				for (size_t pos_ind = 0; pos_ind < oldPos.size(); ++pos_ind)
				{
					batch_elements[fold].push_back(oldPos[pos_ind]);
				}

				//if all elements for the current batch are found, create it
				const size_t next_fold = (fold+1) % numberOfPartitions;
				const size_t next_fold_first = partitionStart[next_fold];
				const size_t next_batch_num = batchNumber + 1;

				if (next_batch_num != next_fold_first && 
					next_batch_num != batchSizes.size() && 
					batch_elements[fold].size() >= batchSizes[batchNumber]) 
				{
					newSet.batch(validationSetStart[fold]) = subBatch(setView,batch_elements[fold]);
					batch_elements[fold].clear();
					++validationSetStart[fold];
				}

				fold = (fold+1) % numberOfPartitions;
			}
		}
		//push non-pushed batches
		for (size_t curr_fold = 0; curr_fold < validationSetStart.size(); ++curr_fold)
		{
			//if all elements for the current batch are found, create it
			if (!batch_elements[curr_fold].empty()) {
				newSet.batch(validationSetStart[curr_fold]) = subBatch(setView,batch_elements[curr_fold]);
				batch_elements[curr_fold].clear();
				++validationSetStart[curr_fold];
			}
		}
		//count batches num
		size_t new_batches_num = 0;
		for (size_t curr_batch = 0; curr_batch < newSet.numberOfBatches(); ++curr_batch)
		{
			if (newSet.batch(curr_batch).size() > 0)
			{
				++new_batches_num;
			}
		}
		LabeledData<I,L> new_new_set(new_batches_num);//set of empty batches

		for (size_t curr_batch = 0, curr_new_batch = 0, curr_fold = 0; curr_batch < newSet.numberOfBatches(); ++curr_batch)
		{
			if (newSet.batch(curr_batch).size() > 0)
			{
				new_new_set.batch(curr_new_batch) = newSet.batch(curr_batch);
				++curr_new_batch;
			}
			else
			{
				for (size_t next_fold = 0; next_fold < partitionStart.size(); ++next_fold)
				{
					if (partitionStart[next_fold] > curr_new_batch)
					{
						partitionStart[next_fold]--;
					}
				}
			}
		}

		for (size_t next_fold = 0; next_fold < partitionStart.size(); )
		{
			if ( next_fold < partitionStart.size() - 1 && partitionStart[next_fold] == partitionStart[next_fold + 1]) 
				
			{
				partitionStart.erase(partitionStart.begin() + next_fold + 1);
			}
			else if ((next_fold == partitionStart.size() - 1) &&
				 (partitionStart[next_fold] >= new_new_set.numberOfBatches()) )
			{
				partitionStart.erase(partitionStart.begin() + next_fold);
			}
			else
			{
				++next_fold;
			}
		}
		//swap old and new set
		swap(set, new_new_set);

		//create folds
		return partitionStart;
}

template<class ElemType, typename DataIDType, class L>
CVFolds<LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> > createCVSameSizeBalancedIndexedElems
	(
	LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set,
	std::size_t numberOfPartitions,
	std::vector< std::vector<std::size_t> > members,
	std::size_t batchSize
	)
{
	auto part = createCVSameSizeBalancedIndexedElemsHelp(set, numberOfPartitions, members, batchSize);
	return CVFolds<LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> >(set, part);
}
template<class ElemType, typename DataIDType, class L>
CVFolds<LabeledData<ElemType,L> > createCVSameSizeBalancedNoIndexedElems
	(
	LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set,
	std::size_t numberOfPartitions,
	std::vector< std::vector<std::size_t> > members,
	std::size_t batchSize
	)
{
	auto part = createCVSameSizeBalancedIndexedElemsHelp(set, numberOfPartitions, members, batchSize);
	auto new_set = CopyFromIndexedToNonIndexed(set);
	return CVFolds<LabeledData<ElemType,L> >(new_set, part);
}

}//nam

//! \brief Create a partition for cross validation
//!
//! Every subset contains (approximately) the same
//! number of elements. For every partition, all
//! but one subset form the training set, while the
//! remaining one is used for validation.
//!
//! \param numberOfPartitions  number of partitions to create
//! \param set the input data from which to draw the partitions
template<class ElemType, typename DataIDType, class L>
CVFolds<LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> > createCVSameSizeBalancedIndexedElems(
	LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set, 
	size_t numberOfPartitions, 
	std::size_t batchSize=Data<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>>::DefaultBatchSize) {
		typedef molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType> I;
	DataView<LabeledData<I,unsigned int> > setView(set);
	std::size_t numInputs = setView.size();
	std::size_t numClasses = numberOfClasses(set);

	//find members of each class
	std::vector< std::vector<std::size_t> > members(numClasses);
	for (std::size_t i = 0; i != numInputs; i++) {
		members[setView[i].label].push_back(i);
	}
	return detail::createCVSameSizeBalancedIndexedElems(set,numberOfPartitions,members,batchSize);

}

template<class ElemType, typename DataIDType, class L>
CVFolds<LabeledData<ElemType, L> > createCVSameSizeBalancedNoIndexedElems(
	LabeledData<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>,L> &set, 
	size_t numberOfPartitions, 
	std::size_t batchSize=Data<molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType>>::DefaultBatchSize) {
		typedef molecule_descriptor::ElemWithIndexAndID<ElemType, DataIDType> I;
		DataView<LabeledData<I,unsigned int> > setView(set);
		std::size_t numInputs = setView.size();
		std::size_t numClasses = numberOfClasses(set);

		//find members of each class
		std::vector< std::vector<std::size_t> > members(numClasses);
		for (std::size_t i = 0; i != numInputs; i++) {
			members[setView[i].label].push_back(i);
		}
		return detail::createCVSameSizeBalancedNoIndexedElems(set,numberOfPartitions,members,batchSize);

}


}