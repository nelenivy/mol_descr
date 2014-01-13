#pragma once
#include <vector>
#include <algorithm>
#include <functional>

namespace molecule_descriptor
{

template <typename T, typename DistType>
void SegmentKMeans(const std::vector<T>& elems, const size_t clust_num, const std::function<DistType(T, T)> dist_calc, const T null_elem,
				   const DistType dist_thresh, const size_t max_iterations,
				   std::vector<size_t>& segmented_elems, std::vector<T>& centers)
{
	if (clust_num == 0)
	{
		return;
	}

	const size_t all_elems_num = elems.size();
	segmented_elems.resize(all_elems_num);
	std::fill(segmented_elems.begin(), segmented_elems.end(), 0);
	const size_t mean_elem_num_in_segm = all_elems_num / clust_num;
#pragma region init seeds
	centers.resize(clust_num + 1);//clust + 1 because clusters numbers go from 1
	std::vector<T> new_centers(clust_num + 1);//clust + 1 because clusters numbers go from 1

	for (size_t ind = 0, curr_clust = 1; ind < elems.size() && curr_clust <= clust_num; ind++)
	{
		if (ind % mean_elem_num_in_segm == 0)
		{
			centers[curr_clust] = elems[ind];
			curr_clust++;
		}
	}
#pragma endregion init seeds

#pragma region main cycle
	size_t changes_number = 1;
	size_t curr_iteration = 0;
	const DistType kMaxDist = 10000;

	std::vector<size_t> curr_clust_sizes(clust_num + 1, 0);

	while (changes_number > 0 && curr_iteration < max_iterations)
	{
		curr_iteration++;
		//assign elements to clusters
		for (int elem_ind = 0; elem_ind < elems.size(); ++elem_ind)
		{
			DistType min_dist = kMaxDist;

			for (int clust_ind = 1; clust_ind < centers.size(); ++clust_ind)
			{
				const DistType curr_dist = dist_calc(elems[elem_ind], centers[clust_ind]);

				if (curr_dist < min_dist)
				{
					min_dist = curr_dist;
					segmented_elems[elem_ind] = clust_ind;
				}
			}
		}

		//calculate new centers as mean of the cluster coordinates
		std::fill(new_centers.begin(), new_centers.end(), null_elem);
		std::fill(curr_clust_sizes.begin(), curr_clust_sizes.end(), 0);

		for (size_t ind = 0; ind < segmented_elems.size(); ++ind)
		{
			const size_t curr_clust_num = segmented_elems[ind];

				new_centers[curr_clust_num] += elems[ind];
				curr_clust_sizes[curr_clust_num]++;			
		}

		for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
		{
			new_centers[curr_clust] *= 1.0 / curr_clust_sizes[curr_clust];
		}
		//find the number of changes in centers coordinates
		changes_number = 0;

		for (size_t curr_clust = 1; curr_clust <= clust_num; ++curr_clust)
		{
			if (dist_calc(centers[curr_clust], new_centers[curr_clust]) > dist_thresh)
			{
				changes_number++;
			}
		}

		centers = new_centers;
	}
}

}