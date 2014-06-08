#pragma once
#include <vector>
#include "opencv2/core/core.hpp"
#include "CommonUtilities\common_functions.h"

namespace molecule_descriptor
{
template <typename T>
size_t CalculateType(const T distance, const std::vector<T>& threshes)
{
	CV_Assert(!threshes.empty());
	int type = -1;

	for (auto iter = threshes.begin(); iter != threshes.end() - 1; ++iter)
	{
		if (distance >= *iter && distance < *(iter + 1))
		{
			type = static_cast<int>((iter - threshes.begin()) + 1);
			break;
		}
	}

	if (type == -1)
	{
		if (distance < threshes.front())
		{
			type = 0;
		}
		else if (distance >= threshes.back())
		{
			type = static_cast<int>(threshes.size());
		}
	}

	return static_cast<size_t>(type);
}

struct IntRange
{
	int min_val;
	int max_val;
	int step;

	int ElemsInRange() const
	{
		return (max_val - min_val) / step + 1;
	}
};

struct TypeWithMax
{
	size_t type;
	IntRange range;

	int Order() const
	{
		return (type - range.min_val) / range.step;
	}
};

inline size_t CalculateTypesComposition(const std::vector<TypeWithMax>& type_with_max)
{
	size_t type = 0;

	for (size_t curr_type = 0; curr_type < type_with_max.size(); ++curr_type)
	{
		if (curr_type > 0)
		{
			type *= type_with_max[curr_type].range.ElemsInRange();
		}

		type += type_with_max[curr_type].Order();
	}

	return type;
}

template <typename T>
void CalculateThresholdsQuantiles(const int intervals, std::vector<T>& arr, std::vector<T>& quantiles)
{
	std::sort(arr.begin(), arr.end());
	quantiles.resize(intervals);
	const double elems_in_block = static_cast<double>(arr.size()) / (intervals + 1);

	for (int block_ind = 1; block_ind <= intervals; block_ind++)
	{
		const size_t curr_ind = std::min(Round(block_ind * elems_in_block), static_cast<int>(arr.size() - 1));
		quantiles[block_ind - 1] = arr[curr_ind];
	}
}

}