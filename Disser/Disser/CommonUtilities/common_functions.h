#pragma once
#include "opencv2/core/core.hpp"

namespace molecule_descriptor
{

template <typename T>
T Sqr(const T num)
{
	return num * num;
}

template <typename T>
int Round(T num)
{
	return num > static_cast<T>(0) ? static_cast<int>(num + static_cast<T>(0.5)) : static_cast<int>(num - static_cast<T>(0.5));
}

template <typename T>
char Sign(T num)
{
	return num > static_cast<T>(0) ? 1 : (num < static_cast<T>(0) ? -1 : 0);
}

struct ReleaseDeleter
{
	template <typename T> 
	void operator()(T* ptr) {
		ptr->Release();
	}
};

template<typename T>
void Point_ToMat_(const cv::Point3_<T>& point, cv::Mat_<T>& mat)
{
	T point_arr[] = {point.x, point.y, point.z};
	cv::Mat_<T> temp_mat(1, 3, point_arr);
	temp_mat.copyTo(mat);
}

template<typename T>
void Point_ToMat_(const cv::Point_<T>& point, cv::Mat_<T>& mat)
{
	T point_arr[] = {point.x, point.y, point.z};
	cv::Mat_<T> temp_mat(1, 2, point_arr);
	temp_mat.copyTo(mat);
}

template<typename T>
void Point_ToMat_Transposed(const cv::Point3_<T>& point, cv::Mat_<T>& mat)
{
	T point_arr[] = {point.x, point.y, point.z};
	cv::Mat_<T> temp_mat(3, 1, point_arr);
	temp_mat.copyTo(mat);
}

template<typename T>
void CopyPoint_ToCol(const cv::Point3_<T>& point, const int col_num, cv::Mat_<T>& mat)
{
	T point_arr[] = {point.x, point.y, point.z};
	for (int i = 0; i < 3; ++i)
	{
		mat(i, col_num) = point_arr[i];
	}
}

template<typename T>
void Point_ToMat_Transposed(const cv::Point_<T>& point, cv::Mat_<T>& mat)
{
	T point_arr[] = {point.x, point.y, point.z};
	cv::Mat_<T> temp_mat(3, 1, point_arr);
	temp_mat.copyTo(mat);
}

template<typename T, int N>
void Mat_ToVec(const cv::Mat_<T>& mat, cv::Vec<T, N>& point)
{
	std::copy(mat[0], mat[0] + N, &(point[0]));
}

template<typename T, int N>
void Vec_ToMat_(const cv::Vec<T, N>& point, cv::Mat_<T>& mat)
{
	mat.create(1, N);
	std::copy(&(point[0]), &(point[0]) + N, mat[0]);
}
template<typename T, int N>
void Vec_ToMat_Transposed(const cv::Vec<T, N>& point, cv::Mat_<T>& mat)
{
	mat.create(N, 1);
	std::copy(&(point[0]), &(point[0]) + N, mat[0]);
}

template<typename IteratorType>
size_t CalculateRangeType(IteratorType types_begin, IteratorType types_end, 
						  IteratorType max_vals_begin, IteratorType max_vals_end)
{
	CV_Assert(types_end - types_begin == max_vals_end - max_vals_begin);
	size_t prev_step_val = *types_begin;

	for (++types_begin, ++max_vals_begin; types_begin != types_end && max_vals_begin != max_vals_end; ++types_begin, ++max_vals_begin)
	{
		prev_step_val = prev_step_val * (*max_vals_begin) + *types_begin;
	}

	return prev_step_val;
}


}