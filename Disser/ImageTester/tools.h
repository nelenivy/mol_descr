#pragma once
#include <string>
#include <stdint.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

namespace im_test
{

template <typename T>
void FindLocalExtr(const cv::Mat_<T>& input, const int rad, cv::Mat_<uint8_t>& output, const bool find_max, const bool find_min)
{
	output.create(input.size());

	for (int y = rad; y < input.rows - rad; ++y)
	{
		for (int x = rad; x < input.cols - rad; ++x)
		{
			bool max = true;
			bool min = true;

			for (int yw = y-rad; yw <= y + rad; ++yw)
			{
				for (int xw = x - rad; xw <= x + rad; ++xw)
				{
					if ((xw == x) && (yw == y))
					{
						continue;
					}
					if (input.at<T>(yw, xw) > input.at<T>(y, x))
					{
						max = false;
					}

					if (input.at<T>(yw, xw) < input.at<T>(y, x))
					{
						min = false;
					}
				}
			}

			if ((max && find_max) || (min &&find_min))
			{
				output.at<uint8_t>(y, x) = 255;
			}
			else
			{
				output.at<uint8_t>(y, x) = 0;
			}
		}
	}
}

template <typename T>
void HessianDeterminant(const cv::Mat_<T> input, const int rad, cv::Mat_<T> output)
{
	output.create(input.size());
	cv::Mat_<T> dxx(input.size()), dyy(input.size()), dxy(input.size()),
		det(input.size());

	cv::Sobel(input, dxx, cv::DataType<T>::type, 2, 0, 2 * rad + 1, 1.0);
	cv::Sobel(input, dyy, cv::DataType<T>::type, 0, 2, 2 * rad + 1, 1.0);
	cv::Sobel(input, dxy, cv::DataType<T>::type, 1, 1, 2 * rad + 1, 1.0);

	for (int y = 0; y < input.rows; ++y)
	{
		for (int x = 0; x < input.cols; ++x)
		{
			output.at<T>(y, x) = dxx.at<T>(y, x) * dyy.at<T>(y, x) - dxy.at<T>(y, x) * dxy.at<T>(y, x);
		}
	}
}

template <typename T>
void GaussianCurv(const cv::Mat_<T> input, const int rad, cv::Mat_<T> output)
{
	output.create(input.size());
	cv::Mat_<T> dx(input.size()), dy(input.size()), det(input.size());

	cv::Sobel(input, dx, cv::DataType<T>::type, 1, 0, 3, 1.0);
	cv::Sobel(input, dy, cv::DataType<T>::type, 0, 1, 3, 1.0);
	HessianDeterminant(input, rad, det);
	for (int y = 0; y < input.rows; ++y)
	{
		for (int x = 0; x < input.cols; ++x)
		{
			const T grad = T(1) + dx.at<T>(y, x) * dx.at<T>(y, x) + dy.at<T>(y, x) * dy.at<T>(y, x);
			output.at<T>(y, x) = std::abs(det.at<T>(y, x) / (grad * grad));
		}
	}
}

template <typename T>
void SobelNorm(const cv::Mat_<T> input, const int rad, cv::Mat_<T> output)
{
	output.create(input.size());
	cv::Mat_<T> dx(input.size()), dy(input.size());

	cv::Sobel(input, dx, cv::DataType<T>::type, 1, 0, 3, 1.0);
	cv::Sobel(input, dy, cv::DataType<T>::type, 0, 1, 3, 1.0);
	for (int y = 0; y < input.rows; ++y)
	{
		for (int x = 0; x < input.cols; ++x)
		{
			output.at<T>(y, x) = sqrt(dx.at<T>(y, x) * dx.at<T>(y, x) + dy.at<T>(y, x) * dy.at<T>(y, x));
		}
	}
}

}
