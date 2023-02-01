#include <stdint.h>
#include <cstdio>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include "InputOutput/params_reader.h"
#include "tools.h"


int main(int argc, char** argv)
{
	std::string input_file;
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-input_file", input_file, std::string("9.jpg"));
	cv::Mat_<uint8_t> input_im = cv::imread(input_file, cv::IMREAD_GRAYSCALE);	
	cv::Mat_<cv::Vec3b> input_im_color = cv::imread(input_file, cv::IMREAD_COLOR);
	//cv::resize(input_im, input_im, cv::Size(input_im.cols / 4, input_im.rows / 4));
	//cv::resize(input_im_color, input_im_color, cv::Size(input_im_color.cols / 4, input_im_color.rows / 4));
	cv::Mat input_im_color_copy = input_im_color.clone();
	
	cv::SiftFeatureDetector detector;
	std::vector<cv::KeyPoint> keypoints;
	detector.detect(input_im, keypoints);
	for (auto p = keypoints.begin(); p != keypoints.end(); ++p)
	{
		p->angle = -1;
	}
	// Add results to image and save.
	cv::Mat output;
	cv::drawKeypoints(input_im_color, keypoints, output, cv::Scalar::all(-1), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
	
	/*
	for (int y = 0; y < input_im_color.rows; ++y)
	{
		for (int x = 0; x < input_im_color.cols; ++x)
		{
			double norm = cv::norm(input_im_color(y, x));
			for (int i = 0; i < 3; ++i)
			{
				input_im_color(y, x)[i] = int(double(input_im_color(y, x)[i]) / norm * 255.0);
			}
		}
	}*/
	/*
	cv::Mat_<float> dst = cv::Mat_<float>::zeros(input_im.size());
	/// Parameters for Shi-Tomasi algorithm
	std::vector<cv::Point2f> corners;
	double qualityLevel = 0.3;
	double minDistance = 1;
	int blockSize = 7, gradientSize = 3;
	bool useHarrisDetector = true;
	double k = 0.04;

	/// Apply corner detection
	cv::goodFeaturesToTrack( input_im,
		corners,
		10000,
		qualityLevel,
		minDistance,
		cv::Mat(),
		blockSize,
		useHarrisDetector,
		k);


	/// Draw corners detected
	
	int radius = 4;
	cv::Mat input_im_color_copy = input_im_color.clone();
	for( size_t i = 0; i < corners.size(); i++ )
	{
		circle( input_im_color_copy, corners[i], 5,  cv::Scalar(0), 2, 8, 0 );
	}
	*/
	//input_im_color_copy.push_back(input_im_color);
	cv::imwrite(std::string("10.jpg"), output);
	//cv::imwrite(std::string("5.jpg"), input_im_color);
	/*
	cv::Mat_<double> input_im_double(input_im.size());
	input_im.convertTo(input_im_double, CV_64F);
	//cv::cvtColor(input_color, input_im, CV_BGR2GRAY);
	cv::Mat_<double> blurred(input_im.size());
	cv::Mat_<uint8_t> maximas(input_im.size());
	cv::Mat_<double> blurred_laplacian_signed(input_im.size());
	cv::Mat blurred_laplacian(input_im.size(), CV_8U);
	cv::Mat blurred_uint(input_im.size(), CV_8U);
	cv::Mat_<double> det(input_im.size()), gauss(input_im.size()), sobel(input_im.size());

	for (double sigma = 0.3; sigma < 100.0; sigma *= 1.5)
	{
		cv::GaussianBlur(input_im_double, blurred, cv::Size(), sigma);
		blurred /= sigma;
		im_test::FindLocalExtr(blurred, 1, maximas, true, true);
		std::stringstream str;
		str << sigma;
		cv::convertScaleAbs(blurred, blurred_uint, 0.05);
		cv::imwrite(std::string("1_") + str.str() + std::string(".bmp"), blurred_uint);
		cv::imwrite(std::string("m1_") + str.str() + std::string(".bmp"), maximas);
		const int rad = 2;//std::min(4, std::max(int(sigma / 3.0), 1));
		cv::Laplacian(blurred, blurred_laplacian_signed, CV_64F, 2 * rad + 1);
		cv::convertScaleAbs(blurred_laplacian_signed, blurred_laplacian, 0.05);
		im_test::FindLocalExtr(blurred_laplacian_signed, 1, maximas, true, true);

		cv::imwrite(std::string("2_") + str.str() + std::string(".bmp"), blurred_laplacian);
		cv::imwrite(std::string("m2_") + str.str() + std::string(".bmp"), maximas);

		im_test::HessianDeterminant(blurred, rad, det);

		//cv::imwrite(std::string("3_") + str.str() + std::string(".bmp"), det);
		im_test::FindLocalExtr(det, 1, maximas, true, true);
		cv::imwrite(std::string("m3_") + str.str() + std::string(".bmp"), maximas);	

		im_test::GaussianCurv(blurred, rad, gauss);
		im_test::FindLocalExtr(gauss, 1, maximas, true, false);
		cv::imwrite(std::string("m4_") + str.str() + std::string(".bmp"), maximas);	

		im_test::SobelNorm(blurred, rad, sobel);
		cv::convertScaleAbs(sobel, blurred_laplacian, 0.5);
		cv::imwrite(std::string("5_") + str.str() + std::string(".bmp"), blurred_laplacian);	
	}*/
}