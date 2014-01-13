#pragma once
#include <math.h>
#include <array>
#include <shark/Core/IParameterizable.h>
//#include <shark/Models/Kernels/GaussianRbfKernel.h>
#include "opencv2/core/core.hpp"

namespace molecule_descriptor
{

using  shark::RealVector;

class GaussianKernelOneDim /*: public shark::IParameterizable*/
{
public:
	typedef double elem_type;

	GaussianKernelOneDim(const double sigma)
		: m_param_vect(1)
	{
		SetSigma(sigma);
	}

	GaussianKernelOneDim()
		: m_param_vect(1)
	{
		SetSigma(1.0);
	}

	inline double operator()(double x, double y)
	{
		const double abs_arg = abs(x - y);

		if (abs_arg > m_threshold)
		{
			return 0.0;
		}
		else
		{
			return std::exp( - abs_arg * abs_arg * m_precalculated_multiplyer);
		}
	}

	void SetSigma(const double sigma)
	{
		m_sigma = sigma;
		m_param_vect[0] = m_sigma;
		m_precalculated_multiplyer = 1.0 / (2.0 * m_sigma * m_sigma);
		m_threshold = 4.0 * sigma;
	}

	/*virtual*/ shark::RealVector 	parameterVector () const {
		 return m_param_vect;
	}

	/*virtual*/ void setParameterVector (shark::RealVector const& new_parameters) {
		CV_Assert(new_parameters.size() == 1);
		SetSigma(new_parameters[0]);
	}

	/*virtual*/ std::size_t 	numberOfParameters () const	{
		return m_param_vect.size();
	}

	double Sigma() { 
		return m_sigma; }
private:
	double m_sigma;
	double m_threshold;
	double m_precalculated_multiplyer;
	shark::RealVector m_param_vect;
};

class MultiplicativeGaussianKernel
{
public:
	MultiplicativeGaussianKernel(const double sigma) : m_sigma(sigma) 
	{
		m_precalculated_multiplyer = 1.0 / (2.0 * m_sigma * m_sigma);
	}

	template <int kArrSize>
	inline double operator()(std::array<double, kArrSize>& x, std::array<double, kArrSize>& y)
	{
		double power = 0;

		for (int ind = 0; ind < kArrSize; power += pow(x[ind] - y[ind], 2), ind++)
		{}

		power *= -m_precalculated_multiplyer;

		return std::exp(power);
	}

	double Sigma() { 
		return m_sigma; }
private:
	double m_sigma;
	double m_precalculated_multiplyer;
};

class TriangleKernel
{
public:
	typedef double elem_type;

	TriangleKernel(const double C)
		: m_C(C) {}
	inline double operator()(double x, double y)
	{
		return std::max(0.0, (m_C - abs(x - y)) / m_C);
	}
private:
	double m_C;
};

class TriangleKernelThresholded  /*	: public shark::IParameterizable*/
{
public:
	typedef double elem_type;
	TriangleKernelThresholded(){

	}
	TriangleKernelThresholded(const double C, const double x_min, const double x_max)
		: m_C(C), 
		m_x_min(x_min), 
		m_x_max(x_max)
	{
	}
	/*virtual*/ shark::RealVector 	parameterVector () const {
		return shark::RealVector();
	}

	/*virtual*/ void setParameterVector (shark::RealVector const& new_parameters) {
	}

	/*virtual*/ std::size_t 	numberOfParameters () const	{
		return 0;
	}
	inline double operator()(double x, double y)
	{
		const double x_thresh = std::max(std::min(x, m_x_max), m_x_min);
		const double y_thresh = std::max(std::min(y, m_x_max), m_x_min);
		return std::max(0.0, (m_C - abs(x_thresh - y_thresh)) / m_C);
	}
private:
	double m_C;
	double m_x_min;
	double m_x_max;
};

inline int DiracKernelFunc(int x, int y)
{
	return x == y ? 1 : 0;
}

template<typename T>
struct DiracKernel
{
	typedef T elem_type;

	inline int operator()(const T x, const T y)
	{
		return x == y ? 1 : 0;
	}

	/*virtual*/ shark::RealVector 	parameterVector () const {
		return shark::RealVector();
	}

	/*virtual*/ void setParameterVector (shark::RealVector const& new_parameters) { }

	/*virtual*/ std::size_t 	numberOfParameters () const	{
		return 0;
	}
};

template<typename T>
struct AlwaysOneKernel /*: public shark::IParameterizable*/
{
	typedef T elem_type;

	inline int operator()(const T x, const T y)
	{
		return 1;
	}

	/*virtual*/ shark::RealVector 	parameterVector () const {
		return shark::RealVector();//fictive parameter
	}

	/*virtual*/ void setParameterVector (shark::RealVector const& new_parameters) { }

	/*virtual*/ std::size_t 	numberOfParameters () const	{
		return 0;
	}
};

}