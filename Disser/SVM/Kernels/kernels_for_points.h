#pragma once
#include "..\stdafx.h"
#include <utility>
#include <type_traits>
#include "shark\Core\IParameterizable.h"
#include "basic_kernels.h"
#include "../Common/singular_point.h"

namespace molecule_descriptor
{

template<typename ContainerType>
struct ScalarProductKernel : public shark::IParameterizable
{
	inline double operator()(const ContainerType& x, const ContainerType& y)
	{
		double res = 0.0;
		SHARK_ASSERT(x.size() == y.size());
		for (size_t ind = 0; ind < x.size(); ++ind)
		{
			res += x[ind] * y[ind];
		}
		return res;
	}

	virtual shark::RealVector parameterVector() const
	{
		return m_params_vect;
	}	 
	/// set the vector of hyper-parameters
	virtual void setParameterVector(RealVector const& new_parameters)
	{
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
private:
	shark::RealVector m_params_vect;
};
template<typename Kernel_1, typename Kernel_2>
struct KernelForPair : public shark::IParameterizable
{
	typedef typename Kernel_1::elem_type elem_type_1;
	typedef typename Kernel_2::elem_type elem_type_2;
	typedef std::pair<elem_type_1, elem_type_2> elem_type;
	KernelForPair(){

	}
	KernelForPair(const Kernel_1 kernel_1, const Kernel_2 kernel_2)
		: m_kernel_1(kernel_1),
		m_kernel_2(kernel_2) {
			CalcParamsVector();
	}

	inline double operator()(const elem_type& x, const elem_type& y)
	{
		return m_kernel_1(x.first, y.first) * m_kernel_2(x.second, y.second);
	}

	virtual shark::RealVector parameterVector() const
	{
		return m_params_vect;
	}	 
	/// set the vector of hyper-parameters
	virtual void setParameterVector(RealVector const& new_parameters)
	{
		SHARK_ASSERT(new_parameters.size() == m_params_vect.size());

		if (std::equal(m_params_vect.begin(), m_params_vect.end(), new_parameters.begin()))
		{
			return;
		}

		shark::init(new_parameters) >> shark::parameters(m_kernel_1), shark::parameters(m_kernel_2);
		CalcParamsVector();
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_kernel_1.numberOfParameters() + m_kernel_2.numberOfParameters());
		shark::init(m_params_vect) << shark::parameters(m_kernel_1), shark::parameters(m_kernel_1);
	}
private:
	Kernel_1 m_kernel_1;
	Kernel_2 m_kernel_2;
	shark::RealVector m_params_vect;
};

template<typename Kernel_1, typename Kernel_2, typename Kernel_3>
struct KernelForTuple : public shark::IParameterizable
{
	typedef typename Kernel_1::elem_type elem_type_1;
	typedef typename Kernel_2::elem_type elem_type_2;
	typedef typename Kernel_3::elem_type elem_type_3;
	typedef std::tuple<elem_type_1, elem_type_2, elem_type_3> elem_type;
	KernelForTuple(){

	}
	KernelForTuple(const Kernel_1 kernel_1, const Kernel_2 kernel_2, const Kernel_3 kernel_3)
		: m_kernel_1(kernel_1),
		m_kernel_2(kernel_2),
		m_kernel_3(kernel_3)
	{
			CalcParamsVector();
	}

	inline double operator()(const elem_type& x, const elem_type& y)
	{
		const double res_1 = m_kernel_1(std::get<0>(x), std::get<0>(y));
		double res_2 = 0;
		res_1 && (res_2 = m_kernel_2(std::get<1>(x), std::get<1>(y)));
		double res_3 = 0;
		res_2 && (res_3 = m_kernel_3(std::get<2>(x), std::get2>(y)));
		return res_1 * res_2 * res_3;
	}

	virtual shark::RealVector parameterVector() const
	{
		return m_params_vect;
	}	 
	/// set the vector of hyper-parameters
	virtual void setParameterVector(RealVector const& new_parameters)
	{
		SHARK_ASSERT(new_parameters.size() == m_params_vect.size());

		if (std::equal(m_params_vect.begin(), m_params_vect.end(), new_parameters.begin()))
		{
			return;
		}

		shark::init(new_parameters) >> shark::parameters(m_kernel_1), shark::parameters(m_kernel_2), shark::parameters(m_kernel_3);
		CalcParamsVector();
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_kernel_1.numberOfParameters() + m_kernel_2.numberOfParameters() + m_kernel_3.numberOfParameters());
		shark::init(m_params_vect) << shark::parameters(m_kernel_1), shark::parameters(m_kernel_2), shark::parameters(m_kernel_3);
	}
private:
	Kernel_1 m_kernel_1;
	Kernel_2 m_kernel_2;
	Kernel_3 m_kernel_3;
	shark::RealVector m_params_vect;
};

template<typename Kernel_1, typename Kernel_2, typename Kernel_3>
struct KernelForParametersSet : public shark::IParameterizable
{
	typedef typename Kernel_1::elem_type elem_type_1;
	typedef typename Kernel_2::elem_type elem_type_2;
	typedef typename Kernel_3::elem_type elem_type_3;
	typedef std::tuple<elem_type_1, elem_type_2, elem_type_3> tuple_type;
	typedef PropertiesSet elem_type;
	//static_assert(std::is_same<tuple_type, PropertiesSet::tuple_type>::value, "wrong_input");
	KernelForParametersSet(){
	}
	KernelForParametersSet(const Kernel_1 kernel_1, const Kernel_2 kernel_2, const Kernel_3 kernel_3)
		: m_kernel_1(kernel_1),
		m_kernel_2(kernel_2),
		m_kernel_3(kernel_3)
	{
		CalcParamsVector();
	}

	inline double operator()(const elem_type& x, const elem_type& y)
	{
		const double res_1 = m_kernel_1(x.SurfaceType(), y.SurfaceType());
		double res_2 = 0;

		if (res_1)
		{
			res_2 = m_kernel_2(x.Charge(), y.Charge());
		}

		double res_3 = 0;

		if (res_2)
		{
			res_3 = m_kernel_3(x.ElectricPotential(), y.ElectricPotential());
		}

		return res_1 * res_2 * res_3;
	}

	virtual shark::RealVector parameterVector() const
	{
		return m_params_vect;
	}	 
	/// set the vector of hyper-parameters
	virtual void setParameterVector(RealVector const& new_parameters)
	{
		SHARK_ASSERT(new_parameters.size() == m_params_vect.size());

		if (std::equal(m_params_vect.begin(), m_params_vect.end(), new_parameters.begin()))
		{
			return;
		}

		init(new_parameters) >> shark::blas::parameters(m_kernel_1), shark::blas::parameters(m_kernel_2), shark::blas::parameters(m_kernel_3);
		CalcParamsVector();
	}	 
	/// return the number of hyper-parameters
	virtual size_t numberOfParameters() const{ 
		return m_params_vect.size();
	}
protected:
	void CalcParamsVector()
	{
		m_params_vect.resize(m_kernel_1.numberOfParameters() + m_kernel_2.numberOfParameters() + m_kernel_3.numberOfParameters());
		init(m_params_vect) << shark::blas::parameters(m_kernel_1), shark::blas::parameters(m_kernel_2), shark::blas::parameters(m_kernel_3);
	}
private:
	Kernel_1 m_kernel_1;
	Kernel_2 m_kernel_2;
	Kernel_3 m_kernel_3;
	shark::RealVector m_params_vect;
};

struct LabelsKernel
{
	inline int operator()(std::pair<int, int>& point_1, std::pair<int, int>& point_2)
	{
		return DiracKernelFunc(point_1.first, point_2.first) && DiracKernelFunc(point_1.second, point_2.second);
	}
};

template<typename PropKernel>
struct SimplePropsKernel
{
	SimplePropsKernel(PropKernel kernel) : m_kernel(kernel) {}

	inline int operator()(std::pair<int, double>& point_1, std::pair<int, double>& point_2)
	{
		if (! DiracKernelFunc(point_1.first, point_2.first))
		{
			return 0;
		}

		return m_kernel(point_1.second, point_2.second);
	}

	PropKernel m_kernel;
};

}
