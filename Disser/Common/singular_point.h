#pragma once

#include <utility>
#include <tuple>
#include <array>
#include <stdint.h>
#include "opencv2/core/core.hpp"

#include "interval_read_write.h"

namespace molecule_descriptor
{

template <typename PropType>
class SingularPoint
{
public:
	typedef PropType property_type;
	SingularPoint()
	{
	}
	SingularPoint(const cv::Point3d& point, const PropType& property)
		: m_singular_point(point, property)
	{
	}

	inline cv::Point3d& Coord()
	{
		return m_singular_point.first;
	}

	inline const cv::Point3d& Coord() const
	{
		return m_singular_point.first;
	}

	inline PropType& Property()
	{
		return m_singular_point.second;
	}

	inline const PropType& Property() const
	{
		return m_singular_point.second;
	}

	inline std::pair<cv::Point3d, PropType>& GetAsPair()//for convinience
	{
		return m_singular_point;
	}
	inline const std::pair<cv::Point3d, PropType>& GetAsPair() const//for convinience
	{
		return m_singular_point;
	}
private:
	std::pair<cv::Point3d, PropType> m_singular_point;
};


class PropertiesSet
{
public:
	typedef std::tuple<size_t, int, float> tuple_type;
	inline size_t& SurfaceType() {
		return std::get<0>(m_properties);
	}
	inline const size_t& SurfaceType() const {
		return std::get<0>(m_properties);
	}

	inline int& Charge() {
		return std::get<1>(m_properties);
	}
	inline const int& Charge() const {
		return std::get<1>(m_properties);
	}

	inline float& ElectricPotential() {
		return std::get<2>(m_properties);
	}
	inline const float& ElectricPotential() const {
		return std::get<2>(m_properties);
	}

	inline std::tuple<size_t, int, float>& GetAsTuple()//for convinience
	{
		return m_properties;
	}
	inline const std::tuple<size_t, int, float>& GetAsTuple() const//for convinience
	{
		return m_properties;
	}
private: 
	std::tuple<size_t, int, float> m_properties;
};

template <typename PropType>
struct WriteElemToFile<SingularPoint<PropType>>;
template <typename PropType>
struct ReadElemNonChecked<SingularPoint<PropType>>;
template <>
struct WriteElemToFile<PropertiesSet>;
template <>
struct ReadElemNonChecked<PropertiesSet>;
//specialization of output
template <typename PropType>
struct WriteElemToFile<SingularPoint<PropType>>
{
	inline void operator()(std::ofstream& file_out, const SingularPoint<PropType>& elem)
	{
		WriteElemToFile<std::pair<cv::Point3d, PropType>>()(file_out, elem.GetAsPair());
	}
};
//specialization of input
template <typename PropType>
struct ReadElemNonChecked<SingularPoint<PropType>>
{
	static inline void Do(std::ifstream& file_in, SingularPoint<PropType>& new_elem)
	{
		ReadElemNonChecked<std::pair<cv::Point3d, PropType>>::Do(file_in, new_elem.GetAsPair());
	}
};
//specialization of output
template <>
struct WriteElemToFile<PropertiesSet>
{
	inline void operator()(std::ofstream& file_out, const PropertiesSet& elem)
	{
		WriteElemToFile<std::tuple<size_t, int, float>>()(file_out, elem.GetAsTuple());
	}
};
//specialization of input
template <>
struct ReadElemNonChecked<PropertiesSet>
{
	static inline void Do(std::ifstream& file_in, PropertiesSet& new_elem)
	{
		ReadElemNonChecked<std::tuple<size_t, int, float>>::Do(file_in, new_elem.GetAsTuple());
	}
};
//singular point with histogram
template <int kHistSize>
class HistogramSingularPoint : public SingularPoint<std::array<uint8_t, kHistSize>>
{

};

//specialization of output
template <int kHistSize>
struct WriteElemToFile<HistogramSingularPoint<kHistSize>>
{
	inline void operator()(std::ofstream& file_out, const HistogramSingularPoint<kHistSize>& elem)
	{
		WriteElemToFile<SingularPoint<std::array<uint8_t, kHistSize>>>()(file_out, elem);
	}
};
//specialization of input
template <int kHistSize>
struct ReadElemNonChecked<HistogramSingularPoint<kHistSize>>
{
	static inline void Do(std::ifstream& file_in, HistogramSingularPoint<kHistSize>& new_elem)
	{
		ReadElemNonChecked<SingularPoint<std::array<uint8_t, kHistSize>>>::Do(file_in, new_elem);
	}
};
typedef SingularPoint<size_t> MarkedSingularPoint;
typedef SingularPoint<PropertiesSet> NonMarkedSingularPoint;

}