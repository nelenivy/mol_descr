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

//mark curvature type
enum CurvatureTypes {kConvexType = 1, kConcaveType = 2, kSaddleType = 3};

class PropertiesSet
{
public:
	typedef std::tuple<size_t, int, float, int, float> tuple_type;
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
	inline int& LennardJonesType() {
		return std::get<3>(m_properties);
	}
	inline const int& LennardJonesType() const {
		return std::get<3>(m_properties);
	}
	inline float& LennardJones() {
		return std::get<4>(m_properties);
	}
	inline const float& LennardJones() const {
		return std::get<4>(m_properties);
	}

	inline tuple_type& GetAsTuple()//for convinience
	{
		return m_properties;
	}
	inline const tuple_type& GetAsTuple() const//for convinience
	{
		return m_properties;
	}
private: 
	tuple_type m_properties;
};

inline bool operator ==(const PropertiesSet& p1, const PropertiesSet& p2)
{
	return (p1.Charge() == p2.Charge()) &&
		(p1.ElectricPotential() == p2.ElectricPotential())
		&& (p1.SurfaceType() == p2.SurfaceType())
		&& (p1.LennardJonesType() == p2.LennardJonesType())
		&& (p1.LennardJones() == p2.LennardJones());
}
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
		WriteElemToFile<PropertiesSet::tuple_type>()(file_out, elem.GetAsTuple());
	}
};
//specialization of input
template <>
struct ReadElemNonChecked<PropertiesSet>
{
	static inline void Do(std::ifstream& file_in, PropertiesSet& new_elem)
	{
		ReadElemNonChecked<PropertiesSet::tuple_type>::Do(file_in, new_elem.GetAsTuple());
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
//////////////////////////////////////////////////////////////////////////
template <typename PropT>
struct SingularPointsPair
{
	PropT elem1;
	PropT elem2;
	double dist;
};

template <typename T>
bool operator ==(const SingularPointsPair<T>& t1, const SingularPointsPair<T>& t2)
{
	return ((t1.elem1 == t2.elem1 && t1.elem2 == t2.elem2) || 
		(t1.elem1 == t2.elem2 && t1.elem2 == t2.elem1)) && (t1.dist == t2.dist);
}
template <typename T>
bool operator !=(const SingularPointsPair<T>& t1, const SingularPointsPair<T>& t2)
{
	return !(t1 == t2);
}
//////////////////////////////////////////////////////////////////////////

typedef SingularPoint<size_t> MarkedSingularPoint;
typedef SingularPoint<PropertiesSet> NonMarkedSingularPoint;

}