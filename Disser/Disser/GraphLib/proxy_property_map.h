#pragma once 

#include <type_traits>
#include "boost/property_map/property_map.hpp"

namespace molecule_descriptor
{

template <typename PropMap, typename TransformFunc>
class ProxyPropMap
	: boost::put_get_helper< 
	typename std::result_of<TransformFunc(typename boost::property_traits<PropMap>::reference)>::type ,
	ProxyPropMap<PropMap, TransformFunc>>
{
public:
	typedef TransformFunc TransformFuncType;
	typedef typename 
		std::result_of<TransformFunc(typename boost::property_traits<PropMap>::reference)>::type OutRef;
	static_assert(std::is_reference<OutRef&>::value, "The result of TransformFunc must be reference");
	typedef typename boost::property_traits<PropMap>::key_type		key_type;
	typedef typename std::remove_reference<OutRef&>::type			value_type;
	typedef OutRef&													reference;
	typedef typename boost::property_traits<PropMap>::category		category;

	ProxyPropMap(PropMap& prop_map, const TransformFunc func) : m_prop_map(prop_map), m_func(func) {}

	OutRef& operator[](const key_type key)
	{
		return m_func(m_prop_map[key]);
	}
	const OutRef& operator[](const key_type key) const
	{
		return const_cast<const OutRef&>(m_func(m_prop_map[key]));
	}
private:
	PropMap& m_prop_map;
	TransformFunc m_func;
};

template <class PropMap, class TransformFunc>
ProxyPropMap<PropMap, TransformFunc> GetProxyPropMap(PropMap& prop_map, const TransformFunc& transform_func)
{
	return ProxyPropMap<PropMap, TransformFunc>(prop_map, transform_func);
}

struct IdenticalTransformFunc
{
	template <typename T>
	T& operator ()(T& elem) const
	{
		return elem;
	}
	template <typename T>
	const T& operator ()(const T& elem) const
	{
		return elem;
	}
};

template <class PropMap>
ProxyPropMap<PropMap, IdenticalTransformFunc> GetPropMapReference(PropMap& prop_map)
{
	return ProxyPropMap<PropMap, IdenticalTransformFunc>(prop_map, IdenticalTransformFunc());
}
}