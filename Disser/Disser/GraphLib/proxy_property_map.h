#pragma once 

#include <type_traits>
#include "boost/property_map/property_map.hpp"

namespace molecule_descriptor
{

template <typename PropMap, typename TransformFunc, typename DataKeeper>
class ProxyPropMapBase
	: boost::put_get_helper< 
	typename std::result_of<TransformFunc(typename boost::property_traits<PropMap>::reference)>::type ,
	ProxyPropMapBase<PropMap, TransformFunc, DataKeeper>>
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

	OutRef& operator[](const key_type key)
	{
		DataKeeper& data_keeper = static_cast<DataKeeper&>(*this);
		return data_keeper.m_func(data_keeper.m_prop_map[key]);
	}
	const OutRef& operator[](const key_type key) const
	{
		const DataKeeper& data_keeper = static_cast<const DataKeeper&>(*this);
		return const_cast<const OutRef&>(data_keeper.m_func(data_keeper.m_prop_map[key]));
	}
private:
};


template <typename PropMap, typename TransformFunc>
class ProxyPropMapVal : public ProxyPropMapBase<PropMap, TransformFunc, ProxyPropMapVal<PropMap, TransformFunc>>
{
public:
	ProxyPropMapVal() {}
	ProxyPropMapVal(PropMap& prop_map, const TransformFunc func) : m_prop_map(prop_map), m_func(func) {}
	friend class ProxyPropMapBase<PropMap, TransformFunc, ProxyPropMapVal<PropMap, TransformFunc>>;
private:
	PropMap m_prop_map;
	TransformFunc m_func;
};
template <typename PropMap, typename TransformFunc>
class ProxyPropMapRef : public ProxyPropMapBase<PropMap, TransformFunc, ProxyPropMapRef<PropMap, TransformFunc>>
{
public:
	ProxyPropMapRef(PropMap& prop_map, const TransformFunc func) : m_prop_map(prop_map), m_func(func) {}
	friend class ProxyPropMapBase<PropMap, TransformFunc, ProxyPropMapRef<PropMap, TransformFunc>>;
private:
	PropMap& m_prop_map;
	TransformFunc m_func;
};

template <class PropMap, class TransformFunc>
ProxyPropMapRef<PropMap, TransformFunc> GetProxyPropMapRef(PropMap& prop_map, const TransformFunc& transform_func)
{
	return ProxyPropMapRef<PropMap, TransformFunc>(prop_map, transform_func);
}

template <class PropMap, class TransformFunc>
ProxyPropMapVal<PropMap, TransformFunc> GetProxyPropMapVal(PropMap& prop_map, const TransformFunc& transform_func)
{
	return ProxyPropMapVal<PropMap, TransformFunc>(prop_map, transform_func);
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
ProxyPropMapRef<PropMap, IdenticalTransformFunc> GetPropMapReference(PropMap& prop_map)
{
	return ProxyPropMapRef<PropMap, IdenticalTransformFunc>(prop_map, IdenticalTransformFunc());
}
}