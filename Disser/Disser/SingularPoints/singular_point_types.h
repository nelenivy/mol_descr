#pragma once
#include "CommonUtilities/attributes_container.h"

namespace molecule_descriptor
{
//mark curvature type
const unsigned char kConvexType = 1;
const unsigned char kConcaveType = 2;
const unsigned char kSaddleType = 3;

ENUM_TYPE(SurfaceTypeProp, unsigned char, 1, 3);
ENUM_TYPE(ElectricSign, char, -1, 1);
SIZE_T_PROP(SingularPointType);
DOUBLE_PROP(ElectricPotential);

//////////////////////////////////////////////////////////////////////////
template <typename T>
size_t GetRange()
{
	return T::kMaxVal - T::kMinVal + 1;
}
//function for property calculation
//template <class ... Props>
//struct CalculateSingularPointType;

//template <class Prop1>
//struct CalculateSingularPointType/*<Prop1>*/
//{
//	size_t operator()(AttributesContainer& attr)
//	{	
//		return static_cast<size_t>(attr.Get<Prop1>() - Prop1::kMinVal);
//	}
//};

template <class Prop1>
size_t CalculateSingularPointType(AttributesContainer& attr)
{	
	return static_cast<size_t>(attr.Get<Prop1>() - Prop1::kMinVal);
}

//template <class Prop1, class ... Props>
//struct CalculateSingularPointType<Prop1, Props...>
//{
//	size_t operator()(AttributesContainer& attr)
//	{
//		return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * 
//		CalculateSingularPointType<Props...>()(attr) +  
//		CalculateSingularPointType<Prop1>()(attr); 
//	}
//};

//Uncomment this if you don't have c++11
template <class Prop1, class Prop2>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5, class Prop6>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5, Prop6>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5, class Prop6, class Prop7>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5, Prop6, Prop7>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5, class Prop6, class Prop7, class Prop8>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5, Prop6, Prop7, Prop8>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5, class Prop6, class Prop7, class Prop8, class Prop9>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5, Prop6, Prop7, Prop8, Prop9>(attr) +  CalculateSingularPointType<Prop1>(attr); }

template <class Prop1, class Prop2, class Prop3, class Prop4, class Prop5, class Prop6, class Prop7, class Prop8, class Prop9, class Prop10>
size_t CalculateSingularPointType(AttributesContainer& attr){	
	return static_cast<size_t>(Prop1::kMaxVal - Prop1::kMinVal + 1) * CalculateSingularPointType<Prop2, Prop3, Prop4, Prop5, Prop6, Prop7, Prop8, Prop9, Prop10>(attr) +  CalculateSingularPointType<Prop1>(attr); }
//////////////////////////////////////////////////////////////////////////
}