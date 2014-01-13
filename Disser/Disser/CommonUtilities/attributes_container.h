#pragma once
#include <stdint.h>
#include <vector>
#include <utility>
#include <typeinfo>
#include <typeindex>
#include <type_traits>
#include <functional>

namespace molecule_descriptor
{

//Attribute must look the following :
//struct SampleProp
//{
//	static int id;
//	static bool registered;
//	typedef int64_t value_type;
//};
//
//int SampleProp::id = -1;
//bool SampleProp::registered = false;

//macros for properties
#define INT64_PROP(PropName)					struct PropName { typedef int64_t	value_type; static const bool kIsBounded = false;}
#define INT_PROP(PropName)						struct PropName { typedef int		value_type; static const bool kIsBounded = false; }
#define SIZE_T_PROP(PropName)					struct PropName { typedef size_t	value_type; static const bool kIsBounded = false; }
#define DOUBLE_PROP(PropName)					struct PropName { typedef double	value_type; static const bool kIsBounded = false; } 
#define BOOL_PROP(PropName)						struct PropName { typedef bool		value_type; static const bool kIsBounded = false; }
#define PTR_PROP(PropName, PtrType)				struct PropName { typedef PtrType*	value_type; static const bool kIsBounded = false; }
#define TEMPLATE_PTR_PROP(PropName)				template <typename PtrType> \
												struct PropName { typedef PtrType*	value_type; static const bool kIsBounded = false; }
#define PTR_TEMPLATE_PROP(PropName, PtrType)	template <typename T> \
												struct PropName { typedef PtrType<T>* value_type; static const bool kIsBounded = false; }

#define ENUM_TYPE(TypeName, ValueType, min_value, max_value) struct TypeName { \
															static_assert(std::is_integral<ValueType>::value, "type must be integral"); \
															typedef ValueType value_type; \
															static const bool kIsBounded = true; \
															static const ValueType kMinVal = min_value; \
															static const ValueType kMaxVal = max_value; }



using namespace std;

//for keeping numerical values in different items

union PropsUnion
{
	int64_t		value_int_64;
	int			value_int;
	char		value_char;
	unsigned char	value_uchar;
	double		value_double;
	void*		value_ptr;
	/*const*/ void *  value_const_ptr;
	bool		value_bool;
	size_t		value_size_t;
};
//getting values is implemented through classes because of absence of partial specialization of member functions
template <typename T> 
struct PropGetter
{
	static inline T& Value(PropsUnion&);
};

template <> 
struct PropGetter<int64_t>
{
	static inline int64_t& Value(PropsUnion& prop) 
	{ 
		return prop.value_int_64;
	}
};

template <> 
struct PropGetter<int>
{
	static inline int& Value(PropsUnion& prop) 
	{ 
		return prop.value_int;
	}
};

template <> 
struct PropGetter<double>
{
	static inline double& Value(PropsUnion& prop) 
	{ 
		return prop.value_double;
	}
};

template <> 
struct PropGetter<bool>
{
	static inline bool& Value(PropsUnion& prop) 
	{ 
		return prop.value_bool;
	}
};

template <> 
struct PropGetter<size_t>
{
	static inline size_t& Value(PropsUnion& prop) 
	{ 
		return prop.value_size_t;
	}
};

template <> 
struct PropGetter<char>
{
	static inline char& Value(PropsUnion& prop) 
	{ 
		return prop.value_char;
	}
};

template <> 
struct PropGetter<unsigned char>
{
	static inline unsigned char& Value(PropsUnion& prop) 
	{ 
		return prop.value_uchar;
	}
};

template <typename T> 
struct PropGetter<T*>
{
	static inline T*& Value(PropsUnion& prop) 
	{ 
		return reinterpret_cast<T*&>(prop.value_ptr);
	}
};

template <typename T> 
struct PropGetter<const T*>
{
	typedef const T const_type;
	typedef const_type* const_ptr_type;
	typedef T* ptr_type;

	static const T*& Value(PropsUnion& prop) 
	{ 
		return const_cast<const T*&>(reinterpret_cast<T*&>(prop.value_const_ptr));
	}
};

//map based on vector because vector is faster
class AttributesContainer
{
public:
	AttributesContainer() : m_inner_size(0) { }
	template <class AttributeType>
	typename AttributeType::value_type& Get();
	template <class AttributeType>
	void Delete();
	template <class AttributeType>
	void Add();

	typedef std::pair<std::type_index, PropsUnion> elem_type;//first element is property id
	typedef std::vector<elem_type> container_type;
private:
	container_type::iterator End() {
		return m_container.begin() + m_inner_size; }
	container_type m_container;
	PropsUnion m_temp;
	size_t m_inner_size;
};

//compare index in vector with property id 
template <class AttributeType, template <typename> class Comparator>
struct TypeIdComparator
{
	inline bool operator()(pair<std::type_index, PropsUnion>& element)
	{
		return m_comparator(element.first, std::type_index(typeid(AttributeType())));
	}

	Comparator<std::type_index> m_comparator; 
};

template <class AttributeType, bool kIsBounded = AttributeType::kIsBounded>
struct InitPropValue;

template <class AttributeType>
struct InitPropValue <AttributeType, true>
{
	inline void operator()(PropsUnion& val)
	{
		PropGetter<typename AttributeType::value_type>::Value(val) = 
			AttributeType::kMinVal;
	}
};

template <class AttributeType>
struct InitPropValue <AttributeType, false>{ 
	inline void operator()(PropsUnion& val){
		typedef typename AttributeType::value_type value_type;
		PropGetter<value_type>::Value(val) = 0;
	}	
};

template <template <typename> class Comparator>
struct CompareFirstInPair
{
	template <typename T1, typename T2>
	bool operator()(const std::pair<T1, T2>& elem_1, const std::pair<T1, T2>& elem_2)
	{
		return Comparator<T1>()(elem_1.first, elem_2.first); 
	}
};
template <class AttributeType>
void AttributesContainer::Add()
{
	auto res = find_if(m_container.begin(), End(), TypeIdComparator<AttributeType, std::equal_to>());

	if (res == End())
	{
		PropsUnion elem_to_push;
		InitPropValue<AttributeType>()(elem_to_push);

		if (End() == m_container.end())
		{
			m_container.push_back(make_pair(std::type_index(typeid(AttributeType())), elem_to_push));
		}
		else
		{
			*res = std::make_pair(std::type_index(typeid(AttributeType())), elem_to_push);
		}

		++m_inner_size;
		std::sort(m_container.begin(), End(), CompareFirstInPair<std::less>());
	}
	else
	{
		//throw std::logic_error("Property is already added");
	}
}

template <class AttributeType, bool kIsBounded = AttributeType::kIsBounded>
struct ThrowIfNotInRange;

template <class AttributeType>
struct ThrowIfNotInRange <AttributeType, true>
{
	inline void operator()(PropsUnion& val)
	{
		typedef typename AttributeType::value_type value_type;
		value_type prop_val = PropGetter<value_type>::Value(val);

		if (prop_val > AttributeType::kMaxVal || prop_val < AttributeType::kMinVal)
		{
			throw std::logic_error("type value out of range");
		}
	}
};

template <class AttributeType>
struct ThrowIfNotInRange <AttributeType, false>{ inline void operator()(PropsUnion& val){ }	};

template <class AttributeType>
typename AttributeType::value_type& AttributesContainer::Get()
{
	auto res = std::lower_bound(m_container.begin(), End(), 
		std::make_pair(std::type_index(typeid(AttributeType())), m_temp),
		CompareFirstInPair<std::less>());

	if (res == End() || TypeIdComparator<AttributeType, std::not_equal_to>()(*res))
	{
		throw std::logic_error("Property is not added. Add property before");
	}
	else
	{
		ThrowIfNotInRange<AttributeType>()(res->second);
		return	PropGetter<typename AttributeType::value_type>::Value(res->second);		
	}
}

template <class AttributeType>
void AttributesContainer::Delete()
{
	//m_container.erase(
	auto res = remove_if(m_container.begin(), End(), TypeIdComparator<AttributeType, std::equal_to>()); 

	if (res != End())
	{
		--m_inner_size;
		std::sort(m_container.begin(), End(), CompareFirstInPair<std::less>());
	}
}

}
