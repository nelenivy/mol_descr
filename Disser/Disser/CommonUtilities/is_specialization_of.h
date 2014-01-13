#pragma once

//test if T is a specilization of template
//template <typename T, template <typename...> class Template>
//struct IsSpecializationOf
//{
//	static const bool value = false;
//};
//
//template < template <typename...> class Template, typename... Types>
//struct IsSpecializationOf<Template<Types...>, Template>
//{
//	static const bool value = true;
//};