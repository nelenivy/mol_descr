#pragma once

#include <vector>
#include "CommonUtilities/attributes_container.h"

namespace molecule_descriptor
{

template <class ElemType>
struct GraphNode
{
	GraphNode() : element(nullptr) { }
	ElemType* element;
	std::vector<GraphNode<ElemType>*> neighbours;
	mutable AttributesContainer attr;
};

template <class PropertyType>
struct CopyPropFromNodeToElem
{
	template <class ElemType>
	void operator()(const GraphNode<ElemType>& node)
	{
		node.element->attr.Add<PropertyType>();
		node.element->attr.Get<PropertyType>() = node.attr.Get<PropertyType>();
	}
};

template <class PropertyType>
struct CopyPropFromElemToNode
{
	template <class ElemType>
	void operator()(const GraphNode<ElemType>& node)
	{
		node.attr.Add<PropertyType>();
		node.attr.Get<PropertyType>() = node.element->attr.Get<PropertyType>();
	}
}; 

}