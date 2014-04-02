#pragma once 

#include <algorithm>
#include <iterator>

#include <boost/mpl/if.hpp>
#include "boost/graph/properties.hpp"
#include "boost/graph/subgraph.hpp"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/property_maps/container_property_map.hpp>

namespace molecule_descriptor
{

template <typename GraphType>
struct is_specialization_of_subgraph
{
	static const bool value = false;
};
template <typename GraphType>
struct is_specialization_of_subgraph<boost::subgraph<GraphType>>
{
	static const bool value = true;
};
enum VertexOrEdge {VERTEX, EDGE};


template <class Graph, VertexOrEdge kVertOrEdge>
struct VertexOrEdgeChoose
{
	typedef typename boost::mpl::if_c<kVertOrEdge == VERTEX, 
		boost::vertex_index_t, boost::edge_index_t>::type			BaseTag;

	typedef typename boost::mpl::if_c<
		is_specialization_of_subgraph<Graph>::value, 
		boost::local_property<BaseTag>,	typename BaseTag>::type		IndexTag;
	typedef typename boost::property_map<Graph, IndexTag>::type		IndexMapType;

	static IndexTag Property(BaseTag tag)
	{
		return IndexTag(tag);
	}
};


template <typename Graph, typename Container, VertexOrEdge kVertOrEdge>
class ContPropMap
	: public boost::iterator_property_map<
	typename Container::iterator, 
	typename VertexOrEdgeChoose<Graph, kVertOrEdge>::IndexMapType>
{
public:
	typedef VertexOrEdgeChoose<Graph, kVertOrEdge>										Chooser;
	typedef typename Chooser::BaseTag													BaseTag;
	typedef typename Chooser::IndexTag													IndexTag;
	typedef typename Chooser::IndexMapType												IndexMapType;

	typedef typename boost::mpl::if_c<kVertOrEdge == VERTEX, 
		typename boost::graph_traits<Graph>::vertex_descriptor,
		typename boost::graph_traits<Graph>::edge_descriptor>::type						Descriptor;
	typedef typename std::iterator_traits<typename Container::iterator>::value_type		ContValueType;
	typedef typename Container::iterator												Iterator;
	typedef typename Container::const_iterator											ConstIterator;
	typedef boost::iterator_property_map<Iterator, IndexMapType>						Base;
	inline ContPropMap(): Base(), m_container() {}

	Iterator begin()
	{
		return m_container.begin();
	}
	Iterator end()
	{
		return m_container.end();
	}
	ConstIterator begin() const
	{
		return m_container.begin();
	}
	ConstIterator end() const
	{
		return m_container.end();
	}
	explicit inline ContPropMap(const Graph& g)
	{
		SetGraph(g);
	}	
	void SetGraph(const Graph& g)
	{
		if (kVertOrEdge == VERTEX)
		{
			m_container.resize(boost::num_vertices(g));
		}
		else
		{
			m_container.resize(boost::num_edges(g));
		}

		Clear();
		const IndexMapType index_map = get(Chooser::Property(BaseTag()), const_cast<Graph&>(g));
		(static_cast<Base*>(this))->operator=(Base(m_container.begin(), index_map));
	}
	void Clear()
	{
		std::fill(m_container.begin(), m_container.end(), ContValueType());
	}
private:
	Container m_container;
};

} // end namespace molecule_descriptor

