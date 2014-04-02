#pragma once
#ifndef	CONNECTED_COMPONENTS_SEGMENTATOR_H
#define CONNECTED_COMPONENTS_SEGMENTATOR_H

#include <vector>
#include <utility>
#include <stack>
#include <type_traits>
#include <boost/graph/graph_traits.hpp>
#include "GraphLib/array_property_map.h"

namespace molecule_descriptor
{
using namespace std;
using namespace cv;
/** 
* Class for connected components segmentation of labeled graph
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
class ConnectedComponentsSegmentator
{
public:
	typedef typename std::remove_const<
		typename std::remove_reference<
		decltype(PropertyMap()[0])>::type
		>::type PropType;
	typedef typename std::remove_const<
		typename std::remove_reference<
		decltype(SegmentsMap()[0])>::type
	>::type SegmType;
	typedef typename boost::graph_traits<GraphType>::vertex_descriptor VertexType;
	ConnectedComponentsSegmentator() : m_segments_num(0){ };
	/** \brief Segment connected regions of predefined value*/ 	
	void SegmentImageValue(const GraphType& graph, const PropertyMap& prop_map, const PropType value,	
		const size_t size_thresh, SegmentsMap& segm_map, size_t& segments_num);
	/** \brief Segment connected regions*/ 
	void SegmentImage(const GraphType& graph, const PropertyMap& prop_map, 	
		const size_t size_thresh, SegmentsMap& segm_map, size_t& segments_num);
	/** \brief Segment regions around seeds which value differ from seed value less than threshold*/ 
	void SegmentImageSeeds(const GraphType& graph, const PropertyMap& prop_map, const size_t size_thresh, 
		const vector<VertexType>& seeds, const PropType color_thresh, SegmentsMap& segm_map, size_t& segments_num);
	/** \brief Delete connected segments which size is less than threshold in binary image*/
	void DeleteSmallSegments(const GraphType& graph, const PropType val, const PropType value, 
		const PropType null_val, const size_t size_thresh, PropertyMap& prop_map);
private:
	ConnectedComponentsSegmentator(const ConnectedComponentsSegmentator&);
	ConnectedComponentsSegmentator& operator=(const ConnectedComponentsSegmentator&);
	/** use only pixels with predefined value or all pixels*/
	enum eSegmentValues { ONE_VALUE, ALL_VALUES };
	/** \brief Implementation of segmentation*/ 	
	void SegmentImageImpl(const GraphType& graph, const PropertyMap& prop_map, eSegmentValues use_value,
		const PropType value, const size_t size_thresh, const PropType color_thresh, 
		const vector<VertexType>& seeds, SegmentsMap& segm_map, size_t& labels_num);
	/** \brief Implementation of one segment distinguishing*/ 	
	void PickOutOneSegment(const VertexType seed, const PropType value, const PropType color_thresh);
	/** \brief Grows region. This function is used in segmentation*/
	void GrowRegion(const VertexType curr_node,	const PropType value, const PropType color_thresh, 
		const SegmType curr_segm);
	/** \brief Deletes region.*/
	void DeleteRegion(const VertexType curr_node, const SegmType curr_segm);

	const GraphType* m_curr_graph;
	const PropertyMap* m_prop_map;
	SegmentsMap* m_segm_map;
	ContPropMap<GraphType, std::vector<bool>, VERTEX> m_visited_map;
	size_t m_segments_num;						/**< current number of segments*/
	vector<VertexType> m_seeds;
	stack<VertexType> m_nodes_stack;	/**< stack for segmentation. 
												 We put neighbors of the current pixel to it. 
												 When it becomes empty we treat current region as processed*/
	size_t	m_size_thresh;						/**< minimum segment size for the current segmentation*/
	ContPropMap<GraphType, std::vector<size_t>, VERTEX> m_internal_segm_map;
};

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		value					only pixels with color = value are considered
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
* \param [out]		segmented_image			segmented image of type CV_32SC1
* \param [out]		segments_num			number of segments
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::SegmentImageValue(
	const GraphType& graph, const PropertyMap& prop_map, const PropType value,	
	const size_t size_thresh, SegmentsMap& segm_map, size_t& segments_num)
{
	typename boost::graph_traits<GraphType>::vertex_iterator beg_vert, last_vert;
	m_seeds.resize(num_vertices(graph));
	boost::tie(beg_vert, last_vert) = vertices(graph);
	std::copy(beg_vert, last_vert, m_seeds.begin());
	const PropType color_thresh = 0;

	SegmentImageImpl(graph, prop_map, ONE_VALUE, value, size_thresh, color_thresh, m_seeds, segm_map, segments_num);
}

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
* \param [out]		segmented_image			segmented image of type CV_32SC1
* \param [out]		segments_num			number of segments
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::SegmentImage(
	const GraphType& graph, const PropertyMap& prop_map, 	
	const size_t size_thresh, SegmentsMap& segm_map, size_t& segments_num)
{	
	typename graph_traits<GraphType>::vertex_iterator beg_vert, last_vert;
	m_seeds.resize(vertices_num(graph));
	boost::tie(beg_vert, last_vert) = vertices(graph);
	std::copy(beg_vert, last_vert, m_seeds.begin());

	const PropType color_thresh = 0;
	const PropType value = 0;
	SegmentImageImpl(graph, ALL_VALUES, value, size_thresh, color_thresh, m_seeds, segments_num);
}

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
* \param [in]		seeds					seeds coordinates
* \param [in]		color_thresh			threshold for color difference between seed pixel and current pixel
* \param [out]		segmented_image			segmented image of type CV_32SC1
* \param [out]		segments_num			number of segments
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::SegmentImageSeeds(
	const GraphType& graph, const PropertyMap& prop_map, const size_t size_thresh, 
	const vector<VertexType>& seeds, const PropType color_thresh, SegmentsMap& segm_map, size_t& segments_num)
{
	const PropType value = 0;
	SegmentImageImpl(graph, ALL_VALUES, value, size_thresh, color_thresh, seeds, segments_num);
}

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::DeleteSmallSegments(
	const GraphType& graph, const PropType value, const PropType val, 
	const PropType null_val, const size_t size_thresh, PropertyMap& prop_map)
{
	int segments_num = 0;
	m_internal_segm_map.SetGraph(graph);
	SegmentImageValue(graph, prop_map, val, size_thresh, m_internal_segm_map, segments_num);

	boost::graph_traits<GraphType>::vertex_iterator curr_node, last;
	for (boost::tie(curr_node, last) = vertices(graph); curr_node != last; ++curr_node)
	{		
		if (m_internal_segm_map[*curr_node] == 0)
		{
			prop_map[*curr_node] = null_val;
		}
	}
}

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		use_value				use parameter value or not
* \param [in]		value					only pixels with color = value are considered
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
* \param [in]		seeds					seeds coordinates
* \param [in]		color_thresh			threshold for color difference between seed pixel and current pixel
* \param [out]		segmented_image			segmented image of type CV_32SC1
* \param [out]		segments_num			number of segments
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::SegmentImageImpl(
	const GraphType& graph, const PropertyMap& prop_map, eSegmentValues use_value,
	const PropType value, const size_t size_thresh, const PropType color_thresh, 
	const vector<VertexType>& seeds, SegmentsMap& segm_map, size_t& labels_num)
{
	m_size_thresh = size_thresh;
	m_segments_num = 0;
	m_curr_graph = &graph;
	m_prop_map = &prop_map;
	m_segm_map = &segm_map;
	m_visited_map.SetGraph(*m_curr_graph);

	//segments numbers are going from 1, not from 0
	//segmentation process
	typedef typename boost::graph_traits<GraphType>::vertex_iterator;
	for (auto curr_vertex = vertices(*m_curr_graph).first, end_vertex = vertices(*m_curr_graph).second; 
		curr_vertex != end_vertex; ++curr_vertex)
	{//mark all nodes as non-visited
		(*m_segm_map)[*curr_vertex] = 0;
		m_visited_map[*curr_vertex] = false;		
	}

	for (auto iter = seeds.begin(); iter != seeds.end(); ++iter)
	{
		PropType curr_value = (*m_prop_map)[*iter];

		if (use_value == ONE_VALUE && curr_value != value)
		{
			continue;
		}

		const VertexType seed = *iter; 
		PickOutOneSegment(seed, curr_value, color_thresh);
	}

	labels_num = m_segments_num;
	m_segments_num = 0;
}

/**
**************************************************************************
* \param [in]		value					only pixels with color = value are considered
* \param [in]		seed					seed coordinate
* \param [in]		color_thresh			threshold for color difference between seed pixel and current pixel
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::PickOutOneSegment(
	const VertexType seed, const PropType value, const PropType color_thresh)
{
	if (m_visited_map[seed] == true)
	{
		return;
	}			
	//init new segment
	m_segments_num++;
	(*m_segm_map)[seed] = m_segments_num;
	m_nodes_stack.push(seed);
	int curr_segm_size = 0;
	//region growing process
	while (!(m_nodes_stack.empty()))
	{
		VertexType curr_node = m_nodes_stack.top();
		m_nodes_stack.pop();
		curr_segm_size++;
		GrowRegion(curr_node, value, color_thresh, m_segments_num);
	}
	//if segment size is less than threshold we delete it
	if (curr_segm_size < m_size_thresh)
	{
		(*m_segm_map)[seed] = 0;
		m_nodes_stack.push(seed);

		while (!(m_nodes_stack.empty()))
		{
			const VertexType curr_node = m_nodes_stack.top();
			m_nodes_stack.pop();
			DeleteRegion(curr_node, m_segments_num);
		}

		m_segments_num--;
	}		
}
/**
**************************************************************************
* \param [in]		value, color_thresh		we consider only pixels with |color - value| <= color_thresh
* \param [in]		curr_segm				a number of the new segment
* \param [in]       curr_pix				current pixel coordinates
**************************************************************************
*/
template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::GrowRegion(	
	const VertexType curr_node,	const PropType value, const PropType color_thresh, const SegmType curr_segm)
{
	//process current pixel
	(*m_segm_map)[curr_node] = m_segments_num;
	m_visited_map[curr_node] = true;

	//process neighbours
	for (auto neighb_it = adjacent_vertices(curr_node, *m_curr_graph).first,
		end_neighb = adjacent_vertices(curr_node, *m_curr_graph).second; 
		neighb_it != end_neighb; ++neighb_it)
	{	
		if (
			abs(static_cast<double>((*m_prop_map)[*neighb_it] - value)) > static_cast<double>(color_thresh) 
			|| m_visited_map[*neighb_it] == true)
		{
			continue;
		}
		
		(*m_segm_map)[*neighb_it] = m_segments_num;
		m_visited_map[*neighb_it] = true;
		m_nodes_stack.push(*neighb_it);	
	}	
}

template <class GraphType, typename PropertyMap, typename SegmentsMap>
void ConnectedComponentsSegmentator<GraphType, PropertyMap, SegmentsMap>::DeleteRegion(
	const VertexType curr_node, const SegmType curr_segm)
{
	//process current pixel
	(*m_segm_map)[curr_node] = 0;

	//process neighbours
	for (auto neighb_it = adjacent_vertices(curr_node, *m_curr_graph).first,
		end_neighb = adjacent_vertices(curr_node, *m_curr_graph).second; 
		neighb_it != end_neighb; ++neighb_it)
	{	
		if ((*m_segm_map)[*neighb_it] !=  curr_segm)
		{
			continue;
		}

		(*m_segm_map)[*neighb_it] = 0;
		m_nodes_stack.push(*neighb_it);	
	}	
}

}

#endif
