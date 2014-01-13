#pragma once
#ifndef	CONNECTED_COMPONENTS_SEGMENTATOR_H
#define CONNECTED_COMPONENTS_SEGMENTATOR_H

#include <vector>
#include <utility>
#include <stack>
#include "CommonUtilities/attributes_container.h"
#include "graph_structures.h"
#include "segmentation_types.h"

namespace molecule_descriptor
{
using namespace std;
using namespace cv;
/** 
* Class for connected components segmentation of labeled graph
*/
template <class NodeType>
class ConnectedComponentsSegmentator
{
public:
	ConnectedComponentsSegmentator() : m_segments_num(0){ };
	/** \brief Segment connected regions of predefined value*/ 	
	template <typename PropertyType>
	void SegmentImageValue(const vector<GraphNode<NodeType>>& graph, const typename PropertyType::value_type value,	
		const size_t size_thresh, size_t& segments_num);
	/** \brief Segment connected regions*/ 
	template <typename PropertyType>
	void SegmentImage(const vector<GraphNode<NodeType>>& graph, const size_t size_thresh, size_t& segments_num);
	/** \brief Segment regions around seeds which value differ from seed value less than threshold*/ 
	template <typename PropertyType>
	void SegmentImageSeeds(const vector<GraphNode<NodeType>>& graph, const size_t size_thresh, 
		const vector<GraphNode<NodeType>>& seeds, const typename PropertyType::value_type color_thresh, size_t& segments_num);
	/** \brief Delete connected segments which size is less than threshold in binary image*/
	template <typename PropertyType>
	void DeleteSmallSegments(const vector<GraphNode<NodeType>>& graph, 
		const typename PropertyType::value_type value, const size_t size_thresh);
private:
	ConnectedComponentsSegmentator(const ConnectedComponentsSegmentator&);
	ConnectedComponentsSegmentator& operator=(const ConnectedComponentsSegmentator&);
	/** use only pixels with predefined value or all pixels*/
	enum eSegmentValues { ONE_VALUE, ALL_VALUES };
	/** \brief Implementation of segmentation*/ 	
	template <typename PropertyType>
	void SegmentImageImpl(const vector<GraphNode<NodeType>>& graph, eSegmentValues use_value,
		const typename PropertyType::value_type value, const size_t size_thresh, 
		const typename PropertyType::value_type color_thresh, const vector<GraphNode<NodeType>>& seeds, 
		size_t& labels_num);
	/** \brief Implementation of one segment distinguishing*/ 	
	template <typename PropertyType>
	void PickOutOneSegment(const GraphNode<NodeType>& seed, const typename PropertyType::value_type value, 
		const typename PropertyType::value_type color_thresh);
	/** \brief Grows region. This function is used in segmentation*/
	template <typename PropertyType>
	void GrowRegion(const typename PropertyType::value_type value, 
		typename PropertyType::value_type color_thresh, const size_t curr_segm, const GraphNode<NodeType>& curr_node);
	/** \brief Deletes region.*/
	void DeleteRegion(const size_t curr_segm, const GraphNode<NodeType>& curr_node);

	size_t m_segments_num;						/**< current number of segments*/
	stack<const GraphNode<NodeType>*> m_nodes_stack;	/**< stack for segmentation. 
												 We put neighbors of the current pixel to it. 
												 When it becomes empty we treat current region as processed*/
	size_t	m_size_thresh;						/**< minimum segment size for the current segmentation*/
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
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::SegmentImageValue(const vector<GraphNode<NodeType>>& graph, const typename PropertyType::value_type value,	
	const size_t size_thresh, size_t& segments_num)
{
	vector<GraphNode<NodeType>> temp_empty_vector;
	const typename PropertyType::value_type color_thresh = 0;

	SegmentImageImpl<PropertyType>(graph, ONE_VALUE, value, size_thresh, color_thresh, temp_empty_vector, segments_num);
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
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::SegmentImage(const vector<GraphNode<NodeType>>& graph, const size_t size_thresh, size_t& segments_num)
{	
	vector<GraphNode<NodeType>> temp_empty_vector;
	const typename PropertyType::value_type color_thresh = 0;
	const typename PropertyType::value_type value = 0;
	SegmentImageImpl<PropertyType>(graph, ALL_VALUES, value, size_thresh, color_thresh, temp_empty_vector, segments_num);
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
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::SegmentImageSeeds(const vector<GraphNode<NodeType>>& graph, const size_t size_thresh, 
	const vector<GraphNode<NodeType>>& seeds, const typename PropertyType::value_type color_thresh, size_t& segments_num)
{
	const typename PropertyType::value_type = 0;
	SegmentImageImpl<PropertyType>(graph, ALL_VALUES, value, size_thresh, color_thresh, seeds, segments_num);
}

/**
**************************************************************************
* \param [in]		image					input image
* \param [in]		connectivity			4- or 8-connectivity
* \param [in]		size_thresh				if size of the region is less then threshold we delete it
**************************************************************************
*/
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::DeleteSmallSegments(const vector<GraphNode<NodeType>>& graph, 
	const typename PropertyType::value_type value, const size_t size_thresh)
{
	int segments_num = 0;
	SegmentImageValue<PropertyType>(graph, value, size_thresh, segments_num);

	for (auto node_it = graph.begin(); node_it != graph.end(); ++node_it)
	{
		if (node_it->attr.Get<SegmentNumProp>() == 0)
		{
			node_it->attr.Get<PropertyType>() = 0;
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
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::SegmentImageImpl(const vector<GraphNode<NodeType>>& graph, eSegmentValues use_value,
	const typename PropertyType::value_type value, const size_t size_thresh, 
	const typename PropertyType::value_type color_thresh, const vector<GraphNode<NodeType>>& seeds, 
	size_t& labels_num)
{
	typedef typename PropertyType::value_type value_type;

	m_size_thresh = size_thresh;
	m_segments_num = 0;
	const vector<GraphNode<NodeType>>& inner_seeds = seeds.empty() ? graph : seeds;
	//segments numbers are going from 1, not from 0
	//segmentation process
	for (auto iter = graph.begin(); iter != graph.end(); ++iter)
	{//mark all nodes as non-visited
		iter->attr.Add<VisitedProp>();
		iter->attr.Get<VisitedProp>() = 0;
		iter->attr.Add<SegmentNumProp>();
	}

	for (auto iter = inner_seeds.begin(); iter != inner_seeds.end(); ++iter)
	{
		value_type curr_value = iter->attr.Get<PropertyType>();

		if (use_value == ONE_VALUE && curr_value != value)
		{
			continue;
		}

		const GraphNode<NodeType>& seed = *iter; 
		PickOutOneSegment<PropertyType>(seed, curr_value, color_thresh);
	}

	labels_num = m_segments_num;
	m_segments_num = 0;

	for (auto iter = graph.begin(); iter != graph.end(); ++iter)
	{//mark all nodes as non-visited
		iter->attr.Delete<VisitedProp>();
	}
}

/**
**************************************************************************
* \param [in]		value					only pixels with color = value are considered
* \param [in]		seed					seed coordinate
* \param [in]		color_thresh			threshold for color difference between seed pixel and current pixel
**************************************************************************
*/
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::PickOutOneSegment(const GraphNode<NodeType>& seed, 
	const typename PropertyType::value_type value, const typename PropertyType::value_type color_thresh)
{
	if (seed.attr.Get<VisitedProp>() != 0)
	{
		return;
	}			
	//init new segment
	m_segments_num++;
	seed.attr.Get<SegmentNumProp>() = m_segments_num;
	m_nodes_stack.push(&seed);
	int curr_segm_size = 0;
	//region growing process
	while (!(m_nodes_stack.empty()))
	{
		const GraphNode<NodeType>& curr_node = *(m_nodes_stack.top());
		m_nodes_stack.pop();
		curr_segm_size++;
		GrowRegion<PropertyType>(value, color_thresh, m_segments_num, curr_node);
	}
	//if segment size is less than threshold we delete it
	if (curr_segm_size < m_size_thresh)
	{
		seed.attr.Get<SegmentNumProp>() = 0;
		m_nodes_stack.push(&seed);

		while (!(m_nodes_stack.empty()))
		{
			const GraphNode<NodeType>& curr_node = *(m_nodes_stack.top());
			m_nodes_stack.pop();
			DeleteRegion(m_segments_num, curr_node);
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
template <typename NodeType>
template <typename PropertyType>
void ConnectedComponentsSegmentator<NodeType>::GrowRegion(const typename PropertyType::value_type value, 
typename PropertyType::value_type color_thresh, const size_t curr_segm, const GraphNode<NodeType>& curr_node)
{
	//process current pixel
	curr_node.attr.Get<SegmentNumProp>() = m_segments_num;
	curr_node.attr.Get<VisitedProp>() = 1;

	//process neighbours
	for (auto neighb_it = curr_node.neighbours.begin(); neighb_it != curr_node.neighbours.end(); ++neighb_it)
	{	
		if (
			abs(static_cast<double>((*neighb_it)->attr.Get<PropertyType>() - value)) > static_cast<double>(color_thresh) 
			|| (*neighb_it)->attr.Get<VisitedProp>() != 0)
		{
			continue;
		}
		
		(*neighb_it)->attr.Get<SegmentNumProp>() = curr_segm;
		(*neighb_it)->attr.Get<VisitedProp>() = 1;	
		m_nodes_stack.push(*neighb_it);	
	}	
}

template <typename NodeType>
void ConnectedComponentsSegmentator<NodeType>::DeleteRegion(const size_t curr_segm, const GraphNode<NodeType>& curr_node)
{
	//process current pixel
	curr_node.attr.Get<SegmentNumProp>() = 0;	

	//process neighbours
	for (auto neighb_it = curr_node.neighbours.begin(); neighb_it != curr_node.neighbours.end(); ++neighb_it)
	{	
		if ((*neighb_it)->attr.Get<SegmentNumProp>() !=  curr_segm)
		{
			continue;
		}

		(*neighb_it)->attr.Get<SegmentNumProp>() = 0;
		m_nodes_stack.push(*neighb_it);	
	}	
}

}

#endif
