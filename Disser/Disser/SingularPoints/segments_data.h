#pragma once

#include "opencv2/core/core.hpp"
#include "boost/graph/graph_traits.hpp"

namespace molecule_descriptor
{
template <typename GraphT>
struct SegmentData
{
	typedef GraphT graph_type;
	graph_type* graph;
	typedef typename boost::graph_traits<GraphT>::vertex_descriptor vertex_descriptor;
	cv::Mat_<double> dist_map;
	std::vector<vertex_descriptor> elements;
	std::vector<vertex_descriptor> borders;

};
}