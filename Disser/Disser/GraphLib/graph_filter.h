#pragma once

#include <vector>
#include <algorithm>
#include <stdint.h>

#include "boost/property_map/property_map.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/subgraph.hpp"
#include "boost/graph/graph_traits.hpp"
#include "SingularPoints/mesh_types.h"
#include "GraphLib/array_property_map.h"
#include "GraphLib/connected_components_segmentator.h"
#include "GraphLib/graph_functions.h"

namespace molecule_descriptor
{

inline double Distance(const cv::Point3d& elem1, const cv::Point3d elem2)
{
	return norm(elem1 - elem2);
}

template <class Elem>
struct DistFunc
{
	inline double operator()(const Elem& elem_1, const Elem& elem_2)
	{
		return Distance(elem_1, elem_2);
	}
};

//////////////////////////////////////////////////////////////////////////
//KERNELS FOR FILTERING
//////////////////////////////////////////////////////////////////////////
template <typename Point3dT, typename PropT>
class GaussianKernel
{
public:
	explicit GaussianKernel(const double sigma) : m_center_vert_set(false), m_sigma(sigma){}
	void SetCentralVertice(const Point3dT& coord, const PropT& prop)
	{
		m_center = coord;
		m_props_sum = prop - prop;
		m_coeffs_sum = 0.0;
		m_center_vert_set = true;
	}
	void AddToSet(const Point3dT& coord, const PropT& prop)
	{
		CV_Assert(m_center_vert_set);
		const double dist = norm(coord - m_center);
		const double coeff = exp( -dist* dist / (2.0 * m_sigma * m_sigma));
		m_coeffs_sum += coeff ;
		m_props_sum += prop * coeff;
	}
	PropT GetResult()
	{
		CV_Assert(m_center_vert_set);
		return m_props_sum * (1.0 / m_coeffs_sum);
	}
	double GetRadius() const
	{
		return 3.0 * m_sigma;
	}
private:
	Point3dT m_center;
	PropT m_props_sum;
	double m_coeffs_sum;
	bool m_center_vert_set;
	double m_sigma;
};

template <typename PropT, typename RadiusType>
class MedianKernel
{
public:
	explicit MedianKernel(const RadiusType radius) : m_center_vert_set(false), m_radius(radius){}
	template <typename Point3dT> 
	void SetCentralVertice(const Point3dT&, const PropT& prop)
	{
		m_props_cont.clear();
		m_center_vert_set = true;
	}
	template <typename Point3dT> 
	void AddToSet(const Point3dT&, const PropT& prop)
	{
		CV_Assert(m_center_vert_set);
		m_props_cont.push_back(prop);
	}
	PropT GetResult()
	{
		CV_Assert(m_center_vert_set);
		CV_Assert(!m_props_cont.empty());
		std::nth_element(m_props_cont.begin(), m_props_cont.begin() + m_props_cont.size() / 2, m_props_cont.end());
		return m_props_cont[m_props_cont.size() / 2];
	}
	RadiusType GetRadius() const
	{
		return m_radius;
	}
private:
	std::vector<PropT> m_props_cont;
	RadiusType m_radius;
	bool m_center_vert_set;
};

template <typename PropT, typename RadiusType>
class AverageKernel
{
public:
	explicit AverageKernel(const RadiusType radius) : m_center_vert_set(false), m_radius(radius){}
	template <typename Point3dT> 
	void SetCentralVertice(const Point3dT&, const PropT& prop)
	{
		m_props_sum = prop - prop;
		m_elems_num = 0;
		m_center_vert_set = true;
	}
	template <typename Point3dT> 
	void AddToSet(const Point3dT&, const PropT& prop)
	{
		CV_Assert(m_center_vert_set);
		m_props_sum += prop;
		++m_elems_num;
	}
	PropT GetResult()
	{
		CV_Assert(m_center_vert_set);
		CV_Assert(m_elems_num > 0.0);
		return m_props_sum * (1.0 / m_elems_num);
	}
	RadiusType GetRadius() const
	{
		return m_radius;
	}
private:
	PropT m_props_sum;
	double m_elems_num;
	RadiusType m_radius;
	bool m_center_vert_set;
};

//////////////////////////////////////////////////////////////////////////
//FILTERS
//////////////////////////////////////////////////////////////////////////
template <class GraphType, class CoordinateMapIn, class Kernel>
class GraphDistFilter
{
public:
	typedef typename boost::property_traits<CoordinateMapIn>::value_type CoordType;
	explicit GraphDistFilter(Kernel kernel = Kernel()) : m_segmentator(DistFunc<CoordType>()), m_kernel(kernel) {}
	void SetKernel(const Kernel& kernel)
	{
		m_kernel = kernel;
	}
	template <class PropMapIn, class PropMapOut>
	void Filter(const GraphType& graph, const CoordinateMapIn& coord_in,
		const std::vector<PropMapIn>& prop_map_in, std::vector<PropMapOut>& prop_map_out)
	{
		CV_Assert(prop_map_in.size() == prop_map_out.size());
		const double max_rad = m_kernel.GetRadius();
		typedef boost::graph_traits<GraphType>::vertex_descriptor vertex_descriptor;
		std::vector<vertex_descriptor> seed(1);
		m_mask.SetGraph(graph);
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second;
			curr_vert != end_vert; ++curr_vert)
		{
			seed[0] = *curr_vert;
			size_t segm_num = 0;
			m_segmentator.SegmentImageSeeds(graph, coord_in, 0, seed, max_rad, m_mask, segm_num);
			CV_Assert(segm_num > 0);

			for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
			{
				m_kernel.SetCentralVertice(coord_in[*curr_vert], prop_map_in[curr_prop][*curr_vert]);
				for (auto curr_vert_coord = vertices(graph).first, end_vert_coord = vertices(graph).second;
					curr_vert_coord != end_vert_coord; ++curr_vert_coord)
				{
					if (m_mask[*curr_vert_coord] == 0)
					{
						continue;
					}
					m_kernel.AddToSet(coord_in[*curr_vert_coord], prop_map_in[curr_prop][*curr_vert_coord]);
				}

				prop_map_out[curr_prop][*curr_vert] = m_kernel.GetResult();
			}
		}	
	}
	template <class PropMapIn, class PropMapOut>
	void Filter(const GraphType& graph, const CoordinateMapIn& coord_in,
		const PropMapIn& prop_map_in, PropMapOut& prop_map_out)
	{
		std::vector<ProxyPropMap<const PropMapIn, IdenticalTransformFunc>> prop_map_in_vect_ref(1, GetPropMapReference(prop_map_in));
		std::vector<ProxyPropMap<PropMapOut, IdenticalTransformFunc>> prop_map_out_vect_ref(1, GetPropMapReference(prop_map_out));
		Filter(graph, coord_in,	prop_map_in_vect_ref, prop_map_out_vect_ref);		
	}
private:
	typedef ContPropMap<GraphType, std::vector<uint8_t>, VERTEX> MaskType;
	MaskType m_mask;
	ConnectedComponentsSegmentator<GraphType, CoordinateMapIn, MaskType, double, DistFunc<CoordType>> m_segmentator;
	Kernel m_kernel;
};

template <class GraphType, class CoordinateMapIn, class Kernel, class PropMapIn, class PropMapOut>
void FilterGraphDist(const Kernel& kernel, const GraphType& graph, const CoordinateMapIn& coord_in,
			const PropMapIn& prop_map_in, PropMapOut& prop_map_out)
{
	GraphDistFilter<GraphType, CoordinateMapIn, Kernel> graph_filter(kernel);
	graph_filter.Filter(graph, coord_in, prop_map_in, prop_map_out);
}
//////////////////////////////////////////////////////////////////////////
template <class GraphType, class Kernel>
class GraphDistEdgesFilter
{
public:
	explicit GraphDistEdgesFilter(Kernel kernel = Kernel()) : m_kernel(kernel) {}
	void SetKernel(const Kernel& kernel)
	{
		m_kernel = kernel;
	}
	template <class PropMapIn, class PropMapOut>
	void Filter(const GraphType& graph, const std::vector<PropMapIn>& prop_map_in, std::vector<PropMapOut>& prop_map_out)
	{
		CV_Assert(prop_map_in.size() == prop_map_out.size());
		const size_t max_rad = m_kernel.GetRadius();
		typedef boost::graph_traits<GraphType>::vertex_descriptor vertex_descriptor;
		m_mask.SetGraph(graph);
		for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second;
			curr_vert != end_vert; ++curr_vert)
		{
			const vertex_descriptor seed = *curr_vert;
			WaveAlgorithm(graph, max_rad, seed, m_mask);

			for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
			{
				m_kernel.SetCentralVertice(*curr_vert, prop_map_in[curr_prop][*curr_vert]);
				for (auto curr_vert_coord = vertices(graph).first, end_vert_coord = vertices(graph).second;
					curr_vert_coord != end_vert_coord; ++curr_vert_coord)
				{
					if (m_mask[*curr_vert_coord] == 0)
					{
						continue;
					}
					m_kernel.AddToSet(*curr_vert_coord, prop_map_in[curr_prop][*curr_vert_coord]);
				}

				prop_map_out[curr_prop][*curr_vert] = m_kernel.GetResult();
			}
		}	
	}
	template <class PropMapIn, class PropMapOut>
	void Filter(const GraphType& graph, const PropMapIn& prop_map_in, PropMapOut& prop_map_out)
	{
		std::vector<ProxyPropMap<const PropMapIn, IdenticalTransformFunc>> prop_map_in_vect_ref(1, GetPropMapReference(prop_map_in));
		std::vector<ProxyPropMap<PropMapOut, IdenticalTransformFunc>> prop_map_out_vect_ref(1, GetPropMapReference(prop_map_out));
		Filter(graph, prop_map_in_vect_ref, prop_map_out_vect_ref);
	}
private:
	typedef ContPropMap<GraphType, std::vector<uint8_t>, VERTEX> MaskType;
	MaskType m_mask;
	Kernel m_kernel;
};

template <class GraphType, class Kernel, class PropMapIn, class PropMapOut>
void FilterGraphEdgeDist(const Kernel& kernel, const GraphType& graph, const PropMapIn& prop_map_in, PropMapOut& prop_map_out)
{
	GraphDistEdgesFilter<GraphType, Kernel> graph_filter(kernel);
	graph_filter.Filter(graph, prop_map_in, prop_map_out);
}
//////////////////////////////////////////////////////////////////////////
template <class Graph, class SegmentsMap, class PropsMap, class Kernel>
void CalculateAveragePropUsingKernel(const Graph& graph, const SegmentsMap& segments_map, const PropsMap& props_map, Kernel& kernel, 
						  std::vector<
						  std::pair<typename boost::property_traits<PropsMap>::value_type, bool>>& average_values)
{
	//find segments number
	int max_segm = - 1;
	for (auto vert_iter = vertices(graph).first, end_vert = vertices(graph).second; vert_iter != end_vert; ++vert_iter)
	{
		CV_Assert(segments_map[*vert_iter] >= 0);
		max_segm = std::max(max_segm, static_cast<int>(segments_map[*vert_iter]));
	}
	average_values.clear();
	if (max_segm <= 0)
	{
		return;
	}

	average_values.resize(max_segm + 1);
	std::vector<Kernel> kernels(max_segm + 1, kernel);
	std::vector<bool> segm_is_processed(max_segm + 1, false);

	for (auto vert_iter = vertices(graph).first, end_vert = vertices(graph).second; vert_iter != end_vert; ++vert_iter)
	{
		typedef typename boost::property_traits<PropsMap>::value_type PropT;
		const int curr_segm = static_cast<int>(segments_map[*vert_iter]);
		const PropT curr_prop = props_map[*vert_iter];

		if (!segm_is_processed[curr_segm])
		{
			segm_is_processed[curr_segm] = true;
			kernels[curr_segm].SetCentralVertice(curr_prop, curr_prop);
		}
		kernels[curr_segm].AddToSet(curr_prop, curr_prop);
	}

	for (size_t curr_segm = 0; curr_segm < average_values.size(); ++curr_segm)
	{
		if (!segm_is_processed[curr_segm])
		{
			average_values[curr_segm].second = false;
		}
		else
		{
			average_values[curr_segm].second = true;
			average_values[curr_segm].first = kernels[curr_segm].GetResult();
		}
	}
}

template <class Graph, class SegmentsMap, class PropsMap>
void CalculateAverageProp(const Graph& graph, const SegmentsMap& segments_map, const PropsMap& props_map, 
									 std::vector<
									 std::pair<typename boost::property_traits<PropsMap>::value_type, bool>>& average_values)
{
	AverageKernel<double, int> kernel(0);
	CalculateAveragePropUsingKernel(graph, segments_map, props_map,kernel, average_values);
}

template <class Graph, class SegmentsMap, class PropsMap>
void CalculateMedianProp(const Graph& graph, const SegmentsMap& segments_map, const PropsMap& props_map, 
						  std::vector<
						  std::pair<typename boost::property_traits<PropsMap>::value_type, bool>>& average_values)
{
	MedianKernel<typename boost::property_traits<PropsMap>::value_type, int> kernel(0);
	CalculateAveragePropUsingKernel(graph, segments_map, props_map,kernel, average_values);
}

template <typename GraphType, typename CoordMap>
class AverageFinder
{
public:
	typedef typename boost::property_traits<CoordMap>::value_type CoordType;
	typedef typename boost::graph_traits<GraphType>::vertex_descriptor vertex_descriptor;
	template <typename PropMap, typename Kernel>
	typename boost::property_traits<PropMap>::value_type
		GetAverageInRad(const GraphType& graph, const PropMap& prop_map, const CoordMap& coord_map, Kernel kernel, const vertex_descriptor vert)
	{
		typedef typename boost::property_traits<PropMap>::value_type PropType;
		std::vector<vertex_descriptor> seed(1, vert);
		m_mask.SetGraph(graph);
		size_t segm_num = 0;
		m_segmentator.SegmentImageSeeds(graph, coord_map, 0, seed, kernel.GetRadius(), m_mask, segm_num);
		CV_Assert(segm_num > 0);
		kernel.SetCentralVertice(coord_map[vert], prop_map[vert]);
		for (auto curr_vert_coord = vertices(graph).first, end_vert_coord = vertices(graph).second;
			curr_vert_coord != end_vert_coord; ++curr_vert_coord)
		{
			if (m_mask[*curr_vert_coord] == 0)
			{
				continue;
			}
			kernel.AddToSet(coord_map[*curr_vert_coord], prop_map[*curr_vert_coord]);
		}

		return kernel.GetResult();
	}
private:
	typedef ContPropMap<GraphType, std::vector<uint8_t>, VERTEX> MaskType;
	MaskType m_mask;
	ConnectedComponentsSegmentator<GraphType, CoordMap, MaskType, double, DistFunc<CoordType>> m_segmentator;
};

template <class Graph, class PropMapIn, class PropMapOut>
void SimpleLaplacian(const Graph& graph, const PropMapIn& prop_map_in, PropMapOut& prop_map_out)
{
	for (auto curr_vert = vertices(graph).first, end_vert = vertices(graph).second;
		curr_vert != end_vert; ++curr_vert)
	{
		for (auto neighb_it = adjacent_vertices(*curr_vert, graph).first,
			end_neighb = adjacent_vertices(*curr_vert, graph).second; neighb_it != end_neighb; ++neighb_it)
		{
			prop_map_out[*curr_vert] = prop_map_in[*curr_vert] - prop_map_in[*neighb_it];
		}

		prop_map_out[*curr_vert] /= static_cast<double>(std::distance(vertices(graph).first, vertices(graph).second));
	}
}
}