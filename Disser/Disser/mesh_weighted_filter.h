#pragma once

#include <vector>
#include <memory>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "boost/graph/graph_traits.hpp"
#include "GraphLib/connected_components_segmentator.h"
#include "GraphLib/array_property_map.h"
#include "mesh_types.h"
namespace molecule_descriptor
{
	//////////////////////////////////////////////////////////////////////////
//KERNELS FOR FILTERING
//////////////////////////////////////////////////////////////////////////
class ExpApproxForGauss
{
public:
	ExpApproxForGauss(const size_t values_num)
	{
		m_values.resize(values_num);
		const double max_val = 4.5;

		for (size_t ind = 0; ind < values_num; ++ind)
		{
			m_values[ind] = exp(- static_cast<double>(ind) / values_num * max_val);
		}
	}
	double GetGauss(const double sqr_dist_divided_by_2_sigma_sqr)
	{
		const double max_val = 4.5;
		const int ind = std::max(Round(std::min(max_val, sqr_dist_divided_by_2_sigma_sqr) / max_val
			* m_values.size()) - 1, 0);
		return m_values[ind];
	}
private:
	std::vector<double> m_values;
};
//
//inline double FastExp(double y) 
//{
//	static_assert(sizeof(double) == 2 * sizeof(int32_t), "");
//
//	union
//	{
//		double d;
//		int32_t i[2];
//	} u;
//	u.i[0] = 0;
//	u.i[1] = static_cast<int32_t>(1512775 * y + 1072632447);
//	return u.d;
//}
template <typename PropT>
class GaussianKernelWeightedDist
{
public:
	explicit GaussianKernelWeightedDist(const double sigma) : m_sigma(sigma){}
	void SetCentralVertice(const PropT& prop)
	{
		m_props_sum = prop - prop;
		m_coeffs_sum = 0.0;
	}
	void AddToSet(const double dist, const double weight, const PropT& prop)
	{
		const double coeff = weight * exp/*FastExp*/( -dist* dist / (2.0 * m_sigma * m_sigma));
		m_coeffs_sum += coeff ;
		m_props_sum += prop * coeff;
	}
	PropT GetResult()
	{
		return m_props_sum * (1.0 / m_coeffs_sum);
	}
	double GetRadius() const
	{
		return 3.0 * m_sigma;
	}
private:
	PropT m_props_sum;
	double m_coeffs_sum;
	double m_sigma;
};

template <typename PropT>
class MedianKernelWeightedDist
{
public:
	explicit MedianKernelWeightedDist(const double radius) : m_radius(radius){}
	void SetCentralVertice(const PropT& prop)
	{
		m_coeffs_sum = 0.0;
		m_collected_props.clear();
	}
	void AddToSet(const double dist, const double weight, const PropT& prop)
	{
		m_coeffs_sum += weight;
		PropWithWeight prop_with_weight = {prop, weight};
		m_collected_props.push_back(PropWithWeight(prop_with_weight));
	}
	PropT GetResult()
	{
		std::sort(m_collected_props.begin(), m_collected_props.end());
		const double half_coeff_sum = m_coeffs_sum / 2.0;
		
		auto curr_prop = m_collected_props.begin();
		for (double curr_coeff_sum = 0.0; curr_coeff_sum <= half_coeff_sum; ++curr_prop)
		{
			curr_coeff_sum += curr_prop->weight;
		}
		return curr_prop->prop;
	}
	double GetRadius() const
	{
		return m_radius;
	}
private:
	struct PropWithWeight
	{
		PropT prop;
		double weight;
		bool operator < (const PropWithWeight& prop_1)
		{
			return prop < prop_1.prop;
		}
	};
	std::vector<PropWithWeight> m_collected_props;
	double m_coeffs_sum;
	double m_radius;
};

template <typename PropT>
class MedianKernelDist
{
public:
	explicit MedianKernelDist(const double radius) : m_radius(radius){}
	void SetCentralVertice(const PropT& prop)
	{
		m_collected_props.clear();
	}
	void AddToSet(const double dist, const double weight, const PropT& prop)
	{
		m_collected_props.push_back(prop);
	}
	PropT GetResult()
	{
		const int mid_index = m_collected_props.size() / 2;
		std::nth_element(m_collected_props.begin(), m_collected_props.begin() + mid_index, m_collected_props.end());		
		return m_collected_props[mid_index];
	}
	double GetRadius() const
	{
		return m_radius;
	}
private:
	std::vector<PropT> m_collected_props;
	double m_radius;
};

template <typename PropT>
class LoGKernelWeightedDist
{
public:
	explicit LoGKernelWeightedDist(const double sigma)
		: m_sigma(sigma), m_sigma_sqr(sigma * sigma), 
		m_norm_coeff(1.0 / (M_PI * pow(sigma,4))){}
	void SetCentralVertice(const PropT& prop)
	{
		m_props_sum = prop - prop;
		m_coeffs_sum = 0.0;
	}
	void AddToSet(const double dist, const double weight, const PropT& prop)
	{
		const double arg = dist* dist / (2.0 * m_sigma_sqr);
		const double coeff = weight * m_norm_coeff* (1.0 - arg) * exp(-arg);
		m_coeffs_sum += abs(coeff) ;
		m_props_sum += prop * coeff;
	}
	PropT GetResult()
	{
		return m_props_sum * (1.0 / m_coeffs_sum);
	}
	double GetRadius() const
	{
		return 3.0 * m_sigma;
	}
private:
	PropT m_props_sum;
	double m_coeffs_sum;
	double m_sigma;
	double m_sigma_sqr;
	double m_norm_coeff;
};

template <typename PropT>
class GaussianKernelWeightedDistTable
{
public:
	explicit GaussianKernelWeightedDistTable(const double sigma, std::shared_ptr<ExpApproxForGauss> table) 
		: m_sigma(sigma), m_table(table){}
	void SetCentralVertice(const PropT& prop)
	{
		m_props_sum = prop-prop;
		m_coeffs_sum = 0.0;
	}
	void AddToSet(const double dist_from_cent, const double weight, const PropT& prop)
	{
		const double coeff = weight * m_table->GetGauss(dist_from_cent* dist_from_cent / (2.0 * m_sigma * m_sigma));
		m_coeffs_sum += coeff ;
		m_props_sum += prop * coeff;
	}
	PropT GetResult()
	{
		return m_props_sum * (1.0 / m_coeffs_sum);
	}
	double GetRadius() const
	{
		return 3.0 * m_sigma;
	}
private:
	PropT m_props_sum;
	double m_coeffs_sum;
	double m_sigma;
	std::shared_ptr<ExpApproxForGauss> m_table;
};

template <typename OutT>
struct SignedDistFunc
{
	template <typename T>
	inline OutT operator()(const T& elem1, const T& elem2)
	{
		return elem1 - elem2;
	}
};
template <typename PropT, typename DistFunc>
class LaplaceBeltramiKernelWeightedDist
{
public:
	LaplaceBeltramiKernelWeightedDist(const double sigma, 
		const DistFunc dist_func) : m_sigma(sigma), m_dist_func(dist_func){}
	void SetCentralVertice(const PropT& prop)
	{
		m_central_prop = prop;
		m_props_sum = prop - prop;
		m_coeffs_sum = 0.0;
	}
	void AddToSet(const double dist, const double weight, const PropT& prop)
	{
		const double coeff = weight * exp( -dist* dist / (2.0 * m_sigma * m_sigma));
		m_coeffs_sum += coeff ;
		m_props_sum += m_dist_func(prop, m_central_prop)* coeff;
	}
	PropT GetResult()
	{
		return m_props_sum * (1.0 / m_coeffs_sum);
	}
	double GetRadius() const
	{
		return 3.0 * m_sigma;
	}
private:
	PropT m_central_prop;
	PropT m_props_sum;
	double m_coeffs_sum;
	double m_sigma;
	DistFunc m_dist_func;
};
//FILTERS
//////////////////////////////////////////////////////////////////////////
template <class TrianglesGraph>
class MeshDistMapWeightedFilter
{
public:
	MeshDistMapWeightedFilter()
	{
		SetThreadsNum(1);
	}
	MeshDistMapWeightedFilter(const int threads_num)
	{
		SetThreadsNum(threads_num);
	}

	void SetThreadsNum(const int threads_num)
	{
		m_threads_num = threads_num;
		m_masks_all.resize(m_threads_num);
		m_segmentator_all.resize(m_threads_num);
		m_triangles_within_dist_all.resize(m_threads_num);
	}
	template <class VerticesGraph, class VertVertDistMap, class VertTrDistMap, 
			class Kernel, class PropMapIn, class PropMapOut>
	void FilterMeshWeighted(const VerticesGraph& vert_graph, const TrianglesGraph& tr_graph, 
							const VertVertDistMap& vert_vert_dist, const VertTrDistMap& vert_tr_dist, 
							const std::vector<PropMapIn>& prop_map_in, const bool use_central_vert,
							Kernel& kernel, std::vector<PropMapOut>& prop_map_out)
	{
		CV_Assert(prop_map_in.size() == prop_map_out.size());
		const double max_rad = kernel.GetRadius();
		typedef boost::graph_traits<VerticesGraph>::vertex_descriptor vertices_descriptor;
		typedef boost::graph_traits<TrianglesGraph>::vertex_descriptor triangles_descriptor;
		const auto trinagles_coords = get(boost::vertex_info_3d, tr_graph);
		std::vector<std::vector<triangles_descriptor>> seeds_all(m_threads_num);
		std::vector<std::vector<Kernel>> prop_kernels_all(m_threads_num);

		for (int curr_thread = 0; curr_thread < m_threads_num; ++curr_thread)
		{
			m_masks_all[curr_thread].SetGraph(tr_graph);
			m_triangles_within_dist_all[curr_thread].SetGraph(tr_graph);
			prop_kernels_all[curr_thread].resize(prop_map_in.size(), kernel);
		}	

		const int prev_threads_num = omp_get_num_threads();
		omp_set_num_threads(m_threads_num);
#pragma omp parallel for
		/*for (auto curr_vert = vertices(vert_graph).first, end_vert = vertices(vert_graph).second;
		curr_vert != end_vert; ++curr_vert)*/
		for (int curr_vert_num = 0;
			curr_vert_num < num_vertices(vert_graph); ++curr_vert_num)
		{
			const auto curr_vert = vertex(curr_vert_num, vert_graph);
			const int curr_thread = omp_get_thread_num();
			auto& seeds = seeds_all[curr_thread];
			auto& prop_kernels = prop_kernels_all[curr_thread];
			auto& mask = m_masks_all[curr_thread];
			auto& triangles_within_dist = m_triangles_within_dist_all[curr_thread];
			auto& segmentator = m_segmentator_all[curr_thread];

			seeds.clear();
			//mark triangles which are within distance
			for (auto curr_tr = vertices(tr_graph).first, end_tr = vertices(tr_graph).second;
				curr_tr != end_tr; ++curr_tr)
			{
				const MeshTriangle& curr_mesh_tr = trinagles_coords[*curr_tr];
				const bool adjacent_triangle = 
					(curr_mesh_tr.GetA() == curr_vert || curr_mesh_tr.GetB() == curr_vert || curr_mesh_tr.GetC() == curr_vert);
				
				if (adjacent_triangle || vert_tr_dist(curr_vert, *curr_tr) < max_rad)
				{
					triangles_within_dist[*curr_tr] = 255;
					if (adjacent_triangle)
					{
						seeds.push_back(*curr_tr);
					}
				}
				else
				{
					triangles_within_dist[*curr_tr] = 0;
				}
			}

			size_t segm_num = 0;
			segmentator.SegmentImageSeeds(tr_graph, triangles_within_dist, 0, seeds, 0, mask, segm_num);

			if (segm_num == 0)
			{
				for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
				{
					prop_map_out[curr_prop][curr_vert] = prop_map_in[curr_prop][curr_vert];
				}
			}
			else
			{
				for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
				{
					prop_kernels[curr_prop].SetCentralVertice(prop_map_in[curr_prop][curr_vert]);
				}

				for (auto curr_tr = vertices(tr_graph).first, end_tr = vertices(tr_graph).second;
					curr_tr != end_tr; ++curr_tr)
				{
					if (mask[*curr_tr] > 0)
					{
						const MeshTriangle curr_mesh_tr = trinagles_coords[*curr_tr];
						const double weight = trinagles_coords[*curr_tr].Area();
						
						for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
						{
							//const double mean_prop = (prop_map_in[curr_prop][curr_mesh_tr.GetA()] + 
							//prop_map_in[curr_prop][curr_mesh_tr.GetB()] + prop_map_in[curr_prop][curr_mesh_tr.GetC()]) / 3.0;

						//prop_kernels[curr_prop].AddToSet(vert_tr_dist(*curr_vert, *curr_tr), weight, mean_prop);
							if (curr_mesh_tr.GetA() != curr_vert || use_central_vert)
							{
								prop_kernels[curr_prop].AddToSet(vert_vert_dist(curr_vert, curr_mesh_tr.GetA()), 
									weight, prop_map_in[curr_prop][curr_mesh_tr.GetA()]);
							}
							if (curr_mesh_tr.GetB() != curr_vert || use_central_vert)
							{
								prop_kernels[curr_prop].AddToSet(vert_vert_dist(curr_vert, curr_mesh_tr.GetB()), 
									weight, prop_map_in[curr_prop][curr_mesh_tr.GetB()]);
							}
							if (curr_mesh_tr.GetC() != curr_vert || use_central_vert)
							{
								prop_kernels[curr_prop].AddToSet(vert_vert_dist(curr_vert, curr_mesh_tr.GetC()), 
									weight, prop_map_in[curr_prop][curr_mesh_tr.GetC()]);
							}
						}
					}
				}

				for (size_t curr_prop = 0; curr_prop < prop_map_in.size(); ++curr_prop)
				{
					prop_map_out[curr_prop][curr_vert] = prop_kernels[curr_prop].GetResult();
				}
			}
		}

		omp_set_num_threads(prev_threads_num);
	}
private:
	int m_threads_num;
	typedef ContPropMap<TrianglesGraph, std::vector<uint8_t>, VERTEX> MaskType;
	std::vector<MaskType> m_triangles_within_dist_all;
	std::vector<MaskType> m_masks_all;
	std::vector<ConnectedComponentsSegmentator<TrianglesGraph, MaskType, MaskType, int>> m_segmentator_all;
};


template <class TrianglesGraph, class VerticesGraph, class VertVertDistMap, class VertTrDistMap, 
class Kernel, class PropMapIn, class PropMapOut>
	void FilterMeshWeightedFunc(const VerticesGraph& vert_graph, const TrianglesGraph& tr_graph, 
	const VertVertDistMap& vert_vert_dist, const VertTrDistMap& vert_tr_dist, 
	const std::vector<PropMapIn>& prop_map_in, const bool use_central_vertice, Kernel& kernel, std::vector<PropMapOut>& prop_map_out)
{
	MeshDistMapWeightedFilter<TrianglesGraph> filter;
	filter.FilterMeshWeighted(vert_graph, tr_graph, vert_vert_dist, vert_tr_dist, 
		prop_map_in, use_central_vertice,kernel, prop_map_out);
}
template <class TrianglesGraph, class VerticesGraph, class VertVertDistMap, class VertTrDistMap, 
class Kernel, class PropMapIn, class PropMapOut>
	void FilterMeshWeightedFuncMultiThread(const VerticesGraph& vert_graph, const TrianglesGraph& tr_graph, 
	const VertVertDistMap& vert_vert_dist, const VertTrDistMap& vert_tr_dist, 
	const std::vector<PropMapIn>& prop_map_in, const bool use_central_vertice, const int threads_num, Kernel& kernel, std::vector<PropMapOut>& prop_map_out)
{
	MeshDistMapWeightedFilter<TrianglesGraph> filter(threads_num);
	filter.FilterMeshWeighted(vert_graph, tr_graph, vert_vert_dist, vert_tr_dist, 
		prop_map_in, use_central_vertice,kernel, prop_map_out);
}

}