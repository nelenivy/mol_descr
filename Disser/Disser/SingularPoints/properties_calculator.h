#pragma once

#include <vector>
#include <utility>
#include "opencv2/core/core.hpp"
#include "mesh_types.h"

namespace molecule_descriptor
{

inline cv::Point3d CalculateElectrostaticForce(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	cv::Point3d force = 0.0;

	for (auto charges_iter = charges.begin(); charges_iter != charges.end(); ++charges_iter)
	{
		const cv::Point3d curr_vect = charges_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.00000001);
		const cv::Point3d curr_direction = curr_vect * inverse_dist;
		force += curr_direction * charges_iter->second * (inverse_dist * inverse_dist);
	}

	return force;
}

inline double CalculatePotential(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	double potential = 0.0;

	for (auto charges_iter = charges.begin(); charges_iter != charges.end(); ++charges_iter)
	{
		const cv::Point3d curr_vect = charges_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.0000000000001);
		potential += charges_iter->second * inverse_dist;
	}

	return potential;
}


inline cv::Point3d CalculateLennardJonesForce(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const double probe_radius)
{
	cv::Point3d force = 0.0;

	for (auto radius_iter = wdv_radii.begin(); radius_iter != wdv_radii.end(); ++radius_iter)
	{
		const cv::Point3d curr_vect = radius_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.00000001);
		const double sigma = (radius_iter->second + probe_radius) / 2.0;
		const double incl = pow(sigma * inverse_dist, 12.0) - pow(sigma * inverse_dist, 6.0);
		const cv::Point3d curr_direction = curr_vect * inverse_dist;
		force += curr_direction * incl * (inverse_dist * inverse_dist);
	}

	return force;
}

inline double CalculateLennardJonesPotential(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const double probe_radius)
{
	double potential = 0.0;

	for (auto radius_iter = wdv_radii.begin(); radius_iter != wdv_radii.end(); ++radius_iter)
	{
		const cv::Point3d curr_vect = radius_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.00000001);
		const double sigma = (radius_iter->second + probe_radius) / 2.0;
		const double incl = pow(sigma * inverse_dist, 12.0) - pow(sigma * inverse_dist, 6.0);
		potential += incl;
	}

	return potential;
}
template <typename PropMap, typename PropMap3d>
void CalculateAllPotentials(const std::vector<std::pair<cv::Point3d, double>>& charges, const Mesh& mesh, PropMap& vertex_charge_map, PropMap3d& vertex_charge_dir_map)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	vertex_charge_map.SetGraph(vertices_graph);
	vertex_charge_dir_map.SetGraph(vertices_graph);
	const double kProbeRadius = 1.53;
	for (auto vertex_iter = vertices(vertices_graph).first, 
		end_iter = vertices(vertices_graph).second; vertex_iter != end_iter; ++vertex_iter)
	{
		const Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Center(); 
		const Point3d curr_norm = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Normal();
		const double curr_potent = CalculatePotential(curr_coord/* + curr_norm * kProbeRadius*/, charges);
		const Point3d curr_force = CalculateElectrostaticForce(curr_coord + curr_norm * kProbeRadius, charges);
		vertex_charge_map[*vertex_iter] = curr_force.x * curr_norm.x + curr_force.y * curr_norm.y + curr_force.z * curr_norm.z;
		
		vertex_charge_dir_map[*vertex_iter] = curr_force;
		vertex_charge_dir_map[*vertex_iter].x /= cv::norm(curr_force);
		vertex_charge_dir_map[*vertex_iter].y /= cv::norm(curr_force);
		vertex_charge_dir_map[*vertex_iter].z /= cv::norm(curr_force);
	}
}

template <typename PropMap, typename PropMap3d>
void CalculateLennardJonesPotentials(const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const Mesh& mesh, PropMap& vertex_lennard_jones_map, 
	PropMap3d& vertex_lennard_jones_dir_map)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	vertex_lennard_jones_map.SetGraph(vertices_graph);
	vertex_lennard_jones_dir_map.SetGraph(vertices_graph);
	const double kProbeRadius = 1.53;

	for (auto vertex_iter = vertices(vertices_graph).first, 
		end_iter = vertices(vertices_graph).second; vertex_iter != end_iter; ++vertex_iter)
	{
		const Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Center(); 
		const Point3d curr_norm = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Normal();
		const cv::Point3d curr_force = CalculateLennardJonesForce(curr_coord + curr_norm * kProbeRadius, wdv_radii, kProbeRadius);
		vertex_lennard_jones_map[*vertex_iter] = CalculateLennardJonesPotential(curr_coord + curr_norm * kProbeRadius, wdv_radii, kProbeRadius);
		vertex_lennard_jones_dir_map[*vertex_iter] = curr_force;
		vertex_lennard_jones_dir_map[*vertex_iter].x /= cv::norm(curr_force);
		vertex_lennard_jones_dir_map[*vertex_iter].y /= cv::norm(curr_force);
		vertex_lennard_jones_dir_map[*vertex_iter].z /= cv::norm(curr_force);
	}
}


}