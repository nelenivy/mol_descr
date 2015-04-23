#pragma once
#include <string>
#include <vector>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "molecule_manager.h"
#include "../../SVM/mgua.h"

namespace molecule_descriptor
{

class SetManager
{
public:
	SetManager() : m_inited(false), m_processed(false) {}
	void Init(int argc, char** argv);
	void ProcessSet();
	void ProcessSingularPoints(const bool calculate);
	void ProcessPairs(const bool calculate);
	void ProcessPairsLevels(const bool calculate);
	void ProcessTriples(const bool calculate);
	void ProcessTriplesLevels(const bool calculate);
	void ProcessKernelSVMPoints();
	void ProcessKernelSVMPointsWithFiltering();
	void ProcessKernelSVMHistograms();
	void ProcessDescriptorsSVM();
	void ProcessSVMClassificationL0();
	void ProcessMGUASVMClassification(const int iter_num, const int descr_num);
	void ProcessMGUASVMRegression();
private:
	void ReadParamsFromCommandLine(int argc, char** argv);
	void FindOutMoleculesNum();
	void CalculateDistThresholdsKMeans();
	void CalculateThresholdLevels();
	std::string MakeMoleculePrefix(const int ind);
private:
	//params
	std::string m_mol_folder;
	std::string m_mol_prefix;
	int m_molecules_num;
	int m_distances_intervals;
	int m_distance_interval_levels;
	int m_distance_quantile_for_thresh;
	int m_mesh_levels_num;
	int m_charges_threshes_num;
	int m_lennard_jones_threshes_num;
	int m_areas_threshes_num;
	bool m_use_levels;
	bool m_use_triples;
	int m_pairs_levels_overlap;
	bool m_write_pairs;
	MoleculeManager m_molecule_manager;
	std::vector<double> m_distances;
	std::vector<double> m_dist_thresholds;

	std::vector<std::vector<double>> m_distances_levels;
	std::vector<std::vector<double>> m_dist_thresholds_levels;
	std::vector<double> m_high_thresholds;

	bool m_inited;
	bool m_processed;
	cv::Mat_<size_t> m_md_matrix;
	cv::Mat_<double> m_md_matrix_double;
	cv::Mat_<float> m_md_matrix_float;
	std::vector<float> m_for_large_matrix;
	std::vector<float> m_for_large_matrix_1;
	std::vector<float> m_for_large_matrix_2;

	MGUATrainer<unsigned int> m_mgua_trainer;
};

}