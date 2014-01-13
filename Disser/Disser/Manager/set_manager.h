#pragma once
#include <string>
#include <vector>

#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "molecule_manager.h"

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
	void ProcessTriples(const bool calculate);
	void ProcessKernelSVM();
	void ProcessDescriptorsSVM();
	void ProcessMGUASVM();
private:
	void ReadParamsFromCommandLine(int argc, char** argv);
	void FindOutMoleculesNum();
	void CalculateDistThresholdsKMeans();
	void CalculateDistThresholdsQuantils();
	std::string MakeMoleculePrefix(const int ind);
private:
	//params
	std::string m_mol_folder;
	std::string m_mol_prefix;
	int m_molecules_num;
	int m_distances_intervals;

	MoleculeManager m_molecule_manager;
	std::vector<double> m_distances;
	std::vector<double> m_dist_thresholds;

	bool m_inited;
	bool m_processed;
	cv::Mat_<size_t> m_md_matrix;
};

}