#include "set_manager.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include <utility>
#include <functional>

#include "opencv2/core/core.hpp"
#include "shark\LinAlg\BLAS\VectorMatrixType.h"
#include "shark/Models/Kernels/LinearKernel.h"

#include "extensions.h"
#include "InputOutput/params_reader.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"
#include "CommonUtilities/classic_k_means.h"

#include "../../SVM/kernel_svm_trainer_manager.h"
#include "../../SVM/descriptor_svm_trainer.h"
#include "../../SVM/cross_validation_svm_trainer.h"
#include "../../SVM/mgua.h"

#include "../../SVM/Kernels/basic_kernels.h"
#include "../../SVM/Kernels/kernels_for_points.h"
namespace molecule_descriptor
{

const bool kCalculate = true;
const bool kNotCalculate = false;

void SetManager::Init(int argc, char** argv)
{
	ReadParamsFromCommandLine(argc, argv);
	FindOutMoleculesNum();
}

void SetManager::ProcessSingularPoints(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	std::cout << "Singular Points Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(calculate);
	}
}

void SetManager::ProcessPairs(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	m_distances.clear();
	m_dist_thresholds.resize(m_distances_intervals);
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(kNotCalculate);
		m_molecule_manager.FindPairs(calculate);
		m_molecule_manager.AppendDistances(m_distances);
	}

	//CalculateDistThresholdsKMeans();
	CalculateDistThresholdsQuantils();
	m_molecule_manager.SetDistThresholds(m_dist_thresholds);
	const int pairs_types_num = static_cast<int>(m_molecule_manager.GetPairsTypeNum());
	m_md_matrix.create(m_molecules_num, pairs_types_num);
	//calculate distance types and write them to the matrix
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindPairs(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypes();
		m_molecule_manager.GetPairsHistogramm(m_md_matrix.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
}

void SetManager::ProcessTriples(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	std::cout << "Triples Processing\n";
	const int triples_types_num = static_cast<int>(m_molecule_manager.GetTriplesTypeNum());
	m_md_matrix.create(m_molecules_num, triples_types_num);

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(kNotCalculate);
		m_molecule_manager.FindTriples(calculate);
		m_molecule_manager.CalculateTriplesWithTypes();
		m_molecule_manager.GetTriplesHistogramm(m_md_matrix.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::TriplesMDMatrix());
}

void SetManager::ProcessKernelSVMPoints()
{
	std::vector<std::vector<NonMarkedSingularPoint>> pts(m_molecules_num);

	for (int ind = 0; ind < m_molecules_num; ind++)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(false);
		m_molecule_manager.GetNonMarkedSingularPoints(pts[ind]);
	}

	std::vector<int> labels;
	ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	typedef KernelForParametersSet<DiracKernel<size_t>, /*AlwaysOneKernel*/DiracKernel<int>, /*AlwaysOneKernel<float>*/ TriangleKernelThresholded> kernel_type;
	TriangleKernelThresholded/*AlwaysOneKernel<float>*/ always_one(0.040, -1.0, 1.0);
	kernel_type tuple_kernel(
		DiracKernel<size_t>(),/*AlwaysOneKernel*/DiracKernel<int>(), always_one);
	KernelSVMTrainerManager<PropertiesSet, GaussianKernelOneDim, kernel_type> svm_trainer;
	svm_trainer.SetData(pts, unsigned_labels);
	std::cout << tuple_kernel.parameterVector();
	svm_trainer.SetKernels(GaussianKernelOneDim(0.1), tuple_kernel);
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix);
}

void SetManager::ProcessKernelSVMHistograms()
{
	std::vector<std::vector<HistogramSingularPoint<9>>> pts(m_molecules_num);

	for (int ind = 0; ind < m_molecules_num; ind++)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(true);
		m_molecule_manager.GetHistogramSingularPoints(pts[ind]);
	}

	std::vector<int> labels;
	ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	typedef std::array<uint8_t, 9> prop_type;
	KernelSVMTrainerManager<prop_type, GaussianKernelOneDim, ScalarProductKernel<prop_type>> svm_trainer;
	std::vector<std::vector<SingularPoint<prop_type>>> pts_in(m_molecules_num);
	for (int ind = 0; ind < pts.size(); ++ind)
	{
		pts_in[ind].resize(pts[ind].size());
		for (int ind1 = 0; ind1 < pts[ind].size(); ++ind1)
		{
			pts_in[ind][ind1] = pts[ind][ind1];
		}
	}
	svm_trainer.SetData(pts_in, unsigned_labels);
	svm_trainer.SetKernels(GaussianKernelOneDim(0.1), ScalarProductKernel<prop_type>());
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix+ "_histo");
}

void SetManager::ProcessDescriptorsSVM()
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}


	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	DescriptorSVMTrainerManager svm_trainer;
	svm_trainer.SetData(m_md_matrix, unsigned_labels);
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix);
}

void SetManager::ProcessMGUASVM()
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	CrossValidationSvmTrainer<shark::RealVector> svm_trainer;
	std::vector<double> param_min(1, -7.0), param_max(1, 3.0);
	std::vector<size_t> param_num(1, 13);
	svm_trainer.SetGrid(param_min, param_max, param_num);
	svm_trainer.SetKernel(shark::LinearKernel<shark::RealVector>());
	svm_trainer.SetOutput(m_mol_folder + m_mol_prefix + "trainer_output.txt");
	MGUATrainer mgua_trainer;
	mgua_trainer.SetData(m_md_matrix, unsigned_labels);
	mgua_trainer.SetParameters(20, 30);
	mgua_trainer.Train(svm_trainer, m_mol_folder + m_mol_prefix);
}

std::string SetManager::MakeMoleculePrefix(const int ind)
{
	std::stringstream curr_file_prefix;
	curr_file_prefix << m_mol_folder << m_mol_prefix << "_" << ind;
	return curr_file_prefix.str();
}

void SetManager::CalculateDistThresholdsKMeans()
{
	const double null_elem = 0.0;
	const double dist_thresh = 0.1;
	const size_t max_iterations = 10000;
	std::vector<size_t> segmented;

	struct AbsMinus
	{
		double operator()(double arg1, double arg2) { 
			return abs(arg1 - arg2); }
	};

	std::function<double(double, double)> abs_minus_for_double = AbsMinus();

	SegmentKMeans(m_distances, m_distances_intervals, abs_minus_for_double, null_elem, 
		dist_thresh, max_iterations, 
		segmented, m_dist_thresholds);
	m_dist_thresholds.erase(m_dist_thresholds.begin());
	std::sort(m_dist_thresholds.begin(), m_dist_thresholds.end(), std::less<double>());
}

void SetManager::CalculateDistThresholdsQuantils()
{
	std::sort(m_distances.begin(), m_distances.end());
	m_dist_thresholds.resize(m_distances_intervals);
	const double elems_in_block = static_cast<double>(m_distances.size()) / (m_distances_intervals + 1);

	for (int block_ind = 1; block_ind <= m_distances_intervals; block_ind++)
	{
		const size_t curr_ind = std::min(Round(block_ind * elems_in_block), static_cast<int>(m_distances.size() - 1));
		m_dist_thresholds[block_ind - 1] = m_distances[curr_ind];
	}
}

void SetManager::ReadParamsFromCommandLine(int argc, char** argv)
{
	ReadParamFromCommandLine(argc, argv, "-mol_folder", m_mol_folder);

	if (!m_mol_folder.empty())
	{
		if (m_mol_folder.back() != '\\')
		{
			m_mol_folder.push_back('\\');
		}
	}

	ReadParamFromCommandLine(argc, argv, "-mol_prefix", m_mol_prefix);
	ReadParamFromCommandLine(argc, argv, "-distances_intervals", m_distances_intervals);
	m_inited = true;
}

void SetManager::FindOutMoleculesNum()
{
	m_molecules_num = 0;
	SurfaceReader surface_reader;
	std::stringstream curr_surf_file;

	do
	{
		curr_surf_file.str("");
		curr_surf_file.clear();
		curr_surf_file << m_mol_folder << m_mol_prefix << "_" << m_molecules_num << Extensions::Surface();
		++m_molecules_num;		
	} 
	while (surface_reader.OpenFile(curr_surf_file.str()));

	m_molecules_num--;
}

}