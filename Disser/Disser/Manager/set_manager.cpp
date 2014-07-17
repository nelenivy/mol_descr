#include "set_manager.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include <utility>
#include <functional>
#include <limits>

#include "opencv2/core/core.hpp"
#include "shark\LinAlg\BLAS\matrix.hpp"
#include "shark/Models/Kernels/LinearKernel.h"

#include "extensions.h"
#include "InputOutput/params_reader.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"
#include "CommonUtilities/classic_k_means.h"
#include "types_calculation.h"

#include "../../SVM/kernel_svm_trainer_manager.h"
#include "../../SVM/descriptor_svm_trainer.h"
#include "../../SVM/cross_validation_svm_trainer.h"
#include "../../SVM/pharm_kernel_naive_configure.h"

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
	m_molecule_manager.SetLevelsNum(m_mesh_levels_num);
}

void SetManager::ProcessSingularPoints(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	std::cout << "Singular Points Processing\n";
	std::vector<double> charges, charges_thresh;//for charges classification
	std::vector<double> lennard_jones, lennard_jones_thresh;//for lennard-jones classification
	std::vector<double> areas, areas_thresh;//for areas classification
	const bool calc_as_average = false;
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(calculate, calc_as_average);
		m_molecule_manager.AppendCharges(charges, false/*true*/);
		m_molecule_manager.AppendLennardJones(lennard_jones, false/*true*/);
		m_molecule_manager.AppendLennardJones(areas, false/*true*/);
	}

	CalculateThresholdsQuantiles(2, charges, charges_thresh);
	CalculateThresholdsQuantiles(2, lennard_jones, lennard_jones_thresh);
	CalculateThresholdsQuantiles(0, areas, areas_thresh);
	m_molecule_manager.SetChargesThresholds(charges_thresh);
	m_molecule_manager.SetLennardJonesThresholds(lennard_jones_thresh);
	m_molecule_manager.SetAreaThresholds(areas_thresh);

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ClassifySingularPoints();
	}
}

void SetManager::ProcessPairsLevels(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	m_distances_levels.clear();
	m_dist_thresholds.resize(m_distances_intervals);
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairsLevels(calculate);
		m_molecule_manager.AppendDistancesLevels(m_distances_levels);
	}

	//CalculateDistThresholdsKMeans();
	CalculateThresholdLevels();
	m_molecule_manager.SetDistLevelsThresholds(m_dist_thresholds_levels, m_high_thresholds);
	std::vector<double> levels_scales(m_dist_thresholds_levels.size());
	levels_scales[0] = 1.0;
	for (int ind = 0; ind < levels_scales.size(); ++ind)
	{
		levels_scales[ind] = 1.0;//static_cast<double>(m_distances_levels[0].size()) / m_distances_levels[ind].size();
	}
	m_molecule_manager.SetLevelsScales(levels_scales);
	const int pairs_types_num = static_cast<int>(m_molecule_manager.GetPairsTypeNumLevels());
	m_md_matrix_double.create(m_molecules_num, pairs_types_num);
	//calculate distance types and write them to the matrix
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairsLevels(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypesLevels();
		m_molecule_manager.GetPairsHistogrammLevels(m_md_matrix_double.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix_double, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrixLevels());
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
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairs(calculate);
		m_molecule_manager.AppendDistances(m_distances);
	}

	//CalculateDistThresholdsKMeans();
	CalculateThresholdsLevels(m_distances_intervals, 0.05, m_distances, m_dist_thresholds);
	//CalculateThresholdsQuantiles(m_distances_intervals, m_distances, m_dist_thresholds);
	m_molecule_manager.SetDistThresholds(m_dist_thresholds);
	const int pairs_types_num = static_cast<int>(m_molecule_manager.GetPairsTypeNum());
	m_md_matrix.create(m_molecules_num, pairs_types_num);
	//calculate distance types and write them to the matrix
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairs(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypes();
		m_molecule_manager.GetPairsHistogramm(m_md_matrix.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
}

void SetManager::CalculateThresholdLevels()
{
	const int levels_num = m_distances_levels.size();
	CV_Assert(levels_num > 1);
	m_high_thresholds.resize(levels_num, std::numeric_limits<double>::max());
	m_dist_thresholds_levels.resize(levels_num);
	m_high_thresholds[levels_num - 1] = std::numeric_limits<double>::max();
	CV_Assert(m_distance_quantile_for_thresh <= m_distance_interval_levels);

	for (int level = levels_num - 1; level >= 0; --level)
	{
		//delete too large distances
		const auto first_elem_higher_than_thresh = std::partition(m_distances_levels[level].begin(), m_distances_levels[level].end(),
			[&](const double dist)->bool
		{
			return dist < m_high_thresholds[level];
		});
		m_distances_levels[level].erase(first_elem_higher_than_thresh, m_distances_levels[level].end());
		std::vector<double> thresh_levels, thresh_quantiles;
		CalculateThresholdsLevels(m_distance_interval_levels, 0.02, m_distances_levels[level], thresh_levels);
		CalculateThresholdsQuantiles(m_distance_interval_levels, m_distances_levels[level], thresh_quantiles);
		const double alpha = 1.0;
		m_dist_thresholds_levels[level].resize(m_distance_interval_levels);
		for (size_t ind = 0; ind < m_distance_interval_levels; ++ind)
		{
			m_dist_thresholds_levels[level][ind] = alpha * thresh_levels[ind] + (1.0  - alpha) * thresh_quantiles[ind];
		}
		if (level > 0 /*&& level < 8*/)
		{
			m_high_thresholds[level - 1] = m_dist_thresholds_levels[level][m_distance_quantile_for_thresh];
		}
	}
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
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
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
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.GetNonMarkedSingularPoints(pts[ind]);
	}

	std::vector<int> labels;
	ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	typedef TriangleKernel/*AlwaysOneKernel<float>*/ PotentialKernel;
	typedef KernelForParametersSet<DiracKernel<size_t>, /*AlwaysOneKernel*/DiracKernel<int>, PotentialKernel> kernel_type;	
	PotentialKernel potential_kernel(0.2);//(0.020, -1.0, 1.0);
	kernel_type tuple_kernel(
		DiracKernel<size_t>(),/*AlwaysOneKernel*/DiracKernel<int>(), potential_kernel);
	KernelSVMTrainerManager<PropertiesSet, GaussianKernelOneDim, kernel_type> svm_trainer;
	svm_trainer.SetData(pts, unsigned_labels);
	std::cout << tuple_kernel.parameterVector();
	svm_trainer.SetKernels(GaussianKernelOneDim(m_dist_thresholds), tuple_kernel);
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix);
}

void SetManager::ProcessKernelSVMPointsWithFiltering()
{	//Calculate good descriptors
	ProcessMGUASVMClassification(1, 1);
	const size_t descriptors_arr[] = {1369, 1542, 1381};//,
	/*1295, 157, 1009, 1048, 118, 38,
	1038, 14, 128, 176, 872, 135,
	901, 37, 271, 1048, 118, 38, 
	1157, 180, 1009, 1048, 118, 38,
	185, 6, 116, 1048, 118, 38,
	180, 157, 1009, 1048, 118, 38 ,
	415, 157, 1009, 1048, 118, 38,
	1053, 14, 128, 176, 872, 135,
	1295, 258, 128, 176, 872, 135};*/
	
	const size_t kGoodSize = sizeof(descriptors_arr) / sizeof(descriptors_arr[0]);
	std::vector<bool> good_descriptors = m_mgua_trainer.GetGoodDescriptors();	
	//ProcessPairs must be called
	//read pairs
	CV_Assert(m_inited);
	//calculate all pairs
	//////////////////////////////////////////////////////////////////////////
	std::cout << "Pairs Processing\n";
	std::vector<NonMarkedSingularPoint> sing_pts;
	std::vector<std::vector<SingularPointsPair<PropertiesSet>>> sing_pts_pairs(m_molecules_num);
	std::vector<std::vector<size_t>> sing_pts_pairs_types(m_molecules_num);
	std::vector<size_t> curr_pairs_types;
	std::vector<SingularPointsPair<PropertiesSet>> curr_non_marked_pairs;
	std::vector<std::vector<size_t>> new_descr_mat_sparsed(m_molecules_num);
	std::vector<std::vector<size_t>> new_descr_mat(m_molecules_num);
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		new_descr_mat_sparsed[ind].resize(kGoodSize, 0);
		new_descr_mat[ind].resize(m_molecule_manager.GetPairsTypeNum(), 0);
		size_t good_pairs = 0;
		size_t bad_pairs = 0;
		std::vector<SingularPointsPair<PropertiesSet>>& curr_pairs = sing_pts_pairs[ind];
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.GetNonMarkedSingularPoints(sing_pts);
		m_molecule_manager.FindPairs(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypes();
		m_molecule_manager.GetPairsTypes(curr_pairs_types);
		m_molecule_manager.GetNonMarkedPairs(curr_non_marked_pairs);

		for (size_t pair_num = 0; pair_num < curr_pairs_types.size(); ++pair_num)
		{
			const size_t curr_descr = curr_pairs_types[pair_num];
			const size_t* pos = std::find(descriptors_arr, descriptors_arr + kGoodSize, curr_descr);
			new_descr_mat[ind][curr_descr]++;
			if (descriptors_arr + kGoodSize == pos)
			//if (!good_descriptors[curr_descr])
			{
				bad_pairs++;
				continue;
			}
			good_pairs++;
			sing_pts_pairs[ind].push_back(curr_non_marked_pairs[pair_num]);
			new_descr_mat_sparsed[ind][pos - descriptors_arr]++;
		}

		std::cout << "molecule num " << ind << " bad_pairs " << bad_pairs << " good_pairs " << good_pairs << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	//get trainer
	typedef /*AlwaysOneKernel<float>*/TriangleKernel PotentialKernel;
	typedef KernelForParametersSet<DiracKernel<size_t>, /*AlwaysOneKernel*/DiracKernel<int>, PotentialKernel> kernel_type;	
	PotentialKernel potential_kernel(0.03);//, -1.0, 1.0);
	kernel_type tuple_kernel(
		DiracKernel<size_t>(),/*AlwaysOneKernel*/DiracKernel<int>(), potential_kernel);
	auto trainer = GetNaiveKernelSVMTrainerForClassification<PropertiesSet>(tuple_kernel, GaussianKernelOneDim(0.3/*m_dist_thresholds*/));
	std::vector<double> min(2);
	std::vector<double> max(2);
	std::vector<size_t> sections(2);
	min[1] = -7; max[1] = 7; sections[1] = 15;  // regularization parameter C
	min[0] = 0.10; max[0] = 1.0; sections[0] = 10;   // kernel parameter gamma
	trainer.SetGrid(min, max, sections);
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}
	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	//convert data to shark format
	shark::Data<std::vector<SingularPointsPair<PropertiesSet>>> input_data = shark::createDataFromRange(sing_pts_pairs);
	shark::Data<unsigned int> input_labels = shark::createDataFromRange(unsigned_labels);
	auto labeled_data = shark::LabeledData<std::vector<SingularPointsPair<PropertiesSet>>, unsigned int>(input_data, input_labels);	
	std::vector<ElemWithIndexAndID<std::vector<SingularPointsPair<PropertiesSet>>>> index_elems;
	CreateSeqWithIndexAndIdFromSeq(sing_pts_pairs.begin(), sing_pts_pairs.end(), 0, index_elems);
	{
		std::ofstream t("m5.txt");
		std::ofstream t1("md_new.txt");
		for (size_t y = 0; y < labeled_data.numberOfElements(); ++y)
		{
			for (int ind = 0; ind < new_descr_mat_sparsed[y].size(); ++ind)
			{
				t << new_descr_mat_sparsed[y][ind] << " ";
			}
			t << "\n";
			for (int ind = 0; ind < new_descr_mat[y].size(); ++ind)
			{
				t1 << new_descr_mat[y][ind] << " ";
			}
			t1 << "\n";
			for (size_t x = 0; x < labeled_data.numberOfElements(); ++x)
			{
				t << y << " " << x << " " << (*trainer.m_kernel)(index_elems[y], index_elems[x]) << " ";
				t << "\n";
			}
		}
	}

	trainer.train(labeled_data, labeled_data.numberOfElements() / 2);
}

void SetManager::ProcessKernelSVMHistograms()
{
	std::vector<std::vector<HistogramSingularPoint<9>>> pts(m_molecules_num);

	for (int ind = 0; ind < m_molecules_num; ind++)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
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
	svm_trainer.SetKernels(GaussianKernelOneDim(m_dist_thresholds), ScalarProductKernel<prop_type>());
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

void SetManager::ProcessMGUASVMClassification(const int iter_num, const int descr_num)
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix_double, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrixLevels());
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	auto svm_trainer = GetSVMTrainerForLinearClassification<shark::RealVector>();
	typedef decltype(svm_trainer) SVMType;
	const int kThreadsNum = 7;
	std::vector<SVMType> svm_vect(kThreadsNum);
	std::vector<double> param_min(1, -2.0), param_max(1, 1.0);
	std::vector<size_t> param_num(1, 4);

	for (size_t ind = 0; ind < svm_vect.size(); ++ind)
	{
		svm_vect[ind] = GetSVMTrainerForLinearClassification<shark::RealVector>();
		svm_vect[ind].SetGrid(param_min, param_max, param_num);
		svm_vect[ind].SetOutput(m_mol_folder + m_mol_prefix + "trainer_output1.txt");
	}

	m_mgua_trainer = MGUATrainer<unsigned int>();
	m_mgua_trainer.SetData(m_md_matrix_double, unsigned_labels, true);
	m_mgua_trainer.SetParameters(iter_num, descr_num);
	//m_mgua_trainer.TrainOneSet(svm_vect[0], m_mol_folder + m_mol_prefix);
	m_mgua_trainer.Train(svm_vect, m_mol_folder + m_mol_prefix,true);
}

void SetManager::ProcessMGUASVMRegression()
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix_double, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
	//read labels
	std::vector<double> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Properties(), labels);
	double mean = 0;
	std::vector<double> abs_labels(labels);
	for (auto it = labels.begin(), it_end = labels.end(); it != it_end; ++it)
	{
		abs_labels[it - labels.begin()] = abs(*it);
	}
	std::nth_element(abs_labels.begin(), abs_labels.begin() + abs_labels.size() / 2, abs_labels.end());
	mean = abs_labels[abs_labels.size() / 2];//labels.size();
	std::cout << "mean = " << mean << "\n";
	auto svm_trainer = GetSVMTrainerForRegression<shark::RealVector>(mean / 10.0);
	typedef decltype(svm_trainer) SVMType;
	const int kThreadsNum = 7;
	std::vector<SVMType> svm_vect(kThreadsNum);
	std::vector<double> param_min(2, -7.0), param_max(2, 3.0);
	std::vector<size_t> param_num(2, 13);
	param_max[1] = log(mean);
	param_min[1] = log(mean / 100.0);
	param_num[1] = 5;

	for (size_t ind = 0; ind < svm_vect.size(); ++ind)
	{
		svm_vect[ind] = GetSVMTrainerForRegression<shark::RealVector>(mean / 10.0);
		svm_vect[ind].SetGrid(param_min, param_max, param_num);
		svm_vect[ind].SetOutput(m_mol_folder + m_mol_prefix + "trainer_output1.txt");
	}

	MGUATrainer<shark::RealVector> mgua_trainer;
	/*cv::Mat_<double> md_mat(m_md_matrix.size());
	std::copy(m_md_matrix.begin(), m_md_matrix.end(), md_mat.begin());*/
	std::vector<shark::RealVector> labels_shark(labels.size(), shark::RealVector(1));
	for (size_t ind = 0; ind < labels_shark.size(); ++ind)
	{
		labels_shark[ind][0] = labels[ind];
	}
	mgua_trainer.SetData(m_md_matrix_double, labels_shark, false);
	mgua_trainer.SetParameters(20, 6);
	mgua_trainer.Train(svm_vect, m_mol_folder + m_mol_prefix, false);
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