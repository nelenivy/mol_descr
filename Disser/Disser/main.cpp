#include "Manager/set_manager.h"
#include "InputOutput/params_reader.h"

int main(int argc, char** argv)
{
	bool use_levels = false, calculate_singular_points = true, calculate_pairs = true;
	int mgua_iterations_num = 30, mgua_descriptors_num = 15;

	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-use_levels", use_levels, true);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-calculate_singular_points", calculate_singular_points, true);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-calculate_pairs", calculate_pairs, true);
	
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-mgua_iterations_num", mgua_iterations_num, 30);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-mgua_descriptors_num", mgua_descriptors_num, 15);
	molecule_descriptor::SetManager manager;
	manager.Init(argc, argv);
	manager.ProcessSingularPoints(calculate_singular_points);

	if (use_levels)
	{
		manager.ProcessPairsLevels(calculate_pairs);
	}
	else
	{
		manager.ProcessPairs(calculate_pairs);
	}

	//manager.ProcessKernelSVMPointsWithFiltering();
	//manager.ProcessTriples(true);
	//manager.ProcessMGUASVMRegression();
	manager.ProcessMGUASVMClassification(mgua_iterations_num, mgua_descriptors_num);
	//manager.ProcessDescriptorsSVM();
	//manager.ProcessKernelSVMPoints();
	//manager.ProcessKernelSVMHistograms();
}