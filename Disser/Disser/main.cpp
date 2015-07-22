#include "Manager/set_manager.h"
#include "InputOutput/params_reader.h"

int main(int argc, char** argv)
{
	bool use_levels = false, 
		calculate_singular_points = true, 
		calculate_pairs = true,
		use_triples = false,
		calculate_descriptors = true;
	int mgua_iterations_num = 30, mgua_descriptors_num = 15;

	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-use_levels", use_levels, true);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-calculate_singular_points", calculate_singular_points, true);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-calculate_pairs", calculate_pairs, true);
	
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-mgua_iterations_num", mgua_iterations_num, 30);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-mgua_descriptors_num", mgua_descriptors_num, 15);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-use_triples", use_triples, false);
	molecule_descriptor::ReadParamFromCommandLineWithDefault(argc, argv, "-calculate_descriptors", calculate_descriptors, true);
	
	molecule_descriptor::SetManager manager;
	manager.Init(argc, argv);

	if (calculate_descriptors)
	{
		manager.ProcessSingularPoints(calculate_singular_points);

		if (!use_triples)
		{
			if (use_levels)
			{
				manager.ProcessPairsLevels(calculate_pairs);
			}
			else
			{
				manager.ProcessPairs(calculate_pairs);
			}
		}
		else
		{
			if (use_levels)
			{
				manager.ProcessTriplesLevels(calculate_pairs);
			}
			else
			{
				manager.ProcessTriples(calculate_pairs);
			}
		}
	}


	//manager.ProcessKernelSVMPointsWithFiltering();
	//manager.ProcessTriples(true);
	//manager.ProcessMGUASVMRegression();
	manager.ProcessSVMClassificationL0();
	//manager.ProcessMGUASVMClassification(mgua_iterations_num, mgua_descriptors_num);
	//manager.ProcessDescriptorsSVM();
	//manager.ProcessKernelSVMPoints();
	//manager.ProcessKernelSVMHistograms();
}