#include "Manager/set_manager.h"

int main(int argc, char** argv)
{
	molecule_descriptor::SetManager manager;
	manager.Init(argc, argv);
	manager.ProcessSingularPoints(true);
	//manager.ProcessPairs(true);
	manager.ProcessPairsLevels(true);
	//manager.ProcessKernelSVMPointsWithFiltering();
	//manager.ProcessTriples(true);
	//manager.ProcessMGUASVMRegression();
	manager.ProcessMGUASVMClassification(30, 15);
	//manager.ProcessDescriptorsSVM();
	//manager.ProcessKernelSVMPoints();
	//manager.ProcessKernelSVMHistograms();
}