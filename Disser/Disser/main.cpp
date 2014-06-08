#include "Manager/set_manager.h"

int main(int argc, char** argv)
{
	molecule_descriptor::SetManager manager;
	manager.Init(argc, argv);
	manager.ProcessSingularPoints(false);
	manager.ProcessPairs(true);
	//manager.ProcessKernelSVMPointsWithFiltering();
	//manager.ProcessTriples(true);
	//manager.ProcessMGUASVMRegression();
	manager.ProcessMGUASVMClassification(20, 10);
	//manager.ProcessDescriptorsSVM();
	//manager.ProcessKernelSVMPoints();
	//manager.ProcessKernelSVMHistograms();
}