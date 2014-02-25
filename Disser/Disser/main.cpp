#include "Manager/set_manager.h"

int main(int argc, char** argv)
{
	molecule_descriptor::SetManager manager;
	manager.Init(argc, argv);
	manager.ProcessSingularPoints(false);
	//manager.ProcessPairs(false);
	//manager.ProcessTriples(true);
	//manager.ProcessMGUASVM();
	//manager.ProcessDescriptorsSVM();
	//manager.ProcessKernelSVMPoints();
	manager.ProcessKernelSVMHistograms();
}