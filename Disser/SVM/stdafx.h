#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <fstream>
#include <array>

#include "opencv2/core/core.hpp"
#include "shark/Models/Kernels/KernelExpansion.h"
#include <shark/Core/IParameterizable.h>

#include "opencv2/core/core.hpp"
#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/Models/Kernels/KernelExpansion.h"
#include "shark/Models/Kernels/NormalizedKernel.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include <shark/Data/CVDatasetTools.h>
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>
#include "shark/Data/Csv.h"
#include "shark/Models/Kernels/AbstractKernelFunction.h"
#include "shark/Models/Kernels/LinearKernel.h"
#include "shark/Models/Kernels/NormalizedKernel.h"

#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/Algorithms/Trainers/EpsilonSvmTrainer.h"
#include "shark/ObjectiveFunctions/Loss/AbstractLoss.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include "shark/ObjectiveFunctions/Loss/AbsoluteLoss.h"
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>