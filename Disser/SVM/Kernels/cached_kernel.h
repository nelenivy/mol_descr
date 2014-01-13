#pragma once
#include "shark/Models/Kernels/AbstractKernelFunction.h"
#include "gramm_matrix_cache.h"
//a wrap for shark kernels, which allow not recalculate kernel while cross-validation training
template <template <typename> class SharkKernelT, typename InputTypeT>
class CachedKernel
{
//private:
	//mutable GrammMatrixCache<SingularPointsSequence> m_mat_cache;
};