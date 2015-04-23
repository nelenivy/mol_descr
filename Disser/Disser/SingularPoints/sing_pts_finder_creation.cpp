#include "i_singular_points_finder.h"
#include "singular_points_finder.h"
#include "sng_pts_finder_scale_space.h"

namespace molecule_descriptor
{

std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder(const SingularPointsAlgorithm alg)
{
	ISingularPointsFinder* instance = nullptr;
	if (alg == SEGMENTATION)
	{
		instance = new SngPtsFinderSegmentation();
	}
	else if (alg == SCALE_SPACE)
	{
		instance = new SngPtsFinderScaleSpace();
	}
	return std::shared_ptr<ISingularPointsFinder>(instance);
}

}