#pragma once
#include "CommonUtilities/attributes_container.h"

namespace molecule_descriptor
{
BOOL_PROP(SubgraphMask);
//for segment numbers
SIZE_T_PROP(SegmentNumProp);
//for marking visited
BOOL_PROP(VisitedProp);

PTR_TEMPLATE_PROP(PtrToSubgraph, GraphNode);
PTR_TEMPLATE_PROP(PtrToGraph, const GraphNode);

}