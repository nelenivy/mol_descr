#pragma once

namespace molecule_descriptor
{
struct Extensions
{
	static const char* Surface() { return ".wrl";}
	static const char* SingPts() { return ".sng";}
	static const char* NonMarkedSingPts() {return ".sngnm";}
	static const char* MeshTriangles() { return ".tr";}
	static const char* Charges() { return ".ch";}
	static const char* SegmSurf() { return ".sgm";}
	static const char* SurfWithTypes() { return ".tps";}
	static const char* Pairs() { return ".prs";}
	static const char* Triples() { return ".trps";}
	static const char* PairsMDMatrix() { return ".pmd";}
	static const char* TriplesMDMatrix() { return ".tmd";}
	static const char* PairKernelMatrix() { return ".pkr";}
	static const char* TriplesKernelMatrix() { return ".tkr";}
	static const char* Labels() { return ".lbs";}
	static const char* CountedLabels() { return ".clbs";}
	static const char* Quality() { return ".qlt";}
};

}