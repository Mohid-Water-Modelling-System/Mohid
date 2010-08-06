#ifndef specieslistH
#define specieslistH

#include "datastructures.h"

class SpeciesList:
	public List<SpeciesInfo>
{
public:
	SpeciesList();
	~SpeciesList();

public:
	bool Sort(Species *s_hplus, int start = 0);
	bool SortByMasterOnly(Species *s_hplus, int start = 0);

private:
	int Compare(SpeciesInfo *a, SpeciesInfo *b);
	int CompareByMasterOnly(SpeciesInfo *a, SpeciesInfo *b);

private:
	Species *s_hplus;

};

#endif