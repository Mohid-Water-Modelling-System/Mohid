#include "datastructures.h"

SpeciesInfo::SpeciesInfo()
{
	Reset();
}

SpeciesInfo::SpeciesInfo(SpeciesInfo *copy)
{
	copy->CopyTo(this);
}

SpeciesInfo::~SpeciesInfo()
{
}

void SpeciesInfo::Reset()
{
	master_s = NULL;
	s = NULL;

	coef = 0.0;
}

/*
SpeciesInfo *SpeciesInfo::Copy()
{
	return new SpeciesInfo(this);
}
*/

void SpeciesInfo::CopyTo(SpeciesInfo *copy)
{
	copy->master_s = master_s;
	copy->s = s;

	copy->coef = coef;
}

