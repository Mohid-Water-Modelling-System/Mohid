#include "datastructures.h"

SurfaceComp::SurfaceComp()
{
	formula_totals = NULL;
	totals = NULL;

	formula_totals = new List<ElementOfSpecies>;
	totals = new List<ElementOfSpecies>;

	Reset();
}

SurfaceComp::SurfaceComp(SurfaceComp *copy)
{
	formula_totals = NULL;
	totals = NULL;

	formula_totals = new List<ElementOfSpecies>;
	totals = new List<ElementOfSpecies>;

	copy->CopyTo(this);
}

SurfaceComp::~SurfaceComp()
{
	delete formula_totals;
	delete totals;
}

void SurfaceComp::Reset()
{
	formula_totals->Clear();
	totals->Clear();

	master = NULL;

	name = "";
	formula = "";
	phase_name = "";
	rate_name = "";

	formula_z = 0.0;
	moles = 0;
	la = 0;
	charge = 0;
	cb = 0;
	phase_proportion = 0;
	Dw = 0;
}

void SurfaceComp::CopyTo(SurfaceComp *copy)
{
	formula_totals->CopyTo(copy->formula_totals);
	totals->CopyTo(copy->totals);

	copy->name = name;
	copy->formula = formula;
	copy->phase_name = phase_name;
	copy->rate_name = rate_name;

	copy->formula_z = formula_z;
	copy->moles = moles;
	copy->la = la;
	copy->charge = charge;
	copy->cb = cb;
	copy->phase_proportion = phase_proportion;
	copy->Dw = Dw;
}

/*
SurfaceComp *SurfaceComp::Copy()
{
	return new SurfaceComp(this);
}
*/
