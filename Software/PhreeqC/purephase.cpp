#include "datastructures.h"

PurePhase::PurePhase()
{
	Reset();
}

PurePhase::PurePhase(PurePhase *copy)
{
	copy->CopyTo(this);
}

PurePhase::~PurePhase()
{
}

void PurePhase::Reset()
{
	phase = NULL;

	name = "";
	add_formula = "";

	si = 0.0;
  moles = 0.0;
  delta = 0.0;
  initial_moles = 0.0;

	force_equality = false;
	dissolve_only = false;
}

void PurePhase::Clear()
{
	phase = NULL;

	name = "";
	add_formula = "";

	si = 0.0;
  moles = 0.0;
  delta = 0.0;
  initial_moles = 0.0;

	force_equality = false;
	dissolve_only = false;
}

/*
PurePhase *PurePhase::Copy()
{
	return new PurePhase(this);
}
*/

void PurePhase::CopyTo(PurePhase *copy)
{
	copy->phase = phase;

	copy->name = name;
	copy->add_formula = add_formula;

	copy->si = si;
  copy->moles = moles;
  copy->delta = delta;
  copy->initial_moles = initial_moles;

	copy->force_equality = force_equality;
	copy->dissolve_only = dissolve_only;
}

