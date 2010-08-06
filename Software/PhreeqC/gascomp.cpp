#include "datastructures.h"

GasComp::GasComp()
{
	Reset();
}

GasComp::GasComp(GasComp *copy)
{
	copy->CopyTo(this);
}

GasComp::~GasComp()
{
}

void GasComp::Clear()
{
  phase = NULL;

  name = "";

  p_read = (LDBLE)NAN;
  moles = (LDBLE)0.0;
  initial_moles = (LDBLE)0.0;
}

void GasComp::Reset()
{
  phase = NULL;

  name = "";

  p_read = 0.0;
  moles = 0.0;
  initial_moles = 0.0;
}

/*
GasComp *GasComp::Copy()
{
	return new GasComp(this);
}
*/

void GasComp::CopyTo(GasComp *copy)
{
	copy->phase = phase;

	copy->name = name;

  copy->p_read = p_read;
  copy->moles = moles;
  copy->initial_moles = initial_moles;
}

void GasComp::PrintToFile(FILE *file, int spaces)
{
}

