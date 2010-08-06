#include "datastructures.h"

GasPhase::GasPhase()
{
	comps = NULL;
	comps = new List<class GasComp>;

	Clear();
}

GasPhase::GasPhase(class GasPhase * copy)
{
	comps = NULL;
	comps = new List<class GasComp>;

	copy->CopyTo(this);
}

GasPhase::~GasPhase()
{
	delete comps;
}

void GasPhase::Clear()
{
  solution_equilibria = false;

  type = PRESSURE;

  total_p = (LDBLE)1.0;
  total_moles = (LDBLE)0.0;
  volume = (LDBLE)1.0;
  temperature = (LDBLE)298.15;

  comps->Clear();
}

void GasPhase::Reset()
{
  solution_equilibria = false;

  type = PRESSURE;

  total_p = (LDBLE)1.0;
  total_moles = (LDBLE)0.0;
  volume = (LDBLE)1.0;
  temperature = (LDBLE)298.15;

  comps->Clear();
}

/*
GasPhase *GasPhase::Copy()
{
	return new GasPhase(this);
}
*/

void GasPhase::CopyTo(GasPhase *copy)
{
  copy->solution_equilibria = solution_equilibria;
  copy->type = type;

  copy->total_p = total_p;
  copy->total_moles = total_moles;
  copy->volume = volume;
  copy->temperature = temperature;

	comps->CopyTo(copy->comps);
}

void GasPhase::PrintToFile(FILE *file, int spaces)
{
}

void GasPhase::AddGasComp(String name, LDBLE p_read)
{
	GasComp *new_item = comps->AddNew();

	new_item->name = name;
	new_item->p_read = p_read;

	comps->Sort(0, true); //sort ignoring case
}