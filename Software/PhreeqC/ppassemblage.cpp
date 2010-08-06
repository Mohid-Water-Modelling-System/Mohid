#include "datastructures.h"

PPAssemblage::PPAssemblage(GlobalData *gd)
{
	this->gd = gd;

	eos_list = NULL;
	pure_phases = NULL;

  eos_list = new List<class ElementOfSpecies>;
  pure_phases = new List<class PurePhase>;

	Reset();
}

PPAssemblage::PPAssemblage(PPAssemblage *copy)
{
	eos_list = NULL;
	pure_phases = NULL;

  eos_list = new List<class ElementOfSpecies>;
  pure_phases = new List<class PurePhase>;

	copy->CopyTo(this);
}

PPAssemblage::~PPAssemblage()
{
	delete eos_list;
	delete pure_phases;
}

void PPAssemblage::Reset()
{
	eos_list->Clear();
	pure_phases->Clear();
}

void PPAssemblage::CopyTo(PPAssemblage *copy)
{
	eos_list->CopyTo(copy->eos_list);
	pure_phases->CopyTo(copy->pure_phases);
	copy->gd = gd;
}

int PPAssemblage::AddPurePhase(PPData &data, int &PurePhaseID)
{
	int index;
	PurePhase * pp = pure_phases->AddNew(index, &data.name);

	pp->name = data.name;
	pp->add_formula = data.add_formula;
	pp->force_equality = data.force_equality;
	pp->dissolve_only = data.dissolve_only;
	pp->moles = data.moles;
	pp->si = data.si;

	PurePhaseID = -1; //Not used for now
	//gd->phase_list.Search(&pp->name, PurePhaseID, true);

	return index;
}

void PPAssemblage::ChangeMoles(String *name, LDBLE new_moles)
{
	PurePhase * pp = pure_phases->Search(name);

	if (pp != NULL)
		pp->moles = new_moles;
}

void PPAssemblage::Print(FILE *file)
{
	fprintf(file, "Pure Phases info:\n\n");
	fprintf(file, "Phases:\n");

	PurePhase *p;
	for (int i = 0; i < pure_phases->Count(); i++)
	{
		p = pure_phases->Element(i);

		fprintf(file, "(%d) Name: %s\n", i, p->name.CharPtr());
		fprintf(file, "(%d) Add Formula: %s\n", i, p->add_formula.CharPtr());
		fprintf(file, "(%d) Moles: %.20e\n", i, p->moles);
		fprintf(file, "(%d) Delta: %.20e\n", i, p->delta);
		fprintf(file, "(%d) SI: %.20e\n\n", i, p->si);
	}
}