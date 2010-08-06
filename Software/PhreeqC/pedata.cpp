#include "datastructures.h"

PEData::PEData()
{
	AllocMemory();
	Reset();
}

PEData::PEData(PEData *copy)
{
	AllocMemory();
	copy->CopyTo(this);
}

PEData::~PEData()
{
	delete rxn;
}

bool PEData::operator == (const PEData &right)
{
	return (element_1 == right.element_1 && element_2 == right.element_2 && valence_1 == right.valence_1 && valence_2 == right.valence_2);
}

void PEData::Reset()
{
	name = "pe";

	element_1 = "";
	element_2 = "";

	valence_1 = 0.0;
	valence_2 = 0.0;

	rxn->Reset();
}

/*
PEData *PEData::Copy()
{
	return new PEData(this);
}
*/

void PEData::CopyTo(PEData *copy)
{
	copy->name = name;

	copy->element_1 = element_1;
	copy->element_2 = element_2;

	copy->valence_1 = valence_1;
	copy->valence_2 = valence_2;

	rxn->CopyTo(copy->rxn);
}

void PEData::AllocMemory()
{
	rxn = NULL;
	rxn = new Reaction;
}