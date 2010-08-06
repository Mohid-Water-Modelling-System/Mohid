#include "datastructures.h"

//===================================================================================
// Contructs
//===================================================================================
Exchange::Exchange()
{
	comps = NULL;

	comps = new List<ExchComp>; //ToDo: Check if delete or no the items of list

	Reset();
}
//-----------------------------------------------------------------------------------	
Exchange::Exchange(Exchange *copy)
{
	copy->CopyTo(this);
}
//-----------------------------------------------------------------------------------	
Exchange::~Exchange()
{
	delete comps;
}

//===================================================================================
// Interface
//===================================================================================
void Exchange::Reset()
{
	name = "";

	number = -1;

	comps->Clear();

	solution_equilibria = false;
	new_def = true;
	related_phases = false;
	related_rate = false;
}

void Exchange::CopyTo(Exchange *copy)
{
	copy->name = name;

	copy->number = number;

	comps->CopyTo(copy->comps);

	copy->solution_equilibria = solution_equilibria;
	copy->new_def = new_def;
	copy->related_phases = related_phases;
	copy->related_rate = related_rate;
}

/*
Exchange *Exchange::Copy()
{
	return new Exchange(this);
}
*/

int Exchange::AddExchanger(ExchCompData *exchanger)
{
	LDBLE conc;
	int index;

	ExchComp *new_comp = comps->AddNew(index, &exchanger->formula);
	new_comp->formula = exchanger->formula;
	new_comp->type = exchanger->type;
	new_comp->amount = exchanger->amount;

	if (exchanger->type == 0) //exchanger/concentration pair
	{
		conc = exchanger->amount;
	}
	else if (exchanger->type == 1) //exchanger/pure-phase pair
	{
		new_comp->phase_name = name;
		related_phases = true;
		new_comp->phase_proportion = exchanger->amount;
	}
	else //exchanger/kinetic pair
	{
		new_comp->rate_name = name;
		related_rate = true;	
		new_comp->phase_proportion = exchanger->amount;
	}

	return index;
}

void Exchange::Print(FILE *file)
{
	fprintf(file, "Exchange info:\n");

	fprintf(file, "Components:\n\n");

	ExchComp *c;
	for (int i = 0; i < comps->Count(); i++)
	{
		c = comps->Element(i);

		fprintf(file, "(%d) Formula: %s\n", i, c->formula.CharPtr());
		fprintf(file, "(%d) Amount: %.20e\n", i, c->amount);
		fprintf(file, "(%d) Type: %d\n", i, c->type);
		fprintf(file, "(%d) Name (phase or kinetic related): %s\n", i, c->phase_name.CharPtr());		
		fprintf(file, "(%d) Moles: %.20e\n\n", i, c->moles);
	}
}

//===================================================================================
// Protected
//===================================================================================


//===================================================================================
// Private
//===================================================================================
