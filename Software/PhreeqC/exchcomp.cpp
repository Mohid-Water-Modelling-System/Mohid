#include "datastructures.h"

//===================================================================================
// Contructs
//===================================================================================
ExchComp::ExchComp()
{
	totals = NULL;
	formula_totals = NULL;

	totals = new List<ElementOfSpecies>; //ToDo: Check if delete or no the items of list
	formula_totals = new List<ElementOfSpecies>; //ToDo: Check if delete or no the items of list

	Reset();
}
//-----------------------------------------------------------------------------------
ExchComp::ExchComp(ExchComp *copy)
{
	totals = NULL;
	formula_totals = NULL;

	totals = new List<ElementOfSpecies>; //ToDo: Check if delete or no the items of list
	formula_totals = new List<ElementOfSpecies>; //ToDo: Check if delete or no the items of list

	copy->CopyTo(this);
}
//-----------------------------------------------------------------------------------
ExchComp::~ExchComp()
{
	delete totals;
	delete formula_totals;
}

//===================================================================================
// Interface
//===================================================================================
void ExchComp::Reset()
{
	formula_totals->Clear();
	totals->Clear();

	formula = "";
	phase_name = "";
	rate_name = "";

	formula_z = 0.0;
	moles = 0.0;
	la = 0.0;
	charge_balance = 0.0;
	phase_proportion = 0.0;

	master = NULL;
}

void ExchComp::CopyTo(ExchComp *copy)
{
	formula_totals->CopyTo(copy->formula_totals);
	totals->CopyTo(copy->totals);

	copy->formula = formula;
	copy->phase_name = phase_name;
	copy->rate_name = rate_name;

	copy->formula_z = formula_z;
	copy->moles = moles;
	copy->la = la;
	copy->charge_balance = charge_balance;
	copy->phase_proportion = phase_proportion;

	copy->master = master;
}

/*
ExchComp *ExchComp::Copy()
{
	return new ExchComp(this);
}
*/
//===================================================================================
// Protected
//===================================================================================


//===================================================================================
// Private
//===================================================================================
