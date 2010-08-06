#include "datastructures.h"

Surface::Surface()
{
	comps = NULL;
	charge = NULL;

	comps = new List<SurfaceComp>;
	charge = new List<SurfaceCharge>;

	Reset();
}

Surface::Surface(Surface *copy)
{
	copy->CopyTo(this);
}

Surface::~Surface()
{
	delete comps;
	delete charge;
}

void Surface::Reset()
{
	comps->Clear();
	charge->Clear();

	dl_type = NO_DL;
	type = DDL;
	sites_units = SITES_ABSOLUTE;

	thickness = (LDBLE)1e-8;
	debye_lengths = 0;
	DDL_viscosity = 1.0;	
	DDL_limit = (LDBLE)0.8;	

	solution_equilibria = false;
	related_phases = false;
	related_rate = false;
	only_counter_ions = false;

	new_def = true;
	number = 0;
}

void Surface::CopyTo(Surface *copy)
{
	comps->CopyTo(copy->comps);
	charge->CopyTo(copy->charge);

	copy->dl_type = dl_type;
	copy->type = type;
	copy->sites_units = sites_units;

	copy->thickness = thickness;
	copy->debye_lengths = debye_lengths;
	copy->DDL_viscosity = DDL_viscosity;	
	copy->DDL_limit = DDL_limit;	

	copy->solution_equilibria = solution_equilibria;
	copy->related_phases = related_phases;
	copy->related_rate = related_rate;
	copy->solution_equilibria = solution_equilibria;

	copy->new_def = new_def;
	copy->number = number;
}

/*
Surface *Surface::Copy()
{
	return new Surface(this);
}
*/

void Surface::SetSurComp(SurCompData data)
{
	SurfaceComp *sc_p = comps->AddNew();

	sc_p->formula = data.formula;

	sc_p->use_moles = data.use_moles;
	sc_p->use_area = data.use_area;
	sc_p->use_grams = data.use_grams;

	sc_p->area = data.area;
	sc_p->grams = data.grams;
	sc_p->phase_name = data.phase_name;
	sc_p->rate_name = data.rate_name;
	sc_p->moles = data.moles;
	sc_p->phase_proportion = data.phase_proportion;
}
