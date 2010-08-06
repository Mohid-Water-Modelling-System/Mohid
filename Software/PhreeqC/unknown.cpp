#include "datastructures.h"

Unknown::Unknown()
{
	master         = NULL;
	comp_unknowns  = NULL;

	master = new ListOfPointers<Master>;
	comp_unknowns = new ListOfPointers<Unknown>;

	Reset();
}

Unknown::Unknown(Unknown *copy)
{
	master         = NULL;
	comp_unknowns  = NULL;

	master = new ListOfPointers<Master>;
	comp_unknowns = new ListOfPointers<Unknown>;

	copy->CopyTo(this);
}

Unknown::~Unknown()
{
	delete master;
	delete comp_unknowns;
}

void Unknown::Reset()
{
	molality_result = 0.0;
	moles_result = 0.0;

	type = 0;

	moles = 0.0;
	ln_moles = 0.0;
	f = 0.0;
	sum = 0.0;
	delta = 0.0;
	la = 0.0;
	si = 0.0;
	related_moles = 0.0;
	mass_water = 0.0;

	dissolve_only = false;
	s_s_in = false;

	number = 0;
	s_s_comp_number = -1;

	name = "";

	master->Clear();
	comp_unknowns->Clear();

	total = NULL;
	p = NULL; 
	surface_comp = NULL;
	pure_phase = NULL;
	gas_phase = NULL;
	s = NULL;
	surface_charge = NULL;
	phase_unknown = NULL;
	s_s_comp = NULL;
}

/*
Unknown *Unknown::Copy()
{
	return new Unknown(this);
}
*/

void Unknown::CopyTo(Unknown *copy)
{
	copy->type = type;

	copy->moles = moles;
	copy->ln_moles = ln_moles;
	copy->f = f;
	copy->sum = sum;
	copy->delta = delta;
	copy->la = la;
	copy->si = si;
	copy->related_moles = related_moles;
	copy->mass_water = mass_water;
	copy->s_s_comp_number = s_s_comp_number;

	copy->dissolve_only = dissolve_only;
	copy->s_s_in = s_s_in;

	copy->number = number;

	copy->name = name;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// this piece of code garantees that the list get only a pointer to the element, not a COPY of it
	// this is necessary because "master" and "comp_unkowns" lists stores POINTERS to the elements, not COPIES of them
	master->CopyTo(copy->master);
	comp_unknowns->CopyTo(copy->comp_unknowns);
	
	/*
	copy->master->Clear();
	for (i = 0; i < master->Count(); i++)
		copy->master->Add((*master)[i]);

	copy->comp_unknowns->Clear();
	for (i = 0; i < comp_unknowns->Count(); i++)
		copy->comp_unknowns->Add((*comp_unknowns)[i]);
	*/
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	copy->total = total;
	copy->p = p; 
	copy->surface_comp = surface_comp;
	copy->pure_phase = pure_phase;
	copy->gas_phase = gas_phase;
	copy->s = s;
	copy->surface_charge = surface_charge;
	copy->phase_unknown = phase_unknown;
	copy->s_s_comp = s_s_comp;
}

void Unknown::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());

	fprintf(file, "%s  number: %d\n", spaces_str.CharPtr(), number);
	fprintf(file, "%s  type: %d\n", spaces_str.CharPtr(), type);
	
	if (s_s_in)
		fprintf(file, "%s  s_s_in: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  s_s_in: FALSE\n", spaces_str.CharPtr());

	if (dissolve_only)
		fprintf(file, "%s  dissolve_only: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  dissolve_only: FALSE\n", spaces_str.CharPtr());

	fprintf(file, "%s  moles: %f\n", spaces_str.CharPtr(), moles);
	fprintf(file, "%s  ln_moles: %f\n", spaces_str.CharPtr(), ln_moles);
	fprintf(file, "%s  f: %f\n", spaces_str.CharPtr(), f);
	fprintf(file, "%s  sum: %f\n", spaces_str.CharPtr(), sum);
	fprintf(file, "%s  delta: %f\n", spaces_str.CharPtr(), delta);
	fprintf(file, "%s  la: %f\n", spaces_str.CharPtr(), la);
	fprintf(file, "%s  si: %f\n", spaces_str.CharPtr(), si);
	fprintf(file, "%s  related_moles: %f\n", spaces_str.CharPtr(), related_moles);
	fprintf(file, "%s  mass_water: %f\n", spaces_str.CharPtr(), mass_water);
	
	if (total != NULL)
		fprintf(file, "%s  total: %s (%p)\n", spaces_str.CharPtr(), total->name.CharPtr(), total);
	else
		fprintf(file, "%s  total: NULL\n", spaces_str.CharPtr()); 

	if (p != NULL)
		fprintf(file, "%s  p: %s (%p)\n", spaces_str.CharPtr(), p->name.CharPtr(), p);
	else
		fprintf(file, "%s  p: NULL\n", spaces_str.CharPtr());  

	if (s != NULL)
		fprintf(file, "%s  s: %s (%p)\n", spaces_str.CharPtr(), s->name.CharPtr(), s);
	else
		fprintf(file, "%s  s: NULL\n", spaces_str.CharPtr());  
}