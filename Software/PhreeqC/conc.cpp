#include "datastructures.h"

Conc::Conc()
{
	m_list = NULL;
	m_list = new ListOfPointers<Master>;
	
	Reset();
}

Conc::Conc(Conc *copy)
{
	m_list = NULL;
	m_list = new ListOfPointers<Master>;

	copy->CopyTo(this);
}

Conc::~Conc()
{
	delete m_list;
}

void Conc::Reset()
{
	type = _unknown_type_;

	charge = false;

	moles = 0.0;
	input_conc = 0.0;
	phase_si = 0.0;
	gfw = 0.0;

	name = "";
	equation_name = "";
	as = "";

	n_pe = 0;

	units = _unknown_unit_;

	m_list->Clear();

	//x = NULL;
}

/*
Conc *Conc::Copy()
{
	return new Conc(this);
}
*/

void Conc::CopyTo(Conc *copy)
{
	copy->name = name;

	copy->type = type;
	copy->charge = charge;

	copy->moles = moles;
	copy->input_conc = input_conc;
	copy->phase_si = phase_si;
	copy->gfw = gfw;

	copy->equation_name = equation_name;
	copy->as = as;

	copy->n_pe = n_pe;

	copy->units = units;

	//copy->x = x;

	m_list->CopyTo(copy->m_list);
}

