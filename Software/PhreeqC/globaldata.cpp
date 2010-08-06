#include "datastructures.h"

GlobalData::GlobalData()
{
	s_hplus   = NULL;
	s_h3oplus = NULL;
	s_eminus  = NULL;
	s_h2o     = NULL;
	s_o2      = NULL;
	e_h_one   = NULL;
}

GlobalData::~GlobalData()
{
}

void GlobalData::ResetToExecuteModel()
{
	int i;

	Master *m;
	for (i = 0; i < master_list.Count(); i++)
	{
		m = master_list[i];
		m->Clear();
	}

	Element *e;
	for (i = 0; i < element_list.Count(); i++)
	{
		e =	element_list[i];
		e->Clear();
	}

	Species *s;
	for (i = 0; i < species_list.Count(); i++)
	{
		s = species_list[i];
		s->Clear();
	}

	Phase *p;
	for (i = 0; i < phase_list.Count(); i++)
	{
		p = phase_list[i];
		p->Clear(); //ToDo: Changed from Clear
	}
}