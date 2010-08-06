#include "datastructures.h"

Phase::Phase()
{
  add_logk = NULL;
	eos_list = NULL;
	rxn = NULL;
	rxn_s = NULL;
	rxn_x = NULL;
	eos_sys_total = NULL;

	rxn = new Reaction;
	rxn_s = new Reaction;
	rxn_x = new Reaction;

	add_logk = new List<NameCoef>;
	eos_list = new List<ElementOfSpecies>;
	eos_sys_total = new List<ElementOfSpecies>;

	Reset();
}

Phase::Phase(Phase *copy)
{
  add_logk = NULL;
	eos_list = NULL;
	rxn = NULL;
	rxn_s = NULL;
	rxn_x = NULL;
	eos_sys_total = NULL;

	rxn = new Reaction;
	rxn_s = new Reaction;
	rxn_x = new Reaction;

	add_logk = new List<NameCoef>;
	eos_list = new List<ElementOfSpecies>;
	eos_sys_total = new List<ElementOfSpecies>;

	copy->CopyTo(this);
}

Phase::~Phase()
{
	delete add_logk;
	delete eos_list;
	delete rxn;
	delete rxn_s;
	delete rxn_x;
	delete eos_sys_total;
}

void Phase::Reset()
{
	name = "";
	formula = "";

	check_equation = true;
	in = false;
	in_system = false;

	original_units = _kjoules_;

	type = 0;

	rxn->Reset();
	rxn_s->Reset();
	rxn_x->Reset();

	add_logk->Clear();
	eos_list->Clear();
	eos_sys_total->Clear();

  for (int i = 0; i < 8; i++)
    logk[i] = 0.0;

	lk = 0.0;
	moles_x = 0.0; 
	p_soln_x = 0.0;
	fraction_x = 0.0;
	log10_lambda = 0.0; 
	log10_fraction_x = 0.0;
	dn = 0.0; 
	dnb = 0.0; 
	dnc = 0.0;
	gn = 0.0; 
	gntot = 0.0;
	gn_n = 0.0; 
	gntot_n = 0.0;
}

void Phase::Clear()
{
	in = false;
	in_system = false;

	if (rxn_x != NULL) rxn_x->Reset();
	if (eos_sys_total != NULL) eos_sys_total->Clear();

	lk = 0.0;
	moles_x = 0.0; 
	p_soln_x = 0.0;
	fraction_x = 0.0;
	log10_lambda = 0.0; 
	log10_fraction_x = 0.0;
	dn = 0.0; 
	dnb = 0.0; 
	dnc = 0.0;
	gn = 0.0; 
	gntot = 0.0;
	gn_n = 0.0; 
	gntot_n = 0.0;
}

/*
Phase *Phase::Copy()
{
	return new Phase(*this);
}
*/

void Phase::CopyTo(Phase *copy)
{
	int i;

	copy->name = name;
	copy->formula = formula;

	copy->check_equation = check_equation;
	copy->in = in;
	copy->in_system = in_system;

	copy->original_units = original_units;

	eos_list->CopyTo(copy->eos_list);
	add_logk->CopyTo(copy->add_logk);
	eos_sys_total->CopyTo(copy->eos_sys_total);

	rxn->CopyTo(copy->rxn);
	rxn_s->CopyTo(copy->rxn_s);
	rxn_x->CopyTo(copy->rxn_x);

  for (i = 0; i < 8; i++)
    copy->logk[i] = logk[i];

	copy->type = type;

	copy->lk = lk;
	copy->moles_x = moles_x; 
	copy->p_soln_x = p_soln_x;
	copy->fraction_x = fraction_x;
	copy->log10_lambda = log10_lambda; 
	copy->log10_fraction_x = log10_fraction_x;
	copy->dn = dn; 
	copy->dnb = dnb; 
	copy->dnc = dnc;
	copy->gn = gn; 
	copy->gntot = gntot;
	copy->gn_n = gn_n; 
	copy->gntot_n = gntot_n;
}

void Phase::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");
	int i;

	fprintf(file, "%s%s - formula: %s\n", spaces_str.CharPtr(), name.CharPtr(), formula.CharPtr());

	fprintf(file, "%s  type: %d\n", spaces_str.CharPtr(), type);
	fprintf(file, "%s  original_units: %d\n", spaces_str.CharPtr(), original_units);
	
	if (in)
		fprintf(file, "%s  in: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  in: FALSE\n", spaces_str.CharPtr());

	if (in_system)
		fprintf(file, "%s  in_system: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  in_system: FALSE\n", spaces_str.CharPtr());

	if (check_equation)
		fprintf(file, "%s  check_equation: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  check_equation: FALSE\n", spaces_str.CharPtr());

	fprintf(file, "%s  lk: %f\n", spaces_str.CharPtr(), lk);
	fprintf(file, "%s  moles_x: %f\n", spaces_str.CharPtr(), moles_x);
	fprintf(file, "%s  p_soln_x: %f\n", spaces_str.CharPtr(), p_soln_x);
	fprintf(file, "%s  fraction_x: %f\n", spaces_str.CharPtr(), fraction_x);
	fprintf(file, "%s  log10_lambda: %f\n", spaces_str.CharPtr(), log10_lambda);
	fprintf(file, "%s  log10_fraction_x: %f\n", spaces_str.CharPtr(), log10_fraction_x);
	fprintf(file, "%s  dn: %f\n", spaces_str.CharPtr(), dn);
	fprintf(file, "%s  dnb: %f\n", spaces_str.CharPtr(), dnb);
	fprintf(file, "%s  dnc: %f\n", spaces_str.CharPtr(), dnc);
	fprintf(file, "%s  gn: %f\n", spaces_str.CharPtr(), gn);
	fprintf(file, "%s  gntot: %f\n", spaces_str.CharPtr(), gntot);
	fprintf(file, "%s  gn_n: %f\n", spaces_str.CharPtr(), gn_n);
	fprintf(file, "%s  gntot_n: %f\n", spaces_str.CharPtr(), gntot_n);
	fprintf(file, "%s  logk: ", spaces_str.CharPtr());
	for (i = 0; i < 8; i++)
		fprintf(file, "%d:(%f)  ",  i, logk[i]);  
	fprintf(file, "\n");

	if(rxn != NULL && rxn->token_list->Count() > 0)
	{
		fprintf(file, "%s  RXN: (%p)\n", spaces_str.CharPtr(), rxn);
		rxn->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  RXN: NULL\n", spaces_str.CharPtr());
	if(rxn_s != NULL && rxn_s->token_list->Count() > 0)
	{
		fprintf(file, "%s  RXN_S: (%p)\n", spaces_str.CharPtr(), rxn_s);
		rxn_s->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  RXN_S: NULL\n", spaces_str.CharPtr());
	if(rxn_x != NULL && rxn_x->token_list->Count() > 0)
	{
		fprintf(file, "%s  RXN_X: (%p)\n", spaces_str.CharPtr(), rxn_x);
		rxn_x->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  RXN_X: NULL\n", spaces_str.CharPtr());

	if (eos_list != NULL && eos_list->Count() > 0)
	{
		fprintf(file, "%s  EOS_LIST: %d (%p)\n", spaces_str.CharPtr(), eos_list->Count(), eos_list);
		for(i = 0; i < eos_list->Count(); i++)
		{
			fprintf(file, "%s    %d:\n", spaces_str.CharPtr(), i);
			(*eos_list)[i]->PrintToFile(file, spaces + 5);
		}
	}
	else
		fprintf(file, "%s  EOS_LIST: NULL\n", spaces_str.CharPtr());
}