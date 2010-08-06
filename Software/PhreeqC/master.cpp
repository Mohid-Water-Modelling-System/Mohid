#include "datastructures.h"

Master::Master()
{
	e = NULL;
	s = NULL;
	u = NULL;

	pe_rxn = NULL;

	rxn_primary = NULL;
	rxn_secondary = NULL;

	rxn_primary = new Reaction;
	rxn_secondary = new Reaction;

	Reset();
}

Master::Master(Master *copy)
{
	e = NULL;
	s = NULL;
	u = NULL;

	pe_rxn = NULL;

	rxn_primary = NULL;
	rxn_secondary = NULL;

	rxn_primary = new Reaction;
	rxn_secondary = new Reaction;

	copy->CopyTo(this);
}

Master::~Master()
{
	delete rxn_primary;
	delete rxn_secondary;
}

void Master::Clear(void)
{
	u = NULL;
	
	pe_rxn = NULL;

	in = false;
	rewrite = false;

	total_primary = 0.0;
  total = 0.0;  
}

void Master::Reset(void)
{
	name = "";
  gfw_formula = "";

  e = NULL;
  s = NULL;
	u = NULL;
	
	rxn_primary->Reset();
	rxn_secondary->Reset();

	pe_rxn = NULL;

  primary = false;
	in = false;
	rewrite = false;

  number = -1;
  type = 0;

	total_primary = 0.0;
	alk = 0.0;
  gfw = 0.0;
  coef = 0.0;
  total = 0.0;  
}

/*
Master *Master::Copy()
{
	return new Master(*this);
}
*/

void Master::CopyTo(Master *copy)
{
	copy->name        = name;
  copy->gfw_formula = gfw_formula;

  copy->e = e;
  copy->s = s;
	copy->u = u;

	rxn_primary->CopyTo(copy->rxn_primary);
	rxn_secondary->CopyTo(copy->rxn_secondary);
	//*(copy->pe_rxn) = *(pe_rxn);
	copy->pe_rxn = pe_rxn;

  copy->primary = primary;
  copy->in      = in;
	copy->rewrite = rewrite;

	copy->number = number;
  copy->type   = type;

	copy->total_primary = total_primary;
	copy->alk           = alk;
  copy->gfw           = gfw;
  copy->coef          = coef;
  copy->total         = total;  
}

void Master::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());

	fprintf(file, "%s  number: %d\n", spaces_str.CharPtr(), number);
	fprintf(file, "%s  type: %d\n", spaces_str.CharPtr(), type);

	if (in)
		fprintf(file, "%s  IN: true\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  IN: false\n", spaces_str.CharPtr());
	if (rewrite)
		fprintf(file, "%s  REWRITE: true\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  REWRITE: false\n", spaces_str.CharPtr());
	if (primary)
		fprintf(file, "%s  PRIMARY: true\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  PRIMARY: false\n", spaces_str.CharPtr());

	fprintf(file, "%s  coef: %f\n", spaces_str.CharPtr(), coef);
	fprintf(file, "%s  total: %f\n", spaces_str.CharPtr(), total);
	fprintf(file, "%s  total_primary: %f\n", spaces_str.CharPtr(), total_primary);
	fprintf(file, "%s  alk: %f\n", spaces_str.CharPtr(), alk);
	fprintf(file, "%s  gfw: %f\n", spaces_str.CharPtr(), gfw);
	fprintf(file, "%s  gfw_formula: %s\n", spaces_str.CharPtr(), gfw_formula.CharPtr());

	if (e != NULL)
		fprintf(file, "%s  e: %s (%p)\n", spaces_str.CharPtr(), e->name.CharPtr(), e);
	else
		fprintf(file, "%s  e: NULL\n", spaces_str.CharPtr());  
	if (u != NULL && !u->name.IsEmpty())
		fprintf(file, "%s  u: %s (%p)\n", spaces_str.CharPtr(), u->name.CharPtr(), u);
	else if (u != NULL)
		fprintf(file, "%s  u: (%p)\n", spaces_str.CharPtr(), u->name.CharPtr(), u);
	else
		fprintf(file, "%s  u: NULL\n", spaces_str.CharPtr());  
	if (s != NULL)
		fprintf(file, "%s  s: %s (%p)\n", spaces_str.CharPtr(), s->name.CharPtr(), s);
	else
		fprintf(file, "%s  s: NULL\n", spaces_str.CharPtr()); 

	if(rxn_primary != NULL)
	{
		fprintf(file, "%s  RXN_PRIMARY: (%p)\n", spaces_str.CharPtr(), rxn_primary);
		rxn_primary->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  RXN_PRIMARY: NULL\n", spaces_str.CharPtr());
	if(rxn_secondary != NULL)
	{
		fprintf(file, "%s  RXN_SECONDARY: (%p)\n", spaces_str.CharPtr(), rxn_secondary);
		rxn_secondary->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  RXN_SECONDARY: NULL\n", spaces_str.CharPtr());
	if(pe_rxn != NULL && *pe_rxn)
	{
		fprintf(file, "%s  PE_RXN: (%p)\n", spaces_str.CharPtr(), *pe_rxn);
		(*pe_rxn)->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  PE_RXN: NULL\n", spaces_str.CharPtr());
	fprintf(file, "=====================================\n");
}