#include "datastructures.h"

Species::Species()
{
	eos_list = NULL;
  e_sec_list = NULL;
	diff_layer = NULL;
  add_logk = NULL;
	rxn = NULL;
	rxn_s = NULL;
	rxn_x = NULL;
	eos_sys_total = NULL;

	primary = NULL;
	secondary = NULL;

	rxn = new Reaction;
	rxn_s = new Reaction;
	rxn_x = new Reaction;

	eos_list = new List<ElementOfSpecies>;
	e_sec_list = new List<ElementOfSpecies>;
	add_logk = new List<NameCoef>;
	diff_layer = new List<SpeciesDiffLayer>;
	eos_sys_total = new List<ElementOfSpecies>;

	Reset();
}

Species::Species(Species *copy)
{
	eos_list = NULL;
  e_sec_list = NULL;
	diff_layer = NULL;
  add_logk = NULL;
	rxn = NULL;
	rxn_s = NULL;
	rxn_x = NULL;
	eos_sys_total = NULL;

	primary = NULL;
	secondary = NULL;

	rxn = new Reaction;
	rxn_s = new Reaction;
	rxn_x = new Reaction;

	eos_list = new List<ElementOfSpecies>;
	e_sec_list = new List<ElementOfSpecies>;
	add_logk = new List<NameCoef>;
	diff_layer = new List<SpeciesDiffLayer>;
	eos_sys_total = new List<ElementOfSpecies>;

	copy->CopyTo(this);
}

Species::~Species()
{
	delete add_logk;
	delete e_sec_list;
	delete eos_list;
	delete rxn;
	delete rxn_s;
	delete rxn_x;
	delete diff_layer;
	delete eos_sys_total;
}

void Species::Reset(void)
{
  int i;

	name = "";
  mole_balance = "";

	if (rxn != NULL) rxn->Reset();
	if (rxn_s != NULL) rxn_s->Reset();
	if (rxn_x != NULL) rxn_x->Reset();

	if (add_logk != NULL) add_logk->Clear();
	if (eos_list != NULL) eos_list->Clear(); 
 	if (e_sec_list != NULL) e_sec_list->Clear();
	if (diff_layer != NULL) diff_layer->Clear();
	if (eos_sys_total != NULL) eos_sys_total->Clear();

	primary   = NULL;
	secondary = NULL;

  gfw    = 0.0;
  z      = 0.0;
  dw     = 0.0;
  alk    = 0.0;
  carbon = 0.0;
  co2    = 0.0;
  h      = 0.0;
  o      = 0.0;
  dha    = 0.0;
  dhb    = 0.0;
  lk     = 0.0;
  lg     = 0.0;
  lm     = 0.0;
  la     = 0.0;
  dg     = 0.0;
  moles  = 0.0;
	equiv  = 0.0;
	lg_pitzer = 0.0;
	dg_total_g = 0.0;
	tot_g_moles = 0.0;
	tot_dh2o_moles = 0.0;

	in             = false;
	check_equation = true;

  for (i = 0; i < 8; i++) logk[i] = 0.0;
  for (i = 0; i < 5; i++) cd_music[i] = 0.0;
  for (i = 0; i < 3; i++) dz[i] = 0.0;

	original_units = _kjoules_;

	number = -1;
	type  = 0;
  gflag = 0;
	exch_gflag = 0;
}

void Species::Clear(void)
{
  mole_balance = "";

	rxn_x->Reset();
	diff_layer->Clear();

	in = false;

	moles = 0.0;
}

/*
Species *Species::Copy()
{
	return new Species(*this);
}
*/

void Species::CopyTo(Species *copy)
{
	int i;

	eos_list->CopyTo(copy->eos_list);
	e_sec_list->CopyTo(copy->e_sec_list);
	add_logk->CopyTo(copy->add_logk);
	diff_layer->CopyTo(copy->diff_layer);
	eos_sys_total->CopyTo(copy->eos_sys_total);

	copy->primary = primary;
	copy->secondary = secondary;

	rxn->CopyTo(copy->rxn);
	rxn_s->CopyTo(copy->rxn_s);
	rxn_x->CopyTo(copy->rxn_x);

	copy->name = name;
	copy->mole_balance = mole_balance;

  copy->gfw = gfw;
  copy->z = z;
  copy->dw = dw;
  copy->alk = alk;
  copy->carbon = carbon;
  copy->co2 = co2;
  copy->h = h;
  copy->o = o;
  copy->dha = dha;
  copy->dhb = dhb;
  copy->lk = lk;
  copy->lg = lg;
  copy->lm = lm;
  copy->la = la;
  copy->dg = dg;
  copy->moles = moles;
	copy->equiv = equiv;
	copy->lg_pitzer = lg_pitzer;
	copy->dg_total_g = dg_total_g;
	copy->tot_g_moles = tot_g_moles;
	copy->tot_dh2o_moles = tot_dh2o_moles;


	for (i = 0; i < 8; i++) 
		copy->logk[i] = logk[i];
  for (i = 0; i < 5; i++) 
		copy->cd_music[i] = cd_music[i];
  for (i = 0; i < 3; i++) 
		copy->dz[i] = dz[i];

	copy->original_units = original_units;

	copy->type = type;
  copy->gflag = gflag;
	copy->number = number;
	copy->exch_gflag = exch_gflag;

  copy->check_equation = check_equation;
	copy->in = in;
}

void Species::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");
	int i;

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());

	fprintf(file, "%s  number: %d\n", spaces_str.CharPtr(), number);
	fprintf(file, "%s  type: %d\n", spaces_str.CharPtr(), type);
	fprintf(file, "%s  gflag: %d\n", spaces_str.CharPtr(), gflag);
	fprintf(file, "%s  exch_gflag: %d\n", spaces_str.CharPtr(), exch_gflag);
	fprintf(file, "%s  original_units: %d\n", spaces_str.CharPtr(), original_units);
	
	if (in)
		fprintf(file, "%s  in: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  in: FALSE\n", spaces_str.CharPtr());

	if (check_equation)
		fprintf(file, "%s  check_equation: TRUE\n", spaces_str.CharPtr());
	else
		fprintf(file, "%s  check_equation: FALSE\n", spaces_str.CharPtr());

	fprintf(file, "%s  gfw: %12.3e\n", spaces_str.CharPtr(), gfw);
	fprintf(file, "%s  z: %12.3e\n", spaces_str.CharPtr(), z);
	fprintf(file, "%s  dw: %12.3e\n", spaces_str.CharPtr(), dw);
	fprintf(file, "%s  equiv: %12.3e\n", spaces_str.CharPtr(), equiv);
	fprintf(file, "%s  alk: %12.3e\n", spaces_str.CharPtr(), alk);
	fprintf(file, "%s  co2: %12.3e\n", spaces_str.CharPtr(), co2);
	fprintf(file, "%s  carbon: %12.3e\n", spaces_str.CharPtr(), carbon);
	fprintf(file, "%s  h: %12.3e\n", spaces_str.CharPtr(), h);
	fprintf(file, "%s  o: %12.3e\n", spaces_str.CharPtr(), o);
	fprintf(file, "%s  dha: %12.3e\n", spaces_str.CharPtr(), dha);
	fprintf(file, "%s  dhb: %12.3e\n", spaces_str.CharPtr(), dhb);
	fprintf(file, "%s  lk: %12.3e\n", spaces_str.CharPtr(), lk);
	fprintf(file, "%s  lg: %12.3e\n", spaces_str.CharPtr(), lg);
	fprintf(file, "%s  lg_pitzer: %12.3e\n", spaces_str.CharPtr(), lg_pitzer);
	fprintf(file, "%s  lm: %12.3e\n", spaces_str.CharPtr(), lm);
	fprintf(file, "%s  la: %12.3e\n", spaces_str.CharPtr(), la);
	fprintf(file, "%s  dg: %12.3e\n", spaces_str.CharPtr(), dg);
	fprintf(file, "%s  dg_total_g: %12.3e\n", spaces_str.CharPtr(), dg_total_g);
	fprintf(file, "%s  moles: %12.3e\n", spaces_str.CharPtr(), moles);
	fprintf(file, "%s  tot_g_moles: %12.3e\n", spaces_str.CharPtr(), tot_g_moles);
	fprintf(file, "%s  tot_dh2o_moles: %12.3e\n", spaces_str.CharPtr(), tot_dh2o_moles);
	fprintf(file, "%s  logk: ", spaces_str.CharPtr());
	for (i = 0; i < 8; i++)
		fprintf(file, "%d:(%12.3e)  ",  i, logk[i]);  
	fprintf(file, "\n");
	fprintf(file, "%s  cd_music: ", spaces_str.CharPtr());
	for (i = 0; i < 5; i++)
		fprintf(file, "%d:(%12.3e)  ",  i, cd_music[i]);  
	fprintf(file, "\n");
	fprintf(file, "%s  dz: ", spaces_str.CharPtr());
	for (i = 0; i < 3; i++)
		fprintf(file, "%d:(%12.3e)  ",  i, dz[i]);  
	fprintf(file, "\n");

	if (primary != NULL)
		fprintf(file, "%s  primary: %s (%p)\n", spaces_str.CharPtr(), primary->name.CharPtr(), primary);
	else
		fprintf(file, "%s  primary: NULL\n", spaces_str.CharPtr());  
	if (secondary != NULL)
		fprintf(file, "%s  secondary: %s (%p)\n", spaces_str.CharPtr(), secondary->name.CharPtr(), secondary);
	else
		fprintf(file, "%s  secondary: NULL\n", spaces_str.CharPtr());  

	if(diff_layer != NULL && diff_layer->Count() > 0)
	{
		fprintf(file, "%s  DIFF_LAYER: (%p)\n", spaces_str.CharPtr(), diff_layer);
		for (i = 0; i < diff_layer->Count(); i++)
			(*diff_layer)[i]->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  DIFF_LAYER: NULL\n", spaces_str.CharPtr());

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