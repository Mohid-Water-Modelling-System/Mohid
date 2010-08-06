#include <math.h>
#include <float.h>
#include "modelengine.h"
#include "conio.h"

//===========================
// Constructors & Destructors
//===========================
ModelEngine::ModelEngine(GlobalData *gd, ModelData *md):
SSAssemblageModel(gd, md)
{
	this->gd = gd;
	this->md = md;

	SetData(md);

	r_temp = NULL;
	r_temp = new Reaction;

	count_unknowns = 0;

	error = NULL;
	error_string[0] = '\0';
}

ModelEngine::~ModelEngine()
{
	delete r_temp;
}

//===========================
// Public functions
//===========================
bool ModelEngine::ConvertUnits(Solution *sol_p)
{
  int i;
  LDBLE sum_solutes;
  Master *m;
  Conc *conc;
	
  sum_solutes = exp(-sol_p->ph * md->LOG_10);

  for (i = 0; i < sol_p->totals->Count(); i++)
	{
    conc = (*sol_p->totals)[i];

    conc->moles = 0.0;
    if (conc->name == "H(1)" || conc->name == "E")
      continue;

		if (conc->input_conc <= 0)
      continue;

    if (conc->gfw <= 0.0)
    {
			if (!conc->as.IsEmpty())
      {
				if (!ComputeGFW(conc->as, conc->gfw))
					return false;

				if (conc->name == "Alkalinity" && conc->as == "CaCO3")
					conc->gfw /= 2.0;
      }
      else
      {
				if ((m = gd->master_list.Search(&conc->name(0, " "), true)) != NULL)
					conc->gfw = m->gfw;
				else
					return false;
      }
    }

    conc->moles = conc->input_conc;

    if (Check_l_(sol_p->units))
      conc->moles *= (LDBLE)1.0 / (sol_p->density);

		conc->moles *= ConversionFactorForMoles(conc->units);

    if (Check_g_kgs_(conc->units) || Check_g_l_(conc->units))
      sum_solutes += conc->moles;
    else if (Check_mol_kgs_(conc->units) || Check_mol_l_(conc->units) || Check_eq_l_(conc->units))
      sum_solutes += (conc->moles * conc->gfw);

		if (Check_g_(conc->units) && conc->gfw != 0.0)
      conc->moles /= conc->gfw;
  }

	if (Check_kgs_(sol_p->units) || Check_l_(sol_p->units))
  {
    md->mass_water_aq_x = (LDBLE)1.0 - ((LDBLE)1e-3 * sum_solutes);
    for (i = 0; i < sol_p->totals->Count(); i++)
      (*sol_p->totals)[i]->moles /= md->mass_water_aq_x;
  }

  md->mass_water_aq_x = sol_p->mass_water;
  for (i = 0; i < sol_p->totals->Count(); i++)
    (*sol_p->totals)[i]->moles *= md->mass_water_aq_x;

	return true;
}


bool ModelEngine::Clear()
{
	int i;
	Solution *sol_p = md->use.sol_p;

	// Clear species solution-dependent data
	for (i = 0; i < gd->species_list.Count(); i++)
	{
    gd->species_list[i]->in = false;
		//gd->species_list[i]->rewrite = false;
	}

	// Set pe structure
	sol_p->pe->CopyTo(md->pe_x);
  md->default_pe_x = sol_p->default_pe;

	// Clear master species solution-dependent data
	Master *m;
	PEData *ped;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];
    
		m->in = false;
		m->rewrite = false;
		m->u = NULL;

		ped = (*md->pe_x)[sol_p->default_pe];
		m->pe_rxn = &(ped->rxn);

		// copy primary reaction to secondary reaction
		m->rxn_primary->CopyTo(m->rxn_secondary);
  }

  if (md->state == INITIAL_SOLUTION)
  {
    gd->s_h2o->secondary->in = true;
		gd->s_h2o->secondary->rewrite = false;
    gd->s_hplus->secondary->in = true;
		gd->s_hplus->secondary->rewrite = false;
  }
  else
  {
    gd->s_h2o->primary->in = true;
		gd->s_h2o->primary->rewrite = false;
    gd->s_hplus->primary->in = true;
		gd->s_hplus->primary->rewrite = false;
  }
  gd->s_eminus->primary->in = true;
	gd->s_eminus->primary->rewrite = false;

	// Set all unknown pointers to NULL
  md->mb_unknown = NULL;
  md->ah2o_unknown = NULL;
  md->mass_hydrogen_unknown = NULL;
  md->mass_oxygen_unknown = NULL;
  md->mu_unknown = NULL;
  md->alkalinity_unknown = NULL;
  md->carbon_unknown = NULL;
  md->ph_unknown = NULL;
  md->pe_unknown = NULL;
  md->charge_balance_unknown = NULL;
  md->solution_phase_boundary_unknown = NULL;
  md->pure_phase_unknown = NULL;
  md->exchange_unknown = NULL;
  md->surface_unknown = NULL;
  md->gas_unknown = NULL;
  md->s_s_unknown = NULL;
 
	// Free arrays used in model   
  FreeModelAllocs(); //free_model_allocs ();

  return true;
}


//===========================
// Private functions
//===========================
bool ModelEngine::GetMasterList(String &list, Master *m, ListOfPointers<Master> *copy)
{
	String *token;	
	Master *m_ptr, *m_1;
	int i;

	tl.SetString(list);
	copy->Clear();

	// Make list of master species pointers
	if (m == m->s->primary)
	{
		// First in list is primary species
    for (i = 0; i < gd->master_list.Count(); i++)
		{
			m_1 = gd->master_list[i];
			if (m_1 == m)
				break;
		}
    i++;		
		
		if (i >= gd->master_list.Count())
		{
			// Element has only one valence
			copy->AddNew(m);
		}
		else
		{
			m_1 = gd->master_list[i];

			if (m != m_1->e->primary)
			{
				// Element has only one valence
				copy->AddNew(m);
			}
			else
			{
				// Element has multiple valences
				if (m->s->secondary == NULL)
					return false;

				copy->AddNew(m->s->secondary);

				while (i < gd->master_list.Count() && (m_ptr = gd->master_list[i])->e->primary == m)
				{
					if (m_ptr->s->primary == NULL)
						copy->AddNew(m_ptr);

					i++;
				}
			}
		}
	}
	else
	{
		// First in list is secondary species, Include all valences from input
    copy->AddNew(m);

		while (!tl.EOTL())
    {
			token = tl.NextToken();
			m_ptr = gd->master_list.Search(token, true);

      if (m_ptr != NULL)
				copy->AddNew(m_ptr);
    }
  }

	return true;
}

bool ModelEngine::RewriteMasterToSecondary(Master *m1, Master *m2)
{
	//
	// Write equation for secondary master species in terms of another secondary master species
	// Store result in rxn_secondary of master_ptr.
	//

  LDBLE coef1, coef2;
  Master *m_p1, *m_p2;

	// Check that the two master species have the same primary master species
	m_p1 = m1->e->primary;
  m_p2 = m2->e->primary;

  if (m_p1 != m_p2 || m_p1 == NULL)
		return false;

	// Find coefficient of primary master in reaction
	if (!GetReactionCoef(m1->rxn_primary, m_p1->s->name, coef1) || !GetReactionCoef(m2->rxn_primary, m_p1->s->name, coef2))
		return false;

  if (Equal(coef1, (LDBLE)0.0, (LDBLE)TOLERANCE) || Equal(coef2, (LDBLE)0.0, (LDBLE)TOLERANCE))
		return false;

	// Rewrite equation to secondary master species
	r_temp->AddTRXN(m1->rxn_primary, 1.0, false);
  r_temp->AddTRXN(m2->rxn_primary, -coef1 / coef2, true);

  return false;
}

bool ModelEngine::GetReactionCoef(Reaction *r, String name, LDBLE &coef, int start)
{
	//
	// Finds coefficient of token in reaction.
	// input: r, pointer to a reaction structure
	//        name, string to find as reaction token
	//        start, search start position in r (default = 1)
	//
	// Return: 0.0, if token not found
	//         coefficient of token, if found.
	//

  coef = 0.0;

	ReactionToken *r_token;	
	for (int i = start; i < r->token_list->Count(); i++)
  {
		r_token = r->token_list->Element(i);

    if (r_token->s->name == name)
    {
      coef = r_token->coef;
      break;
    }
  }

	return true;
}
bool ModelEngine::CheckIn()
{
	//
	// Routine goes through trxn to determine if each master species is in this model.
  // Assumes equation is written in terms of primary and secondary species
  // Checks to see if in is TRUE or REWRITE for each species
  // Returns TRUE if in model FALSE if not
	//

  ReactionToken *t;

	for (int i = 1; i < r_temp->token_list->Count(); i++)
  {
		t = r_temp->token_list->Element(i);

		// Check primary master species in
    if (t->s->primary != NULL && t->s->primary->in)
      continue;

		// Check secondary master species
    if (t->s->secondary != NULL	&& (t->s->secondary->in || t->s->secondary->rewrite))
      continue;

		// Must be primary master species that is out
		return false;
  }

  return true;
}

bool ModelEngine::WriteMassActionEqnX()
{
	//
	// Reduce mass-action equation to the master species that are in the model
	//

  LDBLE coef_e;
  int count;
  int i, count_rxn_orig;


	// Rewrite any secondary master species flagged REWRITE. Replace pe if necessary
	count = 0;
  bool repeat = true;
  while (repeat)
  {
    count++;
    if (count > MAX_ADD_EQUATIONS)
			return false;

    repeat = false;

		ReactionToken *r;
		count_rxn_orig = r_temp->token_list->Count();
    for (i = 1; i < count_rxn_orig; i++)
    {
			r = r_temp->token_list->Element(i);
      if (r->s->secondary == NULL)
				continue;

      if (r->s->secondary->rewrite)
      {
				repeat = true;

				if (!GetReactionCoef(r->s->secondary->rxn_secondary, "e-", coef_e))
					return false;

				r_temp->AddTRXN(r->s->secondary->rxn_secondary,	r->coef, false);

				if (!Equal(coef_e, (LDBLE)0.0, (LDBLE)TOLERANCE))
					r_temp->AddTRXN(*r->s->secondary->pe_rxn, r->coef * coef_e, false);
      }
    }
    
		r_temp->TRXNCombine();
  }

  return true;
}

bool ModelEngine::AddPotentialFactor()
{
	//
	// Add the potential factor to surface mass-action equations.
	// Factor is essentially the activity coefficient, representing
	// the work required to bring charged ions to the surface
	//

  int i;
  LDBLE sum_z;
  Master *master_ptr;
  Unknown *unknown_ptr;

  if (md->use.sur_p->type != DDL)
    return true;

  sum_z = 0.0;
  master_ptr = NULL;

	// Find sum of charge of aqueous species and surface master species
	ReactionToken *r;
	for (i = 1; i < r_temp->token_list->Count(); i++)
  {
		r = r_temp->token_list->Element(i);

    if (r->s->type == AQ || r->s == gd->s_hplus || r->s == gd->s_eminus)
      sum_z += r->s->z * r->coef;

		if (r->s->type == SURF)
      master_ptr = r->s->primary;
  }

	if (master_ptr == NULL)
		return false;
		
	// Find potential unknown for surface species
  unknown_ptr = FindSurfaceChargeUnknown(master_ptr->name, SURF_PSI);

  if (unknown_ptr == NULL)
		return false;

	master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include psi in mass action equation
  if (master_ptr != NULL)
  {
		r = r_temp->token_list->AddNew();

    r->name = master_ptr->s->name;
    r->s = master_ptr->s;
    r->coef = (LDBLE)-2.0 * sum_z;
  }
  else
		return false;

	return true;
}

bool ModelEngine::AddCDMusicFactors(int n)
{
	//
	// Add the potential factors for cd_music to surface mass-action equations.
	// Factors are essentially the activity coefficient, representing
	// the work required to bring charged ions to the three charge layers 
	// of the cd_music model
	//

  int i;
	String name;
  Master *master_ptr;
  Unknown *unknown_ptr;
	ReactionToken *r;

  if (md->use.sur_p->type != CD_MUSIC)
    return true;

  master_ptr = NULL;

	// Find sum of charge of aqueous species and surface master species
	for (i = 1; i < r_temp->token_list->Count(); i++)
	{
		r = r_temp->token_list->Element(i);

    if (r->s->type == SURF)
      master_ptr = r->s->primary;
	}

	if (master_ptr == NULL)
		return false;

	name = master_ptr->e->name;

	// Find potential unknown for surface species
	unknown_ptr = FindSurfaceChargeUnknown (name, SURF_PSI);

  if (unknown_ptr == NULL)
		return false;

	master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include psi in mass action equation
	r = r_temp->token_list->AddNew();

  r->name = master_ptr->s->name;
  r->s = master_ptr->s;
	r->coef = gd->species_list[n]->dz[0];
  
	// Plane 1
  unknown_ptr = FindSurfaceChargeUnknown(name, SURF_PSI1);

  if (unknown_ptr == NULL)
		return false;

  master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include psi in mass action equation
	r = r_temp->token_list->AddNew();

  r->name = master_ptr->s->name;
  r->s = master_ptr->s;
	r->coef = gd->species_list[n]->dz[1];
  
	// Plane 2
  unknown_ptr = FindSurfaceChargeUnknown(name, SURF_PSI2);

  if (unknown_ptr == NULL)
		return false;

  master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include psi in mass action equation
	r = r_temp->token_list->AddNew();

  r->name = master_ptr->s->name;
  r->s = master_ptr->s;
	r->coef = gd->species_list[n]->dz[2];

	return true;
}

Unknown *ModelEngine::FindSurfaceChargeUnknown(String &name, int plane)
{
  int i;
  String token;

	token = name;
	token.Replace ("_", " ");
	token = token(0, " ");

  if (plane == SURF_PSI)
		token += "_CB";
  else if (plane == SURF_PSI1)
		token += "_CBb";
  else if (plane == SURF_PSI2)
		token += "_CBd";

	for (i = 0; i < count_unknowns; i++)
		if (token == (*unknown_list)[i]->name)
      return (*unknown_list)[i];

  return NULL;
}
bool ModelEngine::WriteMBEqnX()
{
	//
	// Rewrite any secondary master species flagged REWRITE
	// Don't add in any pe reactions
	//

  int count;
	bool repeat;
  int i, count_rxn_orig;
  int j, k;
	char *ptr, token[DEFAULT_STRING_LENGTH];
  Master *master_ptr;
	ReactionToken *rxn_token;

  count = 0;
  repeat = true;
  while (repeat)
  {
    count++;
    if (count > MAX_ADD_EQUATIONS)
			return false;

		repeat = false;

		count_rxn_orig = r_temp->token_list->Count();
    for (i = 1; i < count_rxn_orig; i++)
    {
			rxn_token = r_temp->token_list->Element(i);

      if (rxn_token->s->secondary == NULL)
				continue;

			if (rxn_token->s->secondary->rewrite)
      {
				repeat = true;
				r_temp->AddTRXN(rxn_token->s->secondary->rxn_secondary, rxn_token->coef, false);
      }
    }

		r_temp->TRXNCombine();
  }

	eos_list->Clear();
  parent_count = 0;

	ReactionToken *r;
	ElementOfSpecies *eos_p;
  for (i = 1; i < r_temp->token_list->Count(); i++)
  {
		r = r_temp->token_list->Element(i);
    j = eos_list->Count();

		
		r->s->name.Copy(token);
		ptr = &token[0];
		if (!GetElementsInSpecies(&ptr, r->coef))
			return false;

		for (k = j; k < eos_list->Count(); k++)
    {
			eos_p = (*eos_list)[k];
      if (r->s->secondary != NULL)
				master_ptr = r->s->secondary->e->primary;   
      else
      	master_ptr = r->s->primary;
      
      if (eos_p->e == master_ptr->e)
      {
				eos_p->coef = 0.0;
				break;
      }
    }

    if (r->s->secondary == NULL)
			r->s->primary->e->name.Copy(token);
    else
			r->s->secondary->e->name.Copy(token);

		ptr = &token[0];
		if (!GetSecondaryElementsInSpecies(&ptr, r->coef))
			return false;
  }

	CombineElements();

  return true;
}

bool ModelEngine::AddSurfaceChargeBalance()
{
	//
	// Include charge balance in list for mass-balance equations
	//

  int i;
  char *ptr;
	char token[DEFAULT_STRING_LENGTH];

  Master *master_ptr;
  Unknown *unknown_ptr;
	ElementOfSpecies *eos_p;

  if (md->use.sur_p->type != DDL)
    return true;

  master_ptr = NULL;

	// Find master species
	for (i = 0; i < eos_list->Count(); i++)
	{
		eos_p = (*eos_list)[i];

    if (eos_p->e->primary->s->type == SURF)
    {
      master_ptr = eos_p->e->primary;
      break;
    }
	}

	if (i >= eos_list->Count())
		throw exception();

	// Find potential unknown for surface species
	unknown_ptr = FindSurfaceChargeUnknown(master_ptr->e->name, SURF_PSI);

  if (unknown_ptr == NULL)
		throw exception();

	master_ptr = (*unknown_ptr->master)[0];	// potential for surface component 

	master_ptr->e->name.Copy(token);
	ptr = &token[0];
  
	// Include charge balance in list for mass-balance equations
	LDBLE coef = 1.0;
	if (!GetSecondaryElementsInSpecies(&ptr, coef))
		return false;

  return true;
}

bool ModelEngine::AddCDMusicChargeBalances(int n)
{
	//
	// Add the potential factor to surface mass-action equations.
	// Factor is essentially the activity coefficient, representing
	// the work required to bring charged ions to the surface
	//

  int i;
  char *ptr;
	char token[DEFAULT_STRING_LENGTH];

  Master *master_ptr;
  Unknown *unknown_ptr;
	ElementOfSpecies *eos_p;

  if (md->use.sur_p->type != CD_MUSIC)
    return true;

  master_ptr = NULL;

	// Find master species
	for (i = 0; i < eos_list->Count(); i++)
	{
		eos_p = (*eos_list)[i];

		if (eos_p->e->primary->s->type == SURF)
    {
      master_ptr = eos_p->e->primary;
      break;
    }
	}
  
  if (i >= eos_list->Count())
		throw exception();

	// Find potential unknown for plane 0
  unknown_ptr = FindSurfaceChargeUnknown(master_ptr->e->name, SURF_PSI);
  master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include charge balance in list for mass-balance equations
	master_ptr->e->name.Copy(token);
	ptr = &token[0];
	if (!GetSecondaryElementsInSpecies(&ptr, gd->species_list[n]->dz[0]))
		return false;

	// Find potential unknown for plane 1
  unknown_ptr = FindSurfaceChargeUnknown(master_ptr->e->name, SURF_PSI1);
  master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include charge balance in list for mass-balance equations
	master_ptr->e->name.Copy(token);
	ptr = &token[0];
	if (!GetSecondaryElementsInSpecies(&ptr, gd->species_list[n]->dz[1]))
		return false;

	// Find potential unknown for plane 2
  unknown_ptr = FindSurfaceChargeUnknown(master_ptr->e->name, SURF_PSI2);
  master_ptr = (*unknown_ptr->master)[0];	// potential for surface component

	// Include charge balance in list for mass-balance equations
	master_ptr->e->name.Copy(token);
	ptr = &token[0];
	if (!GetSecondaryElementsInSpecies(&ptr, gd->species_list[n]->dz[2]))
		return false;

  return true;
}

bool ModelEngine::MBForSpeciesAQ(int n)
{
	//
	// Make list of mass balance and charge balance equations in which
  //   to insert species n. 
  //
  // count_mb_unknowns - number of equations and summation relations
  // mb_unknowns.unknown - pointer to unknown which contains row number
  // mb_unknowns.source - pointer to the LDBLE number to be multiplied
  //                      by coef, usually moles.
  // mb_unknowns.coef - coefficient of s[n] in equation or relation

  int i, j;
	Species *s;
  Master *master_ptr;
  Unknown *unknown_ptr;

	md->mb_unknowns->Clear();

	s = gd->species_list[n];

	// e- does not appear in any mass balances
	if (s->type == EMINUS)
    return true;

	// Do not include diffuse layer in cb, alk, ah2o, mu
  if (md->charge_balance_unknown != NULL && s->type < H2O)
    StoreMBUnknowns(md->charge_balance_unknown, &s->moles, s->z, &s->dg);

	if (md->alkalinity_unknown != NULL && s->type < H2O)
    StoreMBUnknowns (md->alkalinity_unknown, &s->moles, s->alk, &s->dg);

	if (md->ah2o_unknown != NULL && s->type < H2O)
    StoreMBUnknowns (md->ah2o_unknown, &s->moles, 1.0, &s->dg);

	if (md->mu_unknown != NULL && s->type < H2O)
    StoreMBUnknowns (md->mu_unknown, &s->moles, s->z * s->z, &s->dg);

	// Include diffuse layer in hydrogen and oxygen mass balance
	if (md->mass_hydrogen_unknown != NULL)
  {
    if (md->dl_type_x != NO_DL && md->state >= REACTION)
      StoreMBUnknowns (md->mass_hydrogen_unknown, &s->tot_g_moles, s->h - 2 * s->o, &s->dg_total_g);
    else
      StoreMBUnknowns (md->mass_hydrogen_unknown, &s->moles, s->h - 2 * s->o, &s->dg);
  }

  if (md->mass_oxygen_unknown != NULL)
  {
    if (md->dl_type_x != NO_DL && md->state >= REACTION)
      StoreMBUnknowns (md->mass_oxygen_unknown, &s->tot_g_moles, s->o, &s->dg_total_g);
    else
      StoreMBUnknowns (md->mass_oxygen_unknown, &s->moles, s->o, &s->dg);
  }

	// Sum diffuse layer charge into (surface + DL) charge balance
	SpeciesDiffLayer *sdl_p;
	if (md->use.sur_p != NULL && s->type < H2O && md->dl_type_x != NO_DL)
  {
    j = 0;
		for (i = 0; i < count_unknowns; i++)
    {
			unknown_ptr = (*unknown_list)[i];

      if (unknown_ptr->type == SURFACE_CB)
      {		
				if (md->use.sur_p->type == CD_MUSIC)
					unknown_ptr = (*unknown_list)[i + 2];
	
				sdl_p = (*s->diff_layer)[j];
				StoreMBUnknowns(unknown_ptr, &sdl_p->g_moles, s->z, &sdl_p->dg_g_moles);
				j++;
      }
    }
  }

	// Other mass balances
	ElementOfSpecies *eos;
	for (i = 0; i < eos_list->Count(); i++)
  {
		eos = (*eos_list)[i];

    if (eos->e->master->s->type > AQ && eos->e->master->s->type < SOLID)
      continue;

    master_ptr = eos->e->master;
    
		if (master_ptr->primary && master_ptr->s->secondary != NULL)
			master_ptr = master_ptr->s->secondary;

    if (master_ptr->u == md->ph_unknown)
      continue;
    else if (master_ptr->u == md->pe_unknown)
      continue;
    else if (master_ptr->u == md->charge_balance_unknown)
      continue;
    else if (master_ptr->u == md->alkalinity_unknown)
      continue;
    else if (master_ptr->u->type == SOLUTION_PHASE_BOUNDARY)
      continue;

		if (md->dl_type_x != NO_DL && md->state >= REACTION)
      StoreMBUnknowns(master_ptr->u, &s->tot_g_moles, eos->coef * master_ptr->coef, &s->dg_total_g);
    else
      StoreMBUnknowns(master_ptr->u, &s->moles, eos->coef * master_ptr->coef, &s->dg);
  }

  return true;
}

bool ModelEngine::StoreMBUnknowns(Unknown *unknown_ptr, LDBLE *LDBLE_ptr, LDBLE coef, LDBLE *gamma_ptr)
{
	//
	// Takes an unknown pointer and a coefficient and puts in
	// list of mb_unknowns
	//

  if (Equal(coef, (LDBLE)0.0, (LDBLE)TOLERANCE))
    return true;

	UnknownInfo *m = md->mb_unknowns->AddNew();

  m->u = unknown_ptr;
  m->source = LDBLE_ptr;
  m->gamma_source = gamma_ptr;
  m->coef = coef;

	return true;
}

bool ModelEngine::MBForSpeciesEX(int n)
{
	//
  // Make list of mass balance and charge balance equations in which
  // to insert exchange species n. 
  //
  //      count_mb_unknowns - number of equations and summation relations
  //      mb_unknowns.source - pointer to the LDBLE number to be multiplied
  //                           by coef, usually moles.
  //      mb_unknowns.unknown - pointer to unknown which contains row number
  //      mb_unknowns.coef - coefficient of s[n] in equation or relation
	//

  int i;
  Master *master_ptr;
	Species *s;

 	md->mb_unknowns->Clear();

	// Master species for exchange do not appear in any mass balances
	s = gd->species_list[n];
	if (s->type == EX && s->primary != NULL)
    return true;

	// Include diffuse layer in hydrogen and oxygen mass balance
  if (md->charge_balance_unknown != NULL)
		StoreMBUnknowns(md->charge_balance_unknown, &s->moles, s->z, &s->dg);

	if (md->mass_hydrogen_unknown != NULL)
		StoreMBUnknowns(md->mass_hydrogen_unknown, &s->moles, s->h - 2 * s->o, &s->dg);

	if (md->mass_oxygen_unknown != NULL)
    StoreMBUnknowns(md->mass_oxygen_unknown, &s->moles, s->o, &s->dg);

	// Other mass balances
	ElementOfSpecies *eos_p;
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos_p = (*eos_list)[i];
    if (eos_p->e->master->s->type > AQ && eos_p->e->master->s->type < SOLID)
      continue;

    master_ptr = eos_p->e->master;

		if (master_ptr->primary && master_ptr->s->secondary != NULL)
			master_ptr = master_ptr->s->secondary;
    
		// Special for ph_unknown, pe_unknown, and alkalinity_unknown
    if (master_ptr->u == md->ph_unknown)
      continue;
    else if (master_ptr->u == md->pe_unknown)
      continue;
    else if (master_ptr->u == md->alkalinity_unknown)
      continue;

		// EX, sum exchange species only into EXCH mass balance in initial calculation
		// into all mass balances in reaction calculation
    if (md->state >= REACTION || master_ptr->s->type == EX)
      StoreMBUnknowns(master_ptr->u, &s->moles, eos_p->coef * master_ptr->coef, &s->dg);
  }

  return true;
}

bool ModelEngine::MBForSpeciesSURF(int n)
{
	//
  // Make list of mass balance and charge balance equations in which
  // to insert species n. 
  //
  //      count_mb_unknowns - number of equations and summation relations
  //      mb_unknowns.source - pointer to the LDBLE number to be multiplied
  //                           by coef, usually moles.
  //      mb_unknowns.unknown - pointer to unknown which contains row number
  //      mb_unknowns.coef - coefficient of s[n] in equation or relation
	//

  int i;
  Master *master_ptr;

	md->mb_unknowns->Clear();

	// Include in charge balance, if diffuse_layer_x == false
	Species *s = gd->species_list[n];
  if (md->charge_balance_unknown != NULL && md->dl_type_x == NO_DL)
    StoreMBUnknowns(md->charge_balance_unknown, &s->moles, s->z, &s->dg);

	// Include diffuse layer in hydrogen and oxygen mass balance
  if (md->mass_hydrogen_unknown != NULL)
    StoreMBUnknowns (md->mass_hydrogen_unknown, &s->moles, s->h - 2 * s->o, &s->dg);

	if (md->mass_oxygen_unknown != NULL)
    StoreMBUnknowns (md->mass_oxygen_unknown, &s->moles, s->o, &s->dg);

	// Other mass balances
	ElementOfSpecies *eos_p;
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos_p = (*eos_list)[i];

		// Skip H+, e-, and H2O
    if (eos_p->e->master->s->type > AQ && eos_p->e->master->s->type < SOLID)
      continue;

    master_ptr = eos_p->e->master;
    if (master_ptr->primary && master_ptr->s->secondary != NULL)
			master_ptr = master_ptr->s->secondary;

		// SURF_PSI, sum surface species in (surface + DL) charge balance
    if (master_ptr->s->type == SURF_PSI && md->use.sur_p->type != CD_MUSIC)
    {
      StoreMBUnknowns (master_ptr->u, &s->moles, s->z, &s->dg);
      continue;
    }

    if (master_ptr->s->type == SURF_PSI && md->use.sur_p->type == CD_MUSIC)
    {
      StoreMBUnknowns (master_ptr->u, &s->moles, s->dz[0], &s->dg);
      continue;
    }

    if (master_ptr->s->type == SURF_PSI1)
    {
      StoreMBUnknowns (master_ptr->u, &s->moles, s->dz[1], &s->dg);
      continue;
    }

    if (master_ptr->s->type == SURF_PSI2)
    {
      StoreMBUnknowns (master_ptr->u, &s->moles, s->dz[2], &s->dg);
      continue;
    }

		// Special for ph_unknown, pe_unknown, and alkalinity_unknown
    if (master_ptr->u == md->ph_unknown)
      continue;
    else if (master_ptr->u == md->pe_unknown)
      continue;
    else if (master_ptr->u == md->alkalinity_unknown)
      continue;

		// SURF, sum surface species only into SURFACE mass balance in initial calculation
		// into all mass balances in reaction calculation
    if (md->state >= REACTION || master_ptr->s->type == SURF)
    {
      StoreMBUnknowns (master_ptr->u, &s->moles, eos_p->coef * master_ptr->coef, &s->dg);
    }
  }

  return true;
}

bool ModelEngine::BuildMBSums()
{
	//
	// Function builds lists sum_mb1 and sum_mb2  that describe how to sum molalities
	// to calculate mass balance sums, including activity of water, ionic strength,
	// charge balance, and alkalinity.
	//

  int i;
  //LDBLE *target;
	UnknownInfo *ui;
	int count = md->mb_unknowns->Count();

	for (i = 0; i < count; i++)
  {
		ui = (*md->mb_unknowns)[i];
    StoreMB(ui->source, &(ui->u->f), ui->coef);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.buildjacobiansums_f, "BuildMBSums-1) i:%d %d %s %20.20e\n",
																			i,
																			ui->u->number,
																			ui->u->name.CharPtr(),
																			ui->coef);
#endif
  }

	return true;
}

bool ModelEngine::BuildJacobianSums(int k)
{
	//
	// Function builds lists sum_jacob1 and sum_jacob2 that describe how to sum molalities
	// to form jacobian.
	//

	int i, j, kk;
	int index;
  int count_g;
  LDBLE coef;
  //LDBLE *source, *target;
	UnknownInfo *mb_unknown;
	SpeciesDiffLayer *sdl_p;
	Unknown *u, *uprior, *u_kk;

	Species *s = gd->species_list[k];

	// Calculate jacobian coefficients for each mass balance equation
	for (i = 0; i < md->mb_unknowns->Count(); i++)
  {
		mb_unknown = (*md->mb_unknowns)[i];

		// Store d(moles) for a mass balance equation

		// initial solution only
    if (mb_unknown->u->type == SOLUTION_PHASE_BOUNDARY)
      continue;

		coef = mb_unknown->coef;

		StoreDN(k, mb_unknown->source, mb_unknown->u->number, coef, mb_unknown->gamma_source);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.buildjacobiansums_f, "1) k:%d i:%d u_number:%d: %20.20e\n",
																			k,
																			i, 
																			mb_unknown->u->number,
																			coef);
#endif

		// Add extra terms for change in dg/dx in diffuse layer model
		if (s->type >= H2O || md->dl_type_x == NO_DL)
      continue;
    else if ((mb_unknown->u->type == MB || mb_unknown->u->type == MH || mb_unknown->u->type == MH2O) && md->state >= REACTION)
    {
      if (md->mass_oxygen_unknown != NULL)
      {
				// term for water, sum of all surfaces
				index = mb_unknown->u->number * (count_unknowns + 1) + md->mass_oxygen_unknown->number;
				StoreJacob(&s->tot_dh2o_moles, &arr[index], coef);

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.buildjacobiansums_f, "2) index:%d s_number:%d u_number:%d: %d %20.20e\n",
																					index,
																					s->number, 
																					md->mass_oxygen_unknown->number,
																					coef);
#endif
      }

			// terms for psi, one for each surface
      count_g = 0;
      for (j = 0; j < count_unknowns; j++)
      {
				u = (*unknown_list)[j];
				sdl_p = (*s->diff_layer)[count_g];

				if (u->type != SURFACE_CB)
					continue;

				index = mb_unknown->u->number * (count_unknowns + 1) + u->number;
				StoreJacob(&sdl_p->dx_moles, &arr[index], coef);
				count_g++;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.buildjacobiansums_f, "3) index:%d %d %d %d %20.20e\n",
																					index,
																					mb_unknown->u->number,
																					u->number,
																					count_g, 
																					coef);
#endif

				if (count_g >= md->use.sur_p->charge->Count())
					break;
      }

			// terms for related phases
      count_g = 0;
      for (j = 0; j < count_unknowns; j++)
      {
				u = (*unknown_list)[j];
				sdl_p = (*s->diff_layer)[count_g];

				if (u->type != SURFACE_CB)
					continue;

				uprior = (*unknown_list)[j - 1];

				// has related phase
				if (uprior->surface_comp->phase_name.IsEmpty())
					continue;

				// now find the related phase
				for (kk = count_unknowns - 1; kk >= 0; kk--)
				{
					u_kk = (*unknown_list)[kk];

					if (u_kk->type != PP)
						continue;

					if (u_kk->p->name == uprior->surface_comp->phase_name)
						break;
				}

				if (kk >= 0)
				{
					index = mb_unknown->u->number * (count_unknowns + 1) + (*unknown_list)[kk]->number;
					StoreJacob(&sdl_p->drelated_moles, &arr[index], coef);

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.buildjacobiansums_f, "4) index:%d %d %d %20.20e\n",
																						index,
																						mb_unknown->u->number,
																						(*unknown_list)[kk]->number,
																						coef);
#endif
				}

				count_g++;
				if (count_g >= md->use.sur_p->charge->Count())
					break;
      }
    }
    else if (mb_unknown->u->type == SURFACE_CB)
    {
      count_g = 0;

			for (j = 0; j < count_unknowns; j++)
      {
				u = (*unknown_list)[j];
				sdl_p = (*s->diff_layer)[count_g];

				if (u->type != SURFACE_CB)
					continue;

				if (mb_unknown->u->number == u->number)
				{
					index = mb_unknown->u->number * (count_unknowns + 1) + u->number;
					StoreJacob(&sdl_p->dx_moles, &arr[index], coef);

					uprior = (*unknown_list)[j - 1];

					// term for related phase 
					// has related phase 
					if (!uprior->surface_comp->phase_name.IsEmpty())
					{
						// now find the related phase
						for (kk = count_unknowns - 1; kk >= 0; kk--)
						{
							u_kk = (*unknown_list)[kk];

							if (u_kk->type != PP)
								continue;

							if (u_kk->p->name == uprior->surface_comp->phase_name)
								break;
						}

						if (kk >= 0)
						{
							index = mb_unknown->u->number * (count_unknowns + 1) + u_kk->number;
							StoreJacob(&sdl_p->drelated_moles, &arr[index], coef);
						}
					}

					if (md->mass_oxygen_unknown != NULL)
					{
						// term for water, for same surfaces
						index = mb_unknown->u->number * (count_unknowns + 1) + md->mass_oxygen_unknown->number;
						StoreJacob(&sdl_p->dh2o_moles, &arr[index], coef);
					}
					break;
				}

				count_g++;

				if (count_g >= md->use.sur_p->charge->Count())
					break;
      }
    }
  }

  return true;
}

bool ModelEngine::WriteMBForSpeciesList(int n)
{
	//
	// Sets up data to add to species_list
	// Original secondary redox states are retained
	//

  int i;
	char *ptr, token[DEFAULT_STRING_LENGTH];
	ElementOfSpecies *new_eos;

	r_temp->Reset();
	r_temp->AddTRXN(gd->species_list[n]->rxn_s, 1.0, false);

	//Copy to elt_list
	eos_list->Clear();
  parent_count = 0;

	ReactionToken *r;
  for (i = 1; i < r_temp->token_list->Count(); i++)
  {
		r = r_temp->token_list->Element(i);

    if (r->s->secondary == NULL)
    {
			r->s->primary->e->name.Copy(token);
      ptr = &token[0];
			GetSecondaryElementsInSpecies(&ptr, r->coef);
    }
    else
    {
			r->s->secondary->e->name.Copy(token);
      ptr = &token[0];
			GetSecondaryElementsInSpecies(&ptr, r->coef);
    }
  }

	ElementOfSpecies *eos;
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos = (*eos_list)[i];
    if (eos->e->name == "O(-2)")
    {
			new_eos = eos_list->AddNew();

			new_eos->e = gd->e_h_one;
      new_eos->coef = eos->coef * 2;
    }
  }

	CombineElements();

	eos_list->CopyTo(gd->species_list[n]->eos_sys_total);
	
  return true;
}

bool ModelEngine::BuildSpeciesList(int n)
{
	//
	// Builds a list that includes an entry for each master species in each
	// secondary reaction. Used for summing species of each element and 
	// printing results.
	//

  int j;
	SpeciesInfo *new_s;
	ElementOfSpecies *eos_p;
	Species *s = gd->species_list[n];

	// Treat species made only with H+, e-, and H2O specially 
	if (IsSpecial(s))
  {
		new_s = md->species_info_list->AddNew();

    new_s->master_s = gd->s_hplus;
    new_s->s = s;
    new_s->coef = 0.0;

		return true;
  }

	// Treat exchange species specially
  if (s->type == EX)
  {
    if (s->primary != NULL)
      return true;		// master species has md->zero molality

    for (j = 0; j < eos_list->Count(); j++)
    {
			eos_p = (*eos_list)[j];

      if (eos_p->e->master->s->type != EX)
				continue; 

			new_s = md->species_info_list->AddNew();

			new_s->master_s = eos_p->e->master->s;
			new_s->s = s;
			new_s->coef = eos_p->e->master->coef * eos_p->coef;
    }

    return true;
  }

	// Treat surface species specially
  if (s->type == SURF_PSI)
    return true;

  if (s->type == SURF)
  {
    for (j = 0; j < eos_list->Count(); j++)
    {
			eos_p = (*eos_list)[j];

      if (eos_p->e->master->s->type != SURF)
				continue;

      new_s = md->species_info_list->AddNew();

			new_s->master_s = eos_p->e->master->s;
			new_s->s = s;
			new_s->coef = eos_p->e->master->coef * eos_p->coef;
    }

    return true;
  }

	// Other aqueous species
	Master *master_ptr;
  for (j = 0; j < eos_list->Count(); j++)
  {
		eos_p = (*eos_list)[j];

    if (IsSpecial(eos_p->e->master->s))
      continue;

    if (eos_p->e->master->s->secondary != NULL)
      master_ptr = eos_p->e->master->s->secondary;
    else
      master_ptr = eos_p->e->master->s->primary;

    new_s = md->species_info_list->AddNew();

		new_s->master_s = master_ptr->s;
		new_s->s = s;
		new_s->coef = master_ptr->coef * eos_p->coef;
  }

  return true;
}

bool ModelEngine::StoreMB(LDBLE * source, LDBLE * target, LDBLE coef)
{
	//
	// Adds item to list sum_mb1 or sum_mb2
	// If coef is 1.0, adds to sum_mb1, which does not require a multiply
	// else, adds to sum_mb2, which will multiply by coef
	//

	STCoef *new_item;

  if (Equal(coef, (LDBLE)1.0, (LDBLE)TOLERANCE))
  {
		new_item = md->sum_mb1->AddNew();

		new_item->source = source;
		new_item->target = target;
  }
  else
  {
		new_item = md->sum_mb2->AddNew();

		new_item->source = source;
		new_item->target = target;
		new_item->coef = coef;
  }

  return true;
}

bool ModelEngine::WritePhaseSysTotal(int n)
{
	//
	// Sets up data to add to species_list
	// Original secondary redox states are retained
	//

  int i;
	ReactionToken *r;
	char *ptr, token[DEFAULT_STRING_LENGTH];
	Phase *p = gd->phase_list[n];

	// Start with secondary reaction
	r_temp->Reset();
	r_temp->AddTRXN(p->rxn_s, 1.0, false);

	// Copy to elt_list
	eos_list->Clear();
  parent_count = 0;

	for (i = 1; i < r_temp->token_list->Count(); i++)
  {
		r = r_temp->token_list->Element(i);

    if (r->s->secondary == NULL)
			r->s->primary->name.Copy(token); //master and element have the same name (change made in the load of database)
    else
			r->s->secondary->name.Copy(token);

    ptr = &token[0];
		GetSecondaryElementsInSpecies(&ptr, r->coef);
  }

	ElementOfSpecies *eos, *new_eos;
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos = (*eos_list)[i];

    if (eos->e->name == "O(-2)")
    {
			new_eos = eos_list->AddNew();

      new_eos->e = gd->e_h_one;
      new_eos->coef = eos->coef * 2;
    }
  }

	CombineElements();

	eos_list->CopyTo(p->eos_sys_total);

	return true;
}

bool ModelEngine::BuildSolutionPhaseBoundaries()
{
	//
	// Build into sums the logic to calculate inverse saturation indices for
	// solution phase boundaries
	//

  int i, j;
  Master *master_ptr;
  ReactionToken *rxn_ptr;
	Unknown *x;

  if (md->solution_phase_boundary_unknown == NULL)
    return true;

	// Calculate inverse saturation index
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];

    if (x->type != SOLUTION_PHASE_BOUNDARY)
      continue;

    StoreMB(&x->p->lk, &x->f, 1.0);
    StoreMB(&x->si, &x->f, 1.0);

    if (!x->p->in)
			return false;

		for (j = 1; j < x->p->rxn_x->token_list->Count(); j++)
		{
			rxn_ptr = x->p->rxn_x->token_list->Element(j); 
      StoreMB(&rxn_ptr->s->la, &x->f, -rxn_ptr->coef);
		}
  }
	
	// Put coefficients into array
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];

    if (x->type != SOLUTION_PHASE_BOUNDARY)
      continue;

		for (j = 1; j < x->p->rxn_x->token_list->Count(); j++)
		{
			rxn_ptr = x->p->rxn_x->token_list->Element(j);

			if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
				master_ptr = rxn_ptr->s->secondary;
      else
				master_ptr = rxn_ptr->s->primary;

			if (master_ptr->u == NULL)
				continue;

      StoreJacob0(x->number, master_ptr->u->number, rxn_ptr->coef);
    }
  }

  return true;
}



bool ModelEngine::BuildMinExch()
{
	//
	// Defines proportionality factor between mineral and exchanger to jacob0
	//

  int i, j, k, jj;
  int row;
  ExchComp *comp_ptr;
  Master *master_ptr, *m0;
  Unknown *unknown_ptr, *u_j, *u_k;
	ElementOfSpecies *eos_p;
  LDBLE coef;

	if (md->use.exc_p == NULL)
    return true;
		
  if (!md->use.exc_p->related_phases)
    return true;
		
	for (i = 0; i < md->use.exc_p->comps->Count(); i++)
  {
		comp_ptr = (*md->use.exc_p->comps)[i];

		if (comp_ptr->phase_name.IsEmpty())
      continue;
			
    // find unknown number
    for (j = count_unknowns - 1; j >= 0; j--)
    {
			u_j = (*unknown_list)[j];
			m0 = (*u_j->master)[0];

			if (u_j->type != EXCH)
				continue;
				
      if (m0 == comp_ptr->master)
				break;
    }
		
    for (k = count_unknowns - 1; k >= 0; k--)
    {
			u_k = (*unknown_list)[k];

			if (u_k->type != PP)
				continue;
				
      if (u_k->p->name == comp_ptr->phase_name)
				break;
    }
		
    if (j == -1)
			return false;
		
    if (k == -1)
      continue;
			
		// Build jacobian

    // charge balance
    StoreJacob0(md->charge_balance_unknown->number, u_k->number, comp_ptr->formula_z * comp_ptr->phase_proportion);
		StoreSumDeltas(&delta[k], &md->charge_balance_unknown->delta, -comp_ptr->formula_z * comp_ptr->phase_proportion);
		
    // mole balance balance
    eos_list->Clear();
    parent_count = 0;
    
		CopyToTempEOSList(comp_ptr->formula_totals, 1.0);
    ChangeHydrogenInEOSList(0);
		
    for (jj = 0; jj < eos_list->Count(); jj++)
    {
			eos_p = (*eos_list)[jj];
      master_ptr = (*eos_list)[jj]->e->primary;
			
			if (!master_ptr->in)
				master_ptr = master_ptr->s->secondary;
			
      if (master_ptr == NULL)
				return false;
			
      if (master_ptr->s->type == EX)
      {
				if (!Equal(u_j->moles, u_k->moles * eos_p->coef * comp_ptr->phase_proportion, (LDBLE)5.0 * md->convergence_tolerance))
				  u_j->moles = u_k->moles * eos_p->coef * comp_ptr->phase_proportion;
      }
			
      coef = eos_p->coef;
			
      if (master_ptr->s == gd->s_hplus)
      {
				row = md->mass_hydrogen_unknown->number;
				unknown_ptr = md->mass_hydrogen_unknown;
      }
      else if (master_ptr->s == gd->s_h2o)
      {
				row = md->mass_oxygen_unknown->number;
				unknown_ptr = md->mass_oxygen_unknown;
      }
      else
      {
				row = master_ptr->u->number;
				unknown_ptr = master_ptr->u;
      }
			
      StoreJacob0 (row, u_k->number, coef * comp_ptr->phase_proportion);
      StoreSumDeltas(&delta[k], &unknown_ptr->delta, -coef * comp_ptr->phase_proportion);
    }
  }
	
  return true;
}

bool ModelEngine::BuildMinSurface()
{
	//
	//   Defines proportionality factor between mineral and surface to jacob0
	//

  int i, j, k, jj, row;
  List<ElementOfSpecies> *next_elt;
	ElementOfSpecies *eos_p;
  SurfaceComp *comp_ptr, *sc_p;
  Unknown *unknown_ptr, *u_j, *u_k, *u;
  Master *master_ptr, *m;
  LDBLE coef;

  if (md->use.sur_p == NULL)
    return true;
		
	if (!md->use.sur_p->related_phases)
		return true;
		
	for (i = 0; i < md->use.sur_p->comps->Count(); i++)
  {
		sc_p = (*md->use.sur_p->comps)[i];

		if (sc_p->phase_name.IsEmpty())
      continue;
			
    // find unknown number
    for (j = count_unknowns - 1; j >= 0; j--)
    {
			u_j = (*unknown_list)[j];
			m = (*u_j->master)[0];

      if (u_j->type != SURFACE)
				continue;
				
      if (m == sc_p->master)
				break;
    }
		
    for (k = count_unknowns - 1; k >= 0; k--)
    {
			u_k = (*unknown_list)[k];

      if (u_k->type != PP)
				continue;
				
      if (u_k->p->name == sc_p->phase_name)
				break;
    }
		
    if (j == -1)
			return false;
		
    if (k == -1)
      continue;
			
    comp_ptr = u_j->surface_comp;
		
    if (md->use.sur_p->type == CD_MUSIC)
    {
      // Add formula for CD_MUSIC
      next_elt = comp_ptr->formula_totals;
    }
    else
    {
      // Add master species for non CD_MUSIC
			next_elt = (*u_j->master)[0]->s->eos_list;
    }
		
    // update grams == moles in this case
		u = (*unknown_list)[j + 1];
    if (j < count_unknowns - 1 && u->type == SURFACE_CB)
      StoreSumDeltas (&delta[k], &u->related_moles, -1.0);
		
    // charge balance
    StoreJacob0 (md->charge_balance_unknown->number, u_k->number, comp_ptr->formula_z * comp_ptr->phase_proportion);
    StoreSumDeltas (&delta[k], &md->charge_balance_unknown->delta, -comp_ptr->formula_z * comp_ptr->phase_proportion);
		
    eos_list->Clear();
    parent_count = 0;
		
    CopyToTempEOSList(next_elt, 1.0);
    ChangeHydrogenInEOSList(0);
		
		for (jj = 0; jj < eos_list->Count(); jj++)
    {
			eos_p = (*eos_list)[jj];

      master_ptr = eos_p->e->primary;
      if (!master_ptr->in)
				master_ptr = master_ptr->s->secondary;

			if (master_ptr == NULL)
				return false;

      if (master_ptr->s->type == SURF)
      {
				if (!Equal(u_j->moles, u_k->moles * eos_p->coef * comp_ptr->phase_proportion, (LDBLE)5.0 * md->convergence_tolerance))
				  u_j->moles = u_k->moles * eos_p->coef * comp_ptr->phase_proportion;
      }
			
      coef = eos_p->coef;
			
      if (master_ptr->s == gd->s_hplus)
      {
				row = md->mass_hydrogen_unknown->number;
				unknown_ptr = md->mass_hydrogen_unknown;
      }
      else if (master_ptr->s == gd->s_h2o)
      {
				row = md->mass_oxygen_unknown->number;
				unknown_ptr = md->mass_oxygen_unknown;
      }
      else
      {
				row = master_ptr->u->number;
				unknown_ptr = master_ptr->u;
      }
			
      StoreJacob0(row, u_k->number, coef * comp_ptr->phase_proportion);
      StoreSumDeltas(&delta[k], &unknown_ptr->delta, -coef * comp_ptr->phase_proportion);
    }
  }
	
  return true;
}

bool ModelEngine::BuildGasPhase()
{
	//
	// Put coefficients into lists to sum iaps to test for equilibrium
	// Put coefficients into lists to build jacobian for 
	// sum of partial pressures equation and
	// mass balance equations for elements contained in gases
	//

  int i, j, k;
  int row, col;
  Master *master_ptr;
  ReactionToken *rxn_ptr;
  GasComp *gas_comp_ptr;
  Phase *phase_ptr;
  Unknown *unknown_ptr;
  LDBLE coef, coef_elt;
	ElementOfSpecies *eos_p;

  if (md->gas_unknown == NULL)
    return true;
		
	for (i = 0; i < md->use.gas_p->comps->Count(); i++)
  {
		// Determine elements in gas component
    eos_list->Clear();
    parent_count = 0;

    gas_comp_ptr = (*md->use.gas_p->comps)[i];
    phase_ptr = gas_comp_ptr->phase;
		
		if (phase_ptr->rxn_x->token_list->Count() <= 0)
      continue;
			
    CopyToTempEOSList(phase_ptr->eos_list, 1.0);
    ChangeHydrogenInEOSList(0);
		
		// Build mass balance sums for each element in gas

		// All elements in gas
    for (j = 0; j < eos_list->Count(); j++)
    {
			eos_p = (*eos_list)[j];

      unknown_ptr = NULL;
			if (eos_p->e->name == "H")
				unknown_ptr = md->mass_hydrogen_unknown;
			else if (eos_p->e->name == "O")
				unknown_ptr = md->mass_oxygen_unknown;
      else
      {
				if (eos_p->e->primary->in)
				  unknown_ptr = eos_p->e->primary->u;
				else if (eos_p->e->primary->s->secondary != NULL)
				  unknown_ptr = eos_p->e->primary->s->secondary->u;
      }
			
      if (unknown_ptr != NULL)
      {
				coef = eos_p->coef;
				StoreMB(&(gas_comp_ptr->phase->moles_x), &(unknown_ptr->f), coef);
      }
    }
		
    if (md->use.gas_p->type == PRESSURE)
    {
      // Total pressure of gases
      StoreMB(&(gas_comp_ptr->phase->p_soln_x), &(md->gas_unknown->f), 1.0);
    }
		
		// Build jacobian sums for mass balance equations
    for (j = 0; j < eos_list->Count(); j++)
    {
			eos_p = (*eos_list)[j];

      unknown_ptr = NULL;
			
			if (eos_p->e->name == "H")
				unknown_ptr = md->mass_hydrogen_unknown;
			else if (eos_p->e->name == "O")
				unknown_ptr = md->mass_oxygen_unknown;
      else
      {
				if (eos_p->e->primary->in)
				  unknown_ptr = eos_p->e->primary->u;
				else if (eos_p->e->primary->s->secondary != NULL)
				  unknown_ptr = eos_p->e->primary->s->secondary->u;
      }
			
      if (unknown_ptr == NULL)
				continue;
			
      row = unknown_ptr->number * (count_unknowns + 1);
      coef_elt = eos_p->coef;
			
			for (k = 1; k < phase_ptr->rxn_x->token_list->Count(); k++)
      {
				rxn_ptr = phase_ptr->rxn_x->token_list->Element(k);

				if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
				  master_ptr = rxn_ptr->s->secondary;
				else
				  master_ptr = rxn_ptr->s->primary;
				
				if (master_ptr == NULL)
					return false;
				
				if (master_ptr->u == NULL)
				  continue;
				
				if (!master_ptr->in)
					return false;
				
				col = master_ptr->u->number;
				coef = coef_elt * rxn_ptr->coef;
				StoreJacob (&(gas_comp_ptr->phase->moles_x), &(arr[row + col]), coef);				
      }
			
      if (md->use.gas_p->type == PRESSURE)
      {
				// derivative wrt total moles of gas
				StoreJacob(&(gas_comp_ptr->phase->fraction_x), &(arr[row + md->gas_unknown->number]), coef_elt);
      }
    }
		
		// Build jacobian sums for sum of partial pressures equation
    if (md->use.gas_p->type != PRESSURE)
      continue;
			
    unknown_ptr = md->gas_unknown;
    row = unknown_ptr->number * (count_unknowns + 1);
		
		for (k = 1; k < phase_ptr->rxn_x->token_list->Count(); k++)
    {
			rxn_ptr = phase_ptr->rxn_x->token_list->Element(k);

			if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
				master_ptr = rxn_ptr->s->secondary;
      else
				master_ptr = rxn_ptr->s->primary;
			
      if (master_ptr->u == NULL)
				continue;

			if (!master_ptr->in)
				return false;
			
      col = master_ptr->u->number;
      coef = rxn_ptr->coef;
      StoreJacob (&(gas_comp_ptr->phase->p_soln_x), &(arr[row + col]), coef);
    }
  }
	
  return true;
}

bool ModelEngine::BuildSSAssemblage()
{
	//
	// Put coefficients into lists to sum iaps to test for equilibrium
	// Put coefficients into lists to build jacobian for 
	//    mass action equation for component
	//    mass balance equations for elements contained in solid solutions
	//

  int i, j, k, l, m;
	bool stop;
  int row, col;
  Master *master_ptr, *m_2, *m0;
  ReactionToken *rxn_ptr;
	Unknown *u, *u2;
  SS *s_s_ptr, *s_s_ptr_old;
	//SSComp *ssc_p;
	ElementOfSpecies *eos_p;
  char token[MAX_LENGTH];
  char *ptr;

  if (md->s_s_unknown == NULL)
    return true;
		
  s_s_ptr_old = NULL;
  col = 0;
	
  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

    if (u->type != S_S_MOLES)
      continue;
			
    s_s_ptr = u->s_s;

    if (s_s_ptr != s_s_ptr_old)
    {
      col = u->number;
      s_s_ptr_old = s_s_ptr;
    }
		
		// Calculate function value (inverse saturation index)
    if (u->p->rxn_x->token_list->Count() <= 0)
      continue;
			
		StoreMB(&(u->p->lk), &(u->f), 1.0);
		
#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.buildssassemblage_f, "1) u_number:%d ss:%s u_p:%s\n", u->number, s_s_ptr->name.CharPtr(), u->p->name.CharPtr());
#endif

		for (m = 1; m < u->p->rxn_x->token_list->Count(); m++)
    {
			rxn_ptr = u->p->rxn_x->token_list->Element(m);
      StoreMB(&(rxn_ptr->s->la), &(u->f), -rxn_ptr->coef);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.buildssassemblage_f, "2) s_number:%d %20.20e\n", rxn_ptr->s->number, rxn_ptr->coef);
#endif		
		}
		
    // include mole fraction 
    StoreMB(&(u->p->log10_fraction_x), &(u->f), 1.0);
		
    // include activity coeficient
    StoreMB(&(u->p->log10_lambda), &(u->f), 1.0);
		
		// Put coefficients into mass action equations
    // first IAP terms 
		for (m = 1; m < u->p->rxn_x->token_list->Count(); m++)
    {
			rxn_ptr = u->p->rxn_x->token_list->Element(m);

			if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
				master_ptr = rxn_ptr->s->secondary;
      else
				master_ptr = rxn_ptr->s->primary;
			
      if (master_ptr == NULL || master_ptr->u == NULL)
				continue;
				
      StoreJacob0 (u->number, master_ptr->u->number, rxn_ptr->coef);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.buildssassemblage_f, "3) m_number:%d m_u_number:%d %20.20e\n", master_ptr->number, master_ptr->u->number, rxn_ptr->coef);
#endif		
    }
		
    if (s_s_ptr->a0 != 0.0 || s_s_ptr->a1 != 0.0)
    {
			// For binary solid solution
      // next dnc terms
      row = u->number * (count_unknowns + 1);
			
      if (u->s_s_comp_number == 0)
				col = u->number;
      else
				col = u->number - 1;
			
      StoreJacob(&(u->p->dnc), &(arr[row + col]), -1);
			
      // next dnb terms
      col++;
      StoreJacob(&(u->p->dnb), &(arr[row + col]), -1);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.buildssassemblage_f, "4) row:%d col1:%d col2:%d\n", row, col - 1, col);
#endif	
    }
    else
    {
			// For ideal solid solution
      row = u->number * (count_unknowns + 1);
			
			for (j = 0; j < s_s_ptr->comps_list->Count(); j++)
      {
				if (j != u->s_s_comp_number)
				  StoreJacob(&(u->p->dn), &(arr[row + col + j]), -1.0);
				else
				  StoreJacob(&(u->p->dnb), &(arr[row + col + j]), -1.0);

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.buildssassemblage_f, "5) row:%d col:%d j:%d s_s_comp_number:%d\n", row, col, j, u->s_s_comp_number);
#endif	
      }
    }
		
		// Put coefficients into mass balance equations
    eos_list->Clear();
    parent_count = 0;

    u->p->formula.Copy(token);
    ptr = &token[0];
		GetElementsInSpecies(&ptr, 1.0);
		
		// Go through elements in phase
	  ChangeHydrogenInEOSList(0);
		
    for (j = 0; j < eos_list->Count(); j++)
    {
			eos_p = (*eos_list)[j];

			if (eos_p->e->name == "H" && md->mass_hydrogen_unknown != NULL)
      {
				StoreJacob0(md->mass_hydrogen_unknown->number, u->number, -eos_p->coef);
				StoreSumDeltas (&(delta[i]), &md->mass_hydrogen_unknown->delta, eos_p->coef);

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.buildssassemblage_f, "6) %d %d %20.20e\n", 
																						 i,
																						 md->mass_hydrogen_unknown->number,
																						 eos_p->coef);
#endif	
      }
      else if (eos_p->e->name == "O" && md->mass_oxygen_unknown != NULL)
      {
				StoreJacob0(md->mass_oxygen_unknown->number, u->number, -eos_p->coef);
				StoreSumDeltas(&(delta[i]), &md->mass_oxygen_unknown->delta, eos_p->coef);
      
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.buildssassemblage_f, "7) %d %d %20.20e\n", 
																						 i,
																						 md->mass_oxygen_unknown->number,
																						 eos_p->coef);
#endif				
			}
      else
      {
				master_ptr = eos_p->e->primary;
				
				if (!master_ptr->in)
				  master_ptr = master_ptr->s->secondary;
				
				if (master_ptr == NULL || !master_ptr->in)
					return false;
				else if (master_ptr->in)
				{
					// Master species is in model
				  StoreJacob0 (master_ptr->u->number, u->number, -eos_p->coef);
				  StoreSumDeltas (&delta[i], &master_ptr->u->delta, eos_p->coef);

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.buildssassemblage_f, "8) %d %d %20.20e\n", 
																							 i,
																							 master_ptr->u->number,
																							 eos_p->coef);
#endif	
				}
				else if (master_ptr->rewrite)
				{
					// Master species in equation needs to be rewritten
				  stop = false;
					
				  for (k = 0; k < count_unknowns; k++)
				  {
						u2 = (*unknown_list)[k];						

				    if (u2->type != MB)
				      continue;
							
				    for (l = 0; l < u2->master->Count() != NULL; l++)
				    {
							m_2 = (*u2->master)[l];

				      if (m_2 == master_ptr)
				      {
								m0 = (*u2->master)[0];

								StoreJacob0(m0->u->number, u->number, -eos_p->coef);
								StoreSumDeltas (&delta[i], &m0->u->delta, eos_p->coef);
								
#ifdef DEBUG_MOHID
								if (d.debug_status)
									fprintf(d.buildssassemblage_f, "9) %d %d %20.20e\n", 
																										 i,
																										 m0->u->number,
																										 eos_p->coef);
#endif	
								stop = true;
								break;
				      }
				    }
						
				    if (stop)
				      break;
				  }
				}
      }
    }
  }
	
#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.buildssassemblage_f, "FIM)-------------\n");
#endif	
	return true;
}

bool ModelEngine::SaveModel()
{
	return true;
}

bool ModelEngine::StoreDN(int k, LDBLE * source, int row, LDBLE coef_in, LDBLE * gamma_source)
{
	//
	// Stores the terms for d moles of species k in solution into row, multiplied
	// by coef_in
	//

  int col, j;
  LDBLE coef;
  ReactionToken *rxn_ptr;
  Master *master_ptr;

  if (Equal(coef_in, (LDBLE)0.0, (LDBLE)TOLERANCE))
    return true;
  
	// Gamma term for d molality of species 
	// Note dg includes molality as a factor 
  row = row * (count_unknowns + 1);
	Species *s = gd->species_list[k];

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.buildjacobiansums_f, "StoreDN-1) %d %s\n",
																		row,
																		s->name.CharPtr());
#endif

  if (s->type != SURF && s != gd->s_h2o && gamma_source != NULL)
	{
		if (gamma_source != NULL)
		{
			StoreJacob (gamma_source, &arr[row + md->mu_unknown->number], (LDBLE)-1.0 * coef_in);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.buildjacobiansums_f, "StoreDN-2) %d %20.20e\n",
																				row + md->mu_unknown->number,
																				coef_in);
#endif
		}
	}

	// Mass of water factor
  if (md->mass_oxygen_unknown != NULL && s->type != EX && s->type != SURF)
	{
    StoreJacob (source, &arr[row + md->mass_oxygen_unknown->number], coef_in);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.buildjacobiansums_f, "StoreDN-3) %d %20.20e\n",
																			row + md->mass_oxygen_unknown->number,
																			coef_in);
#endif
	}

	if (s == gd->s_h2o)
    return true;

	for (j = 1; j < s->rxn_x->token_list->Count(); j++)
  {
		rxn_ptr = s->rxn_x->token_list->Element(j);

    if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
      master_ptr = rxn_ptr->s->secondary;
    else
      master_ptr = rxn_ptr->s->primary;

    if (master_ptr->u == NULL)
      continue;

    col = master_ptr->u->number;
    coef = coef_in * rxn_ptr->coef;
    StoreJacob (source, &arr[row + col], coef);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.buildjacobiansums_f, "StoreDN-4) %d %s\n",
																			row + col,
																			coef);
#endif
  }

	return true;
}

bool ModelEngine::StoreJacob(LDBLE * source, LDBLE * target, LDBLE coef)
{
	// Adds a new item to either sum_jacob1 or sum_jacob2
	// If coef is 1.0, adds to sum_jacob1, which does not require a multiply
	// Otherwise, adds to sum_jacob2, which allows multiply by coef

 STCoef *n_i;

  if (Equal(coef, (LDBLE)1.0, (LDBLE)TOLERANCE))
  {
		n_i = md->sum_jacob1->AddNew();

		n_i->source = source;
		n_i->target = target;
  }
  else
  {
		n_i = md->sum_jacob2->AddNew();

		n_i->source = source;
		n_i->target = target;
		n_i->coef = coef;
  }

	return true;
}

bool ModelEngine::IsSpecial(Species *spec)
{
	//
	// Checks to see if a species is composed of only H, O, and e-
	// Returns TRUE if true
	//         FALSE if not
	//

  ReactionToken *token_ptr;

	for (int i = 1; i < spec->rxn_s->token_list->Count(); i++)
  {
		token_ptr = spec->rxn_s->token_list->Element(i);

    if (token_ptr->s != gd->s_hplus && token_ptr->s != gd->s_h2o && token_ptr->s != gd->s_eminus)
			return false;
  }

  return true;
}

bool ModelEngine::StoreJacob0(int row, int column, LDBLE coef)
{
	//
	// Stores in list a constant coef which will be added into jacobian array
	//

	STCoef *n_i = md->sum_jacob0->AddNew();

	n_i->target = &arr[row * (count_unknowns + 1) + column];
	n_i->coef = coef;

  return true;
}

bool ModelEngine::ChangeHydrogenInEOSList(LDBLE charge)
{
  int j;
  int found_h, found_o;

  LDBLE coef_h, coef_o, coef;

	found_h = -1;
  found_o = -1;

	coef_h = 0.0;
  coef_o = 0.0;

	CombineElements();

	ElementOfSpecies *eos_p;
  for (j = 0; j < eos_list->Count(); j++)
  {
		eos_p = (*eos_list)[j];

    if (eos_p->e->name == "H")
    {
      found_h = j;
      coef_h = eos_p->coef;
    }
    else if (eos_p->e->name == "O")
    {
      found_o = j;
      coef_o = eos_p->coef;
    }
  }

  coef = coef_h - 2 * coef_o - charge;

  if (found_h < 0 && found_o < 0)
    return true;

  if (found_h >= 0 && found_o < 0)
    return true;

  if (found_h < 0 && found_o >= 0)
  {
		ElementOfSpecies *neos = eos_list->AddNew();

		neos->name = gd->s_hplus->primary->e->name;
		neos->e = gd->s_hplus->primary->e;
		neos->coef = coef;

		CombineElements();

		return true;
  }

	eos_p = (*eos_list)[found_h];
  eos_p->coef = coef;

  return true;
}

bool ModelEngine::StoreSumDeltas(LDBLE * source, LDBLE * target, LDBLE coef)
{
	//
	// List sum_delta is summed to determine the change in the mass of 
	// each element due to mass transfers of minerals, changes show up
	// in x[i]->delta. These may be multiplied by a factor under some
	// situations where the entire calculated step is not taken
	//

	STCoef * nstc = md->sum_delta->AddNew();

  nstc->source = source;
  nstc->target = target;
  nstc->coef = coef;

	return true;
}

CONVERGE_RESULT ModelEngine::ExecuteModel()
{
	//
  // ExecuteModel is called after the equations have been set up by prep
  // and initial guesses have been made in set.
  //
  // Here is the outline of the calculation sequence:
  // residuals--residuals are calculated, if small we are done
  // sum_jacobian--jacobian is calculated
  // ineq--inequality solver is called
  // reset--estimates of unknowns revised, if changes are small solution
  // has been found, usually convergence is found in residuals.
  // gammas--new activity coefficients
  // molalities--calculate molalities
  // mb_sums--calculate mass-balance sums
  // mb_gases--decide if gas_phase exists
  // mb_s_s--decide if solid_solutions exists
  // switch_bases--check to see if new basis species is needed
  // reprep--rewrite equations with new basis species if needed
  // revise_guesses--revise unknowns to get initial mole balance
  // check_residuals--check convergence one last time
  // sum_species--calculate sums of elements from species concentrations
  //
  //	  An additional pass through may be needed if unstable phases still exist
  //		 in the phase assemblage.
 
  int kode, return_kode;
  CONVERGE_RESULT r;
  int count_infeasible, count_basis_change;
  bool mass_water_switch_save, result;
	String filename;

  mass_water_switch_save = md->mass_water_switch;
	
  if (!mass_water_switch_save && md->delay_mass_water) //delay_mass_water will be always FALSE (used to debug in original version)
    md->mass_water_switch = true;

  md->pe_step_size_now = md->pe_step_size;
  md->step_size_now  = md->step_size;

  md->iterations = 0;
  count_basis_change = 0;
	count_infeasible = 0;

  md->stop_program = false;
  md->remove_unstable_phases = false;

  for (int i = 0;; i++)
  {
    if (!MBGases()) 
			return ERROR_CR;
    if (!MBSS()) 
			return ERROR_CR;

    kode = 1;


		while ((r = Residuals()) != CONVERGED_CR || md->remove_unstable_phases)
    {
			if (r == ERROR_CR) 
				return ERROR_CR;

			md->iterations++;

      if (md->iterations > md->itmax)
			{
				// Iterations exceeded
				md->stop_program = true;
				break;
			}

			if (!JacobianSums ()) 
				return ERROR_CR;

			if (!NumericalJacobian ()) 
				return ERROR_CR;

			// Full matrix with pure phases
			if (r == OK_CR || md->remove_unstable_phases)
      {
				return_kode = Ineq(kode);
				
				if (return_kode != OK)
				  count_infeasible++;

				if (return_kode == 2)
				{
					Ineq(0);
				}

				Reset();
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.gammas_f, "called_by_ExecuteModel-1)\n");
#endif
			result = Gammas(md->mu_x);
      if (!result) 
				return ERROR_CR;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.molalities_f, "called_by_ExecuteModel-1)\n");
#endif
      if (!Molalities(false))
			{
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reviseguesses_f, "called_by_ExecuteModel-1)\n");
#endif
				if (!ReviseGuesses ()) 
					return ERROR_CR;
			}

			if (md->use.sur_p != NULL && md->use.sur_p->dl_type != NO_DL && md->use.sur_p->related_phases)
			{
				if (!InitialSurfaceWater()) 
					return ERROR_CR;
			}
      
			if (!MBSums())
				return ERROR_CR;

			if (!MBGases())
				return ERROR_CR;

			if (!MBSS())
				return ERROR_CR;

			// Switch bases if necessary
      if (SwitchBases())
      {
				count_basis_change++;

				if (!Reprep ())
					return ERROR_CR;

				result = Gammas(md->mu_x);
				if (!result)
					return ERROR_CR;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.molalities_f, "called_by_ExecuteModel-2)\n");
#endif
				if (!Molalities (true))
					return ERROR_CR;

				if (md->use.sur_p != NULL && md->use.sur_p->dl_type != NO_DL && md->use.sur_p->related_phases)
					if (!InitialSurfaceWater ()) 
						return ERROR_CR;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reviseguesses_f, "called_by_ExecuteModel-2)\n");
#endif
				if (!ReviseGuesses ())
					return ERROR_CR;

				if (!MBSums())
					return ERROR_CR;
				
				if (!MBGases())
					return ERROR_CR;
				
				if (!MBSS())
					return ERROR_CR;
			}

			if (md->stop_program)
				break;
    }

    if (md->stop_program)
      break;

    if (!CheckResiduals())
    {
      md->stop_program = true;
      break;
    }

    if (!md->remove_unstable_phases && !mass_water_switch_save &&	md->mass_water_switch)
    {
      md->mass_water_switch = false;
      continue;
    }

    if (!md->remove_unstable_phases)
      break;
  }

	//d.PrintLAToFile("model-26", count_unknowns, &gd->x);
	//d.PrintSpeciesInfoToFile("ExecuteModel-26", md->species_info_list, md, gd);

	if (md->stop_program)
		return ERROR_CR;

  return OK_CR;
}

bool ModelEngine::MBGases()
{
  md->gas_in = false;

	if (md->gas_unknown == NULL || md->use.gas_p == NULL)
    return true;

  if (md->use.gas_p->type == PRESSURE && (md->gas_unknown->f > (md->use.gas_p->total_p + 1e-7) || md->gas_unknown->moles > MIN_TOTAL))
      md->gas_in = true;

	return true;
}

bool ModelEngine::MBSS()
{
  int i, j, k;
  LDBLE lp, log10_iap, total_moles;
  LDBLE iapc, iapb, kc, kb, lc, lb, xcaq, xbaq, xb, xc;
  LDBLE sigmapi_aq, sigmapi_solid;
  LDBLE total_p;
  SS *s_s_ptr;
  ReactionToken *rxn_ptr;
  Phase *phase_ptr;
	SSComp *ssc_p, *ssc_p0, *ssc_p1;

	// Determines whether solid solution equation is needed
  if (md->s_s_unknown == NULL || md->use.ssa_p == NULL)
    return true;

	for (i = 0; i < md->use.ssa_p->ss_list->Count(); i++)
  {
    s_s_ptr = (*md->use.ssa_p->ss_list)[i];

    total_moles = 0;

		for (j = 0; j < s_s_ptr->comps_list->Count(); j++)
    {
			ssc_p = (*s_s_ptr->comps_list)[j];
      total_moles += ssc_p->moles;
    }

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "MBSS-1) %20.20e\n", total_moles);
#endif

    if (total_moles > 1e-13)
      s_s_ptr->s_s_in = true;
    else if (s_s_ptr->a0 != 0.0 || s_s_ptr->a1 != 0.0)
    {
			ssc_p0 = (*s_s_ptr->comps_list)[0];
			ssc_p1 = (*s_s_ptr->comps_list)[1];

      //  Calculate IAPc and IAPb
			if (ssc_p0->phase->rxn_x->token_list->Count() > 0)
      {
				log10_iap = 0;
		
				//ToDo: Check if the "translation" of the original code above was done right
				//for (rxn_ptr = ss0->phase->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
				for (k = 1; k < ssc_p0->phase->rxn_x->token_list->Count(); k++)
				{
					rxn_ptr = ssc_p0->phase->rxn_x->token_list->Element(k);
					log10_iap += rxn_ptr->s->la * rxn_ptr->coef;
				}

				iapc = exp (log10_iap * md->LOG_10);
      }
      else
      {
				iapc = (LDBLE)1e-99;
      }

			if (ssc_p1->phase->rxn_x->token_list->Count() > 0)
      {
				log10_iap = 0;

				//ToDo: Check if the "translation" of the original code above was done right
				//for (rxn_ptr = ss1->phase->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
				for (k = 1; k < ssc_p1->phase->rxn_x->token_list->Count(); k++)
				{
					rxn_ptr = ssc_p1->phase->rxn_x->token_list->Element(k);
					log10_iap += rxn_ptr->s->la * rxn_ptr->coef;
				}

				iapb = exp (log10_iap * md->LOG_10);
      }
      else
      {
				iapb = (LDBLE)1e-99;
      }
      
      // Calculate sigma pi, aq
      sigmapi_aq = iapc + iapb;

      // Calculate xc,aq and xb, aq
      xcaq = iapc / (iapb + iapc);
      xbaq = iapb / (iapb + iapc);
      
      // Get Kc and Kb
      kc = exp (ssc_p0->phase->lk * md->LOG_10);
      kb = exp (ssc_p1->phase->lk * md->LOG_10);
      
			// Solve for xb
      xb = SSRoot(s_s_ptr->a0, s_s_ptr->a1, kc, kb, xcaq, xbaq);
      
			// Calculate lambdac and lambdab
      xc = 1 - xb;
      lc = exp ((s_s_ptr->a0 - s_s_ptr->a1 * (-4 * xb + 3)) * xb * xb);
      lb = exp ((s_s_ptr->a0 + s_s_ptr->a1 * (4 * xb - 1)) * xc * xc);
      
			// Calculate sigma pi, solid
      sigmapi_solid = xb * lb * kb + xc * lc * kc;
      
			// If Sigma pi, solid < sigma pi, aq, then use eqns
      if (sigmapi_solid < sigmapi_aq)
				s_s_ptr->s_s_in = true;
      else
				s_s_ptr->s_s_in = false;
    }
    else
    {
      // Calculate total mole fraction from solution activities
      total_p = 0;

      for (j = 0; j < s_s_ptr->comps_list->Count(); j++)
      {
				phase_ptr = (*s_s_ptr->comps_list)[j]->phase;

				if (phase_ptr->in)
				{
					lp = -phase_ptr->lk;
	  
					//ToDo: Check if the "translation" of the original code above was done right
					//for (rxn_ptr = s_s_ptr->comps[j].phase->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
					for (k = 1; k < phase_ptr->rxn_x->token_list->Count(); k++)
					{
						rxn_ptr = phase_ptr->rxn_x->token_list->Element(k);
						lp += rxn_ptr->s->la * rxn_ptr->coef;
					}
	  
					total_p += exp (lp * md->LOG_10);
				}
      }

      if (total_p > 1.0)
				s_s_ptr->s_s_in = true;
      else
				s_s_ptr->s_s_in = false;	    
		}

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "MBSS-2) %d\n", s_s_ptr->s_s_in);
#endif
  }

	Unknown *u;
	for (i = md->s_s_unknown->number; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

    if (u->type != S_S_MOLES)
      break;

    u->s_s_in = u->s_s->s_s_in;
  }

  return true;
}

CONVERGE_RESULT ModelEngine::Residuals()
{
	//
	// Calculates residuals for all equations
	//

  int i, j;
  bool converge;

  LDBLE toler,
        sum_residual,
				sinh_constant,
				sum, 
				sum1,
				*res,
				g_moles;

  Master *master_ptr, 
		     *master_ptr1, 
				 *master_ptr2,
				 *m_p;

	Unknown *x, *u2;
	Species *s_p;
  
	LDBLE sigmaddl, 
		    negfpsirt;

  sum_residual = 0.0;
	sigmaddl = 0;
  sum = 0;
  converge = true;
  toler = md->convergence_tolerance;

	int index;

	// Calculate residuals
	for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];
		res = &residual[i];

#ifdef DEBUG_MOHID
		if (d.debug_status)
		{
			if(md->mass_oxygen_unknown != NULL)
			{
				fprintf(d.set_f, "Residuals-1) %d %d %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %d\n", 
													md->mass_water_switch,
													x->type,
													x->moles,
													x->f,
													md->LOG_10,
													toler,
													md->mu_x,
													md->mass_water_aq_x,
													gd->s_h2o->la,
													md->mass_oxygen_unknown->moles,
													md->mass_oxygen_unknown->f,
													md->iterations);
			}
			else
			{
				fprintf(d.set_f, "Residuals-1) %d %d %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %d\n", 
													md->mass_water_switch,
													x->type,
													x->moles,
													x->f,
													md->LOG_10,
													toler,
													md->mu_x,
													md->mass_water_aq_x,
													gd->s_h2o->la,
													md->iterations);
			}
		}
#endif

    if (x->type == MB)
    {
      *res = x->moles - x->f;

			if (fabs(*res) > toler * x->moles && x->moles > MIN_TOTAL)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-2) %20.20e\n", *res);
#endif
    }
    else if (x->type == ALK)
    {
      *res = x->moles - x->f;

      if (fabs(*res) > toler * x->moles)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-3) %20.20e\n", *res);
#endif
		}
    else if (x->type == SOLUTION_PHASE_BOUNDARY)
    {
      *res = x->f * md->LOG_10;

      if (fabs(*res) > toler)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-4) %20.20e\n", *res);
#endif
    }
    else if (x->type == CB)
    {
      *res = -x->f;

      if (md->ph_unknown == md->charge_balance_unknown) 
				*res += x->moles;
      
      if (fabs(*res) >= toler * md->mu_x * md->mass_water_aq_x)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-5) %20.20e\n", *res);
#endif
    }
    else if (x->type == MU) //&& pitzer_model == FALSE
    {
      *res = md->mass_water_aq_x * md->mu_x - (LDBLE)0.5 * x->f;

      if (fabs(*res) > toler * md->mu_x * md->mass_water_aq_x)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-6) %20.20e\n", *res);
#endif
    }
    else if (x->type == AH2O) //&& pitzer_model == FALSE
    {
      *res = md->mass_water_aq_x * exp (gd->s_h2o->la * md->LOG_10) - md->mass_water_aq_x + (LDBLE)0.017 * x->f;

			if (fabs(*res) > toler)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-7) %20.20e\n", *res);
#endif
    }
    else if (x->type == MH) // && (pitzer_model == FALSE || pitzer_pe == TRUE))
    {
      *res = x->moles - x->f;

      if (md->mass_water_switch)
				*res -= 2 * (md->mass_oxygen_unknown->moles - md->mass_oxygen_unknown->f);

      if (fabs(*res) > toler * (x->moles + 2 * md->mass_oxygen_unknown->moles))
				converge = FALSE;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-8) %20.20e\n", *res);
#endif
    }
    else if (x->type == MH2O)
    {
      if (md->mass_water_switch)
				continue;
	     
			*res = x->moles - x->f;

      if (fabs(*res) > 0.01 * toler * x->moles)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-9) %20.20e\n", *res);
#endif
    }
    else if (x->type == PP)
    {
      *res = x->f * md->LOG_10;

			if (x->pure_phase->add_formula.IsEmpty())
      {
				if (x->dissolve_only)
				{
					if ((*res > toler && x->moles > 0.0) || (*res < -toler && (x->pure_phase->initial_moles - x->moles) > 0))
						converge = false;
				}
				else
				{
					if (*res < -toler || md->iterations < 1)
						converge = false;
				}

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.set_f, "Residuals-10) %20.20e %20.20e\n", *res, x->pure_phase->initial_moles);
#endif
      }
      else
      {
				if (*res < -toler || md->iterations < 1)
					converge = FALSE;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.set_f, "Residuals-11) %20.20e\n", *res);
#endif
      }
    }
    else if (x->type == GAS_MOLES)
    {
      *res = x->gas_phase->total_p - x->f;

      if (fabs(*res) > toler && md->gas_in)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-12) %20.20e %20.20e %d\n", 
													*res,
													x->gas_phase->total_p,
													md->gas_in);
#endif
    }
    else if (x->type == S_S_MOLES)
    {
      *res = x->f * md->LOG_10;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-13) %20.20e\n", *res);
#endif

      if (x->moles <= MIN_TOTAL_SS && md->iterations > 2)
				continue;

      if (fabs(*res) > toler && x->s_s_in)
				converge = false;
    }
    else if (x->type == EXCH)
    {
      *res = x->moles - x->f;

      if (x->moles <= MIN_RELATED_SURFACE)
      {
				if (fabs(*res) > toler)
					converge = false;
      }
      else if (fabs(*res) > toler * x->moles)
				converge = FALSE;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-14) %20.20e\n", *res);
#endif
    }
    else if (x->type == SURFACE)
    {
      *res = x->moles - x->f;

      if (x->moles <= MIN_RELATED_SURFACE)
      {
				if (fabs(*res) > toler)
					converge = false;
      }
      else if (fabs(*res) > toler * x->moles)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-15) %20.20e\n", *res);
#endif
    }
    else if (x->type == PITZER_GAMMA)
    {
      if (!md->full_pitzer)
				continue;

      *res = x->s->lg - x->s->lg_pitzer;

      if (fabs(*res) > toler)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-16) %20.20e\n", *res);
#endif
    }
    else if (x->type == SURFACE_CB && md->use.sur_p->type == DDL)
    {
      sinh_constant = sqrt((LDBLE)8 * (LDBLE)EPSILON * (LDBLE)EPSILON_ZERO * ((LDBLE)R_KJ_DEG_MOL * (LDBLE)1000) * md->tk_x * (LDBLE)1000);

			if (x->surface_charge->grams == 0)
				*res = 0.0;
      else if (md->dl_type_x != NO_DL)
				*res = -x->f;
      else
				*res = sinh_constant * sqrt (md->mu_x) * sinh ((*x->master)[0]->s->la * md->LOG_10) - x->f * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);

      if (x->surface_charge->grams > MIN_RELATED_SURFACE && fabs(*res) > toler)
				converge = false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-17) %20.20e\n", *res);
#endif
    }
    else if (x->type == SURFACE_CB && md->use.sur_p->type == CD_MUSIC)
    {
      if (x->surface_charge->grams == 0)
				*res = 0.0;
      else
      {
				// sum is in moles of charge
				master_ptr  = SurfaceGetPSIMaster(x->surface_charge->name, SURF_PSI);
				master_ptr1 = SurfaceGetPSIMaster(x->surface_charge->name, SURF_PSI1);
				master_ptr2 = SurfaceGetPSIMaster(x->surface_charge->name, SURF_PSI2);

				x->surface_charge->psi  = -(master_ptr->s->la * md->LOG_10) * (LDBLE)R_KJ_DEG_MOL * md->tk_x / (LDBLE)F_KJ_V_EQ;
				x->surface_charge->psi1 = -(master_ptr1->s->la * md->LOG_10) * (LDBLE)R_KJ_DEG_MOL * md->tk_x / (LDBLE)F_KJ_V_EQ;
				x->surface_charge->psi2 = -(master_ptr2->s->la * md->LOG_10) * (LDBLE)R_KJ_DEG_MOL * md->tk_x / (LDBLE)F_KJ_V_EQ;

				sum = 0;

				for (j = 0; j < x->comp_unknowns->Count(); j++)
				{
					u2 = (*x->comp_unknowns)[j];
					m_p = (*u2->master)[0];
					sum += u2->moles * m_p->s->z;
				}

				x->surface_charge->sigma0 = (x->f + sum) * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);

				// f is in moles
				// eqns A-3
				*res = x->surface_charge->sigma0 - x->surface_charge->capacitance[0] * (x->surface_charge->psi - x->surface_charge->psi1);
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-18) %20.20e\n", *res);
#endif
      if (x->surface_charge->grams > MIN_RELATED_SURFACE && fabs(*res) > toler)
				converge = false;
    }
    else if (x->type == SURFACE_CB1)
    {
      if (x->surface_charge->grams == 0)
				*res = 0.0;
      else
      {
				// eqns A-4
				x->surface_charge->sigma1 =	x->f * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);
				*res = (x->surface_charge->sigma0 + x->surface_charge->sigma1) - x->surface_charge->capacitance[1] * (x->surface_charge->psi1 - x->surface_charge->psi2);
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-19) %20.20e\n", *res);
#endif
			if (x->surface_charge->grams > MIN_RELATED_SURFACE && fabs (*res) > toler)
				converge = false;
    }
    else if (x->type == SURFACE_CB2)
    {
      if (x->surface_charge->grams == 0)
				*res = 0.0;
      else if (md->dl_type_x != NO_DL)
      {
				sum = 0;
				sum1 = 0;

				for (j = 0; j < md->s_x->Count(); j++)
				{
					s_p = (*md->s_x)[j];

					if (s_p->type == SURF)
						sum += Under(s_p->lm) * s_p->dz[2];

					if (s_p->type < H2O)
					{
						g_moles = (*s_p->diff_layer)[0]->g_moles;
						sum1 += s_p->z * g_moles;
					}
				}

				x->surface_charge->sigma2 = sum * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);
				x->surface_charge->sigmaddl = (x->f - sum) * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);

				*res = x->f + (x->surface_charge->sigma0 + x->surface_charge->sigma1) * (x->surface_charge->specific_area * x->surface_charge->grams) / (LDBLE)F_C_MOL;
      }
      else
      {
				// eqns A-6 and A-7
				sinh_constant = sqrt ((LDBLE)8 * (LDBLE)EPSILON * (LDBLE)EPSILON_ZERO * ((LDBLE)R_KJ_DEG_MOL * (LDBLE)1000) * md->tk_x * (LDBLE)1000);
				master_ptr2 = SurfaceGetPSIMaster(x->surface_charge->name, SURF_PSI2);
				negfpsirt = master_ptr2->s->la * md->LOG_10;
				
				sum = 0;
				sum1 = 0;

				for (j = 0; j < md->s_x->Count(); j++)
				{
					s_p = (*md->s_x)[j];

					if (s_p->type < H2O)
					{
						sum += Under(s_p->lm) * (exp(s_p->z * negfpsirt) - 1);
						sum1 += Under(s_p->lm) * s_p->z;
					}
				}

				// add fictitious monovalent ion that balances charge
				sum += sum1 * (exp(-sum1 / fabs (sum1) * negfpsirt) - 1);

				if (sum < 0)
				{
					sum = -sum;
					converge = false;
				}

				x->surface_charge->sigma2 = x->f * (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);

				if ((negfpsirt) < 0)
					sigmaddl = (LDBLE)-0.5 * sinh_constant * sqrt(sum);
				else
					sigmaddl = (LDBLE)0.5 * sinh_constant * sqrt(sum);

				x->surface_charge->sigmaddl = sigmaddl;
	
				*res = (x->surface_charge->sigma0 + x->surface_charge->sigma1 + x->surface_charge->sigma2) + sigmaddl;
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "Residuals-20) %20.20e\n", *res);
#endif
      if (x->surface_charge->grams > MIN_RELATED_SURFACE && fabs(*res) > toler)
				converge = false;
    }

		index = (i + 1) * (count_unknowns + 1) - 1;
    arr[index] = *res;
    sum_residual += fabs(*res);	
	}

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.set_f, "Residuals-21) %20.20e\n", sum_residual);
#endif
			
	if (converge)
		return (CONVERGED_CR);

	return (OK_CR);
}

bool ModelEngine::KTemp(LDBLE tempc)
{
	//
	// Calculates log k's for all species and pure_phases
	//

  int i;
  LDBLE tempk;

  tempk = tempc + (LDBLE)273.15;

	// Calculate log k for all aqueous species
	Species *s;
  for (i = 0; i < md->s_x->Count(); i++)
	{
		s = (*md->s_x)[i];
    s->lk = KCalc(s->rxn_x->logk, tempk);

#ifdef DEBUG_MOHID
		if (d.debug_status)
		{
			fprintf(d.set_f, "KTemp-1) %s %d ", s->name.CharPtr(), i);
			for (int debug_i = 0; debug_i < 7; debug_i++)
				fprintf(d.set_f, "%20.20e ", s->rxn_x->logk[debug_i]);
			fprintf(d.set_f, "%20.20e\n", tempk);
		}
#endif
	}

	Phase *p;
	for (i = 0; i < gd->phase_list.Count(); i++)
	{
		p = gd->phase_list[i];

    if (p->in)
      p->lk = KCalc(p->rxn_x->logk, tempk);
	}

	// Calculate miscibility gaps for solid solutions
	SS *ss_p;
	if (md->use.ssa_p != NULL)
  {
		for (i = 0; i < md->use.ssa_p->ss_list->Count(); i++)
		{
			ss_p = (*md->use.ssa_p->ss_list)[i];

      if (fabs(tempk - ss_p->tk) > 0.01)
				SSPrep(tempk, ss_p);
		}
  }

  return true;
}

bool ModelEngine::Set(bool initial)
{
	//
	// Sets initial guesses for unknowns if initial == TRUE
	// Revises guesses whether initial is true or not
	//

  int i;
  
  md->iterations = -1;
	Solution *sol_p = md->use.sol_p;

	// Set initial log concentrations to zero
	Species *s;
  for (i = 0; i < md->s_x->Count(); i++)
	{	
		s = (*md->s_x)[i];

    s->lm = LOG_ZERO_MOLALITY;
    s->lg = 0.0;
  }

	// Set master species activities
  md->tc_x = sol_p->tc;
  md->tk_x = md->tc_x + (LDBLE)273.15;

	// H+, e-, H2O
  md->mass_water_aq_x = sol_p->mass_water;
  md->mu_x = sol_p->mu;
	//d.PrintSpeciesInfoToFile("Set", md->species_info_list, md, gd);
  gd->s_h2o->moles = md->mass_water_aq_x / md->gfw_water;
  gd->s_h2o->la = log10(sol_p->ah2o);
  gd->s_hplus->la = -sol_p->ph;
  gd->s_hplus->lm = gd->s_hplus->la;
  gd->s_hplus->moles = exp(gd->s_hplus->lm * md->LOG_10) * md->mass_water_aq_x;
  gd->s_eminus->la = -sol_p->solution_pe;
  
#ifdef DEBUG_MOHID
	if (d.debug_status)
	{
		fprintf(d.set_f, "%i %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e\n", 
											initial, md->tc_x, md->tk_x, md->mass_water_aq_x, md->mu_x, gd->s_h2o->moles, gd->s_h2o->la, 
											gd->s_hplus->la, gd->s_hplus->lm, gd->s_hplus->moles, gd->s_eminus->la);
	}
#endif 

	if (initial)
    InitialGuesses();

  if (md->dl_type_x != NO_DL)
    InitialSurfaceWater ();
   
#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.reviseguesses_f, "called_by_set: %d\n", d.reviseguesses_count++);
#endif

	ReviseGuesses ();

  return true;
}

bool ModelEngine::CheckResiduals()
{
	//
	// Checks for convergence of all equations, prints any nonconvergence
	// Sets variable remove_unstable_phases if a phase is present,
	// but undersaturated (i.e. aragonite in calcite-saturated solution).

  int i;
  LDBLE epsilon, r;
	bool return_value;

  epsilon = md->convergence_tolerance;

	return_value = true;

  if (md->stop_program)
		return false;

	Unknown *x;
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];
		r = fabs((LDBLE) residual[i]);

		switch(x->type)
		{
		case MB:
		case ALK:
      if (r >= (epsilon * x->moles) && x->moles > MIN_TOTAL)
			{
				sprintf (error_string, "%20s has not converged. Total: %e\tCalculated: %e\tResidual: %e\n", 
																x->name, (double) x->moles, (double) x->f, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);

				if (x->type == ALK)
				{
					if(error != NULL) error->SaveMessage("Is non-carbonate alkalinity greater than total alkalinity?\n");
				}

				return_value = false;
			}
			break;

		case SOLUTION_PHASE_BOUNDARY:
      if (r >= epsilon)
      {
				sprintf (error_string, "%20s solution phase boundary has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
      }
			break;
		case AH2O:
      if (r >= epsilon)
			{
				sprintf (error_string, "%20s Activity of water has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case CB:
      if (r >= epsilon * md->mu_x * md->mass_water_aq_x)
      {
				sprintf (error_string, "%20s Charge balance has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
      }
			break;
		case MU:
      if (r >= (epsilon * md->mu_x * md->mass_water_aq_x))
			{
				sprintf (error_string, "%20s Ionic strength has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case MH:
      if (r > (epsilon * (x->moles + 2 * md->mass_oxygen_unknown->moles)))
			{
				sprintf (error_string, "%20s Mass of hydrogen has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);

				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case MH2O:
      if (md->mass_water_switch)
				continue;

      if (r >= (0.01 * epsilon * x->moles))
			{
				sprintf (error_string, "%20s Mass of oxygen has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case PP:
			if (x->pure_phase->add_formula.IsEmpty())
      {
				if (x->dissolve_only)
				{
					if ((residual[i] > epsilon && x->moles > 0.0)	|| (residual[i] < -epsilon	&& (x->pure_phase->initial_moles - x->moles) > 0))
					{
						sprintf (error_string, "%20s Dissolve_only pure phase has not converged. \tResidual: %e\n",
																		x->name, (double) residual[i]);
						if(error != NULL) error->SaveMessage(error_string);
					}
				}
				else 
				{
					if (residual[i] >= epsilon && x->moles > 0.0)
					{
						md->remove_unstable_phases = true;
						sprintf (error_string, "%20s Pure phase has not converged. \tResidual: %e\n",
																		x->name, (double) residual[i]);
						if(error != NULL) error->SaveMessage(error_string);
					}
					else if (residual[i] <= -epsilon)
					{
						sprintf (error_string, "%20s Pure phase has not converged. \tResidual: %e\n", 
																		x->name, (double) residual[i]);
						if(error != NULL) error->SaveMessage(error_string);
					}
				}
      }
      else
			{
				if (r >= epsilon && x->moles > 0.0)
				{
					sprintf (error_string, "%s, Pure phase with add formula has not converged.\n\t SI may be a local minimum.\tResidual: %e\n", 
																	x->name, (double) residual[i]);
					if(error != NULL) error->SaveMessage(error_string);
				}	
			}
			break;
		case EXCH:
			if ((x->moles <= MIN_RELATED_SURFACE && r > epsilon) || (x->moles > MIN_RELATED_SURFACE && (r > epsilon * x->moles)))
			{
				sprintf (error_string, "%20s Exchanger mass balance has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case SURFACE:
			if ((x->moles <= MIN_RELATED_SURFACE && r > epsilon) || (x->moles > MIN_RELATED_SURFACE && (r > epsilon * x->moles)))
			{
				sprintf (error_string, "%20s Surface mass balance has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
      if (x->surface_charge->grams > MIN_RELATED_SURFACE && r > epsilon)
			{
				sprintf (error_string, "%20s Surface charge/potential has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case GAS_MOLES:
      if (!md->gas_in)
				continue;

      if (residual[i] >= epsilon || residual[i] <= -epsilon)
			{
				sprintf (error_string, "%20s Total moles in gas phase has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case PITZER_GAMMA:
      if (r > epsilon)
			{
				sprintf (error_string, "%20s log gamma not converged.\tResidual: %e\n",
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		case S_S_MOLES:
      if (!x->s_s_in)
				continue;

      if (x->moles <= MIN_TOTAL_SS)
				continue;

      if (residual[i] >= epsilon || residual[i] <= -epsilon)
			{
				sprintf (error_string, "%20s Total moles in solid solution has not converged. \tResidual: %e\n", 
																x->name, (double) residual[i]);
				if(error != NULL) error->SaveMessage(error_string);
			}
			break;
		}
  }

  if (md->remove_unstable_phases)
  {
    sprintf (error_string, "%20sRemoving unstable phases, iteration %d.", " ", md->iterations);
		if(error != NULL) error->SaveMessage(error_string);
  }

	return return_value;
}

bool ModelEngine::SumSpecies()
{
	//
	// Calculates total alk, total carbon, total co2, electrical balance,
	// total hydrogen, and total oxygen.
	// Sorts species for summing and printing based on valence state and
	// concentrations.
	//
	// Sums total valence states and stores in master[i]->total.
	//

  int i, j;
  Master *master_ptr;

	// Set global variables
	md->ph_x = -gd->s_hplus->la;
  md->solution_pe_x = -gd->s_eminus->la;
  md->ah2o_x = exp(gd->s_h2o->la * md->LOG_10);
  md->density_x = 1.0;

  if (gd->s_o2 != NULL)
    gd->s_o2->moles = Under(gd->s_o2->lm) * md->mass_water_aq_x;

  if (gd->s_h2 != NULL)
    gd->s_h2->moles = Under(gd->s_h2->lm) * md->mass_water_aq_x;

	// Calculate sums
  md->total_alkalinity = 0.0;
  md->total_carbon = 0.0;
  md->total_co2 = 0.0;
  md->cb_x = 0.0;
  md->total_ions_x = 0.0;
  md->total_o_x = 0.0;
  md->total_h_x = 0.0;

	Species *sx;
  for (i = 0; i < md->s_x->Count(); i++)
  {
		sx = (*md->s_x)[i];

    if (sx->type == EX || sx->type == SURF)
      continue;

    md->cb_x += sx->z * sx->moles;
    md->total_ions_x += fabs(sx->z * sx->moles);
    md->total_alkalinity += sx->alk * sx->moles;
    md->total_carbon += sx->carbon * sx->moles;
    md->total_co2 += sx->co2 * sx->moles;

    md->total_h_x += sx->h * sx->moles;
    md->total_o_x += sx->o * sx->moles;

		if (md->use.sur_p != NULL)
    {
      if (md->use.sur_p->debye_lengths > 0 && md->state >= REACTION && sx->type == H2O)
      {
				md->total_h_x -= 2 * md->mass_water_surfaces_x / md->gfw_water;
				md->total_o_x -= md->mass_water_surfaces_x / md->gfw_water;
      }
    }
  }

	// Sum valence states, put in master->total
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		master_ptr = gd->master_list[i];
    master_ptr->total = 0.0;
    master_ptr->total_primary = 0.0;
  }

	SpeciesInfo *si;
  for (i = 0; i < md->species_info_list->Count(); i++)
  {
		si = (*md->species_info_list)[i];

    if (si->master_s->secondary != NULL)
      master_ptr = si->master_s->secondary;
    else
      master_ptr = si->master_s->primary;

		master_ptr->total += si->s->moles * si->coef;
  }

	// Calculate mass-balance sums
	Unknown *u;
  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

    if (u->type == MB || u->type == SOLUTION_PHASE_BOUNDARY || u->type == EXCH ||
				u->type == SURFACE ||	(u->type == CB && u != md->ph_unknown && u != md->pe_unknown))
    {
      u->sum = 0.0;

			for (j = 0; j < u->master->Count(); j++)
			{
				master_ptr = (*u->master)[j];
				u->sum += master_ptr->total;
			}
    }
    else if (u->type == ALK)
      u->sum = md->total_co2;
  }

	// Calculate total element concentrations
  for (i = 0; i < gd->master_list.Count(); i++)
	{
		master_ptr = gd->master_list[i];
    master_ptr->e->primary->total_primary += master_ptr->total;
	}

	// Calculate isotope ratios
  //CalculateValues (); //Probably isn't used until isotopes are included on program

  return true;
}

bool ModelEngine::NumericalJacobian()
{
	LDBLE *base;

  LDBLE d_, d1, d2;
  int i, j;

	if (md->use.sur_p == NULL || md->use.sur_p->type != CD_MUSIC)
    return true;

	//d.PrintArrayToFile(arr, arr_capacity, "numerical_jacobian-0");
	base = NULL;

	try
	{
		md->calculating_deriv = true;

		Gammas(md->mu_x);
		//d.PrintGammasToFile("NumericalJacobian-1", md->s_x);
		Molalities(true);
		MBSums();
		Residuals();
		//d.PrintArrayToFile(arr, arr_capacity, "NumericalJacobian-1");

		// Clear array, note residuals are in array[i, count_unknowns+1]
		for (i = 0; i < count_unknowns; i++)
			arr[i] = 0.0;

		for (i = 1; i < count_unknowns; i++)
			memcpy ((void *) &arr[i * (count_unknowns + 1)], (void *) &arr[0], (size_t) count_unknowns * sizeof (LDBLE));

		base = new LDBLE [count_unknowns];

		for (i = 0; i < count_unknowns; i++)
			base[i] = residual[i];

		d_  = (LDBLE)1e-6;
		d1 = d_ * log((LDBLE)10.0);
		d2 = 0;

		Unknown *x;
		Master *m;
		for (i = 0; i < count_unknowns; i++)
		{
			x = (*unknown_list)[i];
			m = (*x->master)[0];

			switch (x->type)
			{
			case MB:
			case ALK:
			case CB:
			case SOLUTION_PHASE_BOUNDARY:
			case EXCH:
			case SURFACE:
			case SURFACE_CB:
			case SURFACE_CB1:
			case SURFACE_CB2:
				m->s->la += d_;
				d2 = d1;
				break;
			case MH:
				gd->s_eminus->la += d_;
				d2 = d1;
				break;
			case AH2O:
				m->s->la += d_;
				d2 = d1;
				break;
			case PITZER_GAMMA:
				x->s->lg += d_;
				d2 = d_;
				break;
			case MH2O:
				md->mass_water_aq_x *= ((LDBLE)1.0 + d_);
				m->s->moles = md->mass_water_aq_x / md->gfw_water;
				d2 = log ((LDBLE)1.0 + d_);
				break;
			case MU:
				d2 = d_ * md->mu_x;
				md->mu_x += d2;
				//d.PrintSpeciesInfoToFile("NumericalJacobian-1", md->species_info_list, md, gd);
				Gammas(md->mu_x);
				//d.PrintGammasToFile("NumericalJacobian-2", md->s_x);
				break;
			case PP:
				for (j = 0; j < count_unknowns; j++)
					delta[j] = 0.0;

				d2 = (LDBLE)-1e-8;
				delta[i] = d2;
				Reset();
				//d.PrintDeltaToFile("Reset-2", delta, count_unknowns);
				d2 = delta[i];
				break;
			case S_S_MOLES:
				if (!x->s_s_in)
					continue;

				for (j = 0; j < count_unknowns; j++)
					delta[j] = 0.0;

				d2 = -d_ * x->moles;
				d2 = (LDBLE)-.1 * x->moles;
				delta[i] = d2;
				Reset ();
				//d.PrintDeltaToFile("Reset-3", delta, count_unknowns);
				d2 = delta[i];
				break;
			case GAS_MOLES:
				if (!md->gas_in)
					continue;

				d2 = d_ * x->moles;

				if (d2 < 1e-14)
					d2 = (LDBLE)1e-14;

				x->moles += d2;
				break;
			}

			Molalities(true);
			MBSums();
			Residuals ();
			//d.PrintArrayToFile(arr, arr_capacity, "NumericalJacobian-2");

			for (j = 0; j < count_unknowns; j++)
				arr[j * (count_unknowns + 1) + i] = -(residual[j] - base[j]) / d2;

			switch (x->type)
			{
			case MB:
			case ALK:
			case CB:
			case SOLUTION_PHASE_BOUNDARY:
			case EXCH:
			case SURFACE:
			case SURFACE_CB:
			case SURFACE_CB1:
			case SURFACE_CB2:
			case AH2O:
				m->s->la -= d_;
				break;
			case MH:
				gd->s_eminus->la -= d_;

				if (arr[i * (count_unknowns + 1) + i] == 0)
					arr[i * (count_unknowns + 1) + i] = exp (gd->s_h2->lm * md->LOG_10) * 2;
				break;
			case PITZER_GAMMA:
				x->s->lg -= d_;
				break;
			case MH2O:
				md->mass_water_aq_x /= (1 + d_);
				m->s->moles = md->mass_water_aq_x / md->gfw_water;
				break;
			case MU:
				md->mu_x -= d2;
				//d.PrintSpeciesInfoToFile("NumericalJacobian-2", md->species_info_list, md, gd);
				Gammas(md->mu_x);
				//d.PrintGammasToFile("NumericalJacobian-3", md->s_x);
				break;
			case PP:
				delta[i] = -d2;
				Reset ();
				//d.PrintDeltaToFile("Reset-4", delta, count_unknowns);
				break;
			case S_S_MOLES:
				delta[i] = -d2;
				Reset ();
				//d.PrintDeltaToFile("Reset-5", delta, count_unknowns);
				break;
			case GAS_MOLES:
				x->moles -= d2;
				break;
			}
		}

		Molalities(true);
		MBSums();
		MBGases();
		MBSS();
		Residuals();
		//d.PrintArrayToFile(arr,arr_capacity, "numerical_jacobian-3");
	}
	catch(...)
	{
		delete base;
		throw;
	}

	delete base;

	md->calculating_deriv = false;

  return true;
}

bool ModelEngine::JacobianSums()
{
	//
	// Fills in jacobian array, uses arrays sum_jacob0, sum_jacob1, and	sum_jacob2.
	//

  int i, 
		  j, 
			k,
			index;

  LDBLE sinh_constant;
	STCoef *stc_p;

	// Clear array, note residuals are in array[i, count_unknowns+1]
  for (i = 0; i < count_unknowns; i++)
    arr[i] = 0.0;

  for (i = 1; i < count_unknowns; i++)
    memcpy ((void *) &arr[i * (count_unknowns + 1)], (void *) &arr[0], (size_t) count_unknowns * sizeof (LDBLE));

	// Add constant terms
	for (k = 0; k < md->sum_jacob0->Count(); k++)
	{
		stc_p = (*md->sum_jacob0)[k];
		*stc_p->target += stc_p->coef;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-1) %20.20e %20.20e\n", *stc_p->target, stc_p->coef); 
#endif
	}

	// Add terms with coefficients of 1.0
  for (k = 0; k < md->sum_jacob1->Count(); k++)
	{
		stc_p = (*md->sum_jacob1)[k];
		*stc_p->target += *stc_p->source;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-2) %20.20e %20.20e\n", *stc_p->target, *stc_p->source); 
#endif
	}

	// Add terms with coefficients != 1.0
  for (k = 0; k < md->sum_jacob2->Count(); k++)
	{
		stc_p = (*md->sum_jacob2)[k];
		*stc_p->target += *stc_p->source * stc_p->coef;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-3) %20.20e %20.20e %20.20e\n", *stc_p->target, *stc_p->source, stc_p->coef); 
#endif
	}

	// Make final adustments to jacobian array

	// Ionic strength
	if (md->mu_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
		{
			index = md->mu_unknown->number * (count_unknowns + 1) + i;
      arr[index] *= 0.5;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "JacobianSums-4) %d %20.20e\n", index, arr[index]); 
#endif
		}

		index = md->mu_unknown->number * (count_unknowns + 1) + md->mu_unknown->number;
		arr[index] -= md->mass_water_aq_x;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-5) %d %20.20e\n", index, arr[index]); 
#endif
  }


	// Activity of water
  if (md->mass_oxygen_unknown != NULL && md->mu_unknown != NULL)
	{
		index = md->mu_unknown->number * (count_unknowns + 1) + md->mass_oxygen_unknown->number;
    arr[index] -= md->mu_x * md->mass_water_aq_x;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-6) %d %20.20e\n", index, arr[index]); 
#endif
	}

  if (md->ah2o_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
		{
			index = md->ah2o_unknown->number * (count_unknowns + 1) + i;
      arr[index] *= (LDBLE)-0.017;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.set_f, "JacobianSums-7) %d %20.20e\n", index, arr[index]); 
#endif
		}

		index = md->ah2o_unknown->number * (count_unknowns + 1) + md->ah2o_unknown->number;
		arr[index] -= md->mass_water_aq_x * exp (gd->s_h2o->la * md->LOG_10);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-8) %d %20.20e\n", index, arr[index]); 
#endif
  }

  if (md->mass_oxygen_unknown != NULL && md->ah2o_unknown != NULL)
	{
		index = md->ah2o_unknown->number * (count_unknowns + 1) + md->mass_oxygen_unknown->number;
    arr[index] -= (exp(gd->s_h2o->la * md->LOG_10) - 1) * md->mass_water_aq_x;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "JacobianSums-9) %d %20.20e\n", index, arr[index]); 
#endif
	}

	// Surface charge balance
  if (md->surface_unknown != NULL && md->dl_type_x == NO_DL)
  {
    sinh_constant = sqrt ((LDBLE)8 * (LDBLE)EPSILON * (LDBLE)EPSILON_ZERO * ((LDBLE)R_KJ_DEG_MOL * (LDBLE)1000) * md->tk_x * (LDBLE)1000);

		Unknown *x;
		Master *m;
		for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];			

      if (x->type == SURFACE_CB && x->surface_charge->grams > 0)
      {
				for (j = 0; j < count_unknowns; j++)
				{
					index = x->number * (count_unknowns + 1) + j;
					arr[index] *= (LDBLE)F_C_MOL / (x->surface_charge->specific_area * x->surface_charge->grams);

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.set_f, "JacobianSums-10) %d %20.20e\n", index, arr[index]); 
#endif
				}

				m = (*x->master)[0];
				index = x->number * (count_unknowns + 1) + x->number;
				arr[index] -= sinh_constant * sqrt(md->mu_x) * cosh(m->s->la * md->LOG_10);

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.set_f, "JacobianSums-11) %d %20.20e\n", index, arr[index]); 
#endif

				if (md->mu_unknown != NULL)
				{
					index = x->number * (count_unknowns + 1) + md->mu_unknown->number;
					arr[index] -= (LDBLE)0.5 * sinh_constant / sqrt(md->mu_x) * sinh(m->s->la * md->LOG_10);

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.set_f, "JacobianSums-12) %d %20.20e\n", index, arr[index]); 
#endif
				}
      }
    }
  }

  return true;
}

int ModelEngine::Ineq(int in_kode)
{
/*
 *	Sets up equations and inequalities for Cl1.
 *	Scales columns if necessary.
 *	Eliminates equations that are not necessary because
 *		gas_phase, s_s, or phase equation is not needed
 *	Mallocs space
 *	Calls Cl1
 *	Rescales results if necessary
 */

  int i, j;
  int return_code;
  int count_rows;
  int count_optimize, count_equal;
  int max_row_count, max_column_count; //eram extern
  int k, l, m, n;
  int klmd, nklmd, n2d;
  int iter;
  LDBLE error;
  LDBLE max;
  int kode;
  LDBLE min;
	Unknown *x;
	LDBLE r, to_fill_LDBLE = 0.0;
	int to_fill_int = 0;

#ifdef DEBUG_MOHID
	if (d.debug_status)
		for(int debug_i = 0; debug_i < count_unknowns; debug_i++)
			fprintf(d.ineq_f, "ineq-0) %s %d %d %20.20e %20.20e\n", 
																(*unknown_list)[debug_i]->name.CharPtr(),
																(*unknown_list)[debug_i]->type,
																(*unknown_list)[debug_i]->number,
																(*unknown_list)[debug_i]->moles,
																(*unknown_list)[debug_i]->f);

#endif

  if (md->remove_unstable_phases)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];
			r = residual[i];

			if (x->type == PP && r > 0e-8 && x->moles > 0 && x->pure_phase->add_formula.IsEmpty() && !x->dissolve_only)
				delta[i] = x->moles;
      else
				delta[i] = 0.0;    

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-1) %d %20.20e %20.20e %20.20e\n", i, x->moles, residual[i], delta[i]);
#endif
    }

    md->remove_unstable_phases = false;

    return OK;
  }

  IneqInit(3 * count_unknowns, 3 * count_unknowns);

	NormalNewMax(count_unknowns);
	Fill(normal, 1.0, normal_max());

  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];
    max = 0.0;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.ineq_f, "Ineq-2) %d %d %20.20e %s\n", 
																i, 
																x->type,
																x->moles,
																x->name.CharPtr());
#endif

    if (x->type == MB || x->type == ALK || x->type == EXCH ||  x->type == SURFACE || x->type == SURFACE_CB || SURFACE_CB1	|| SURFACE_CB2) 
		{
			if (x->moles <= MIN_RELATED_SURFACE && (x->type == EXCH || x->type == SURFACE))
				continue;

      for (j = 0; j < count_unknowns; j++)
      {
				if (x->type == SURFACE && (*unknown_list)[j]->type == SURFACE_CB)
					continue;
				if (x->type == SURFACE_CB1 && (*unknown_list)[j]->type == SURFACE_CB2)
					continue;

				if (fabs(arr[j * (count_unknowns + 1) + i]) > max)
				{
					max = fabs(arr[j * (count_unknowns + 1) + i]);

					if (max > md->min_value)
					{
#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.ineq_f, "Ineq-3) %d %d %20.20e %20.20e %s\n", 
																				j, 
																				(*unknown_list)[j]->type,
																				max,
																				md->min_value,
																				(*unknown_list)[j]->name);
#endif

						break;
					}
				}
      }

      if (md->diagonal_scale)
      {
				if (fabs (arr[i * (count_unknowns + 1) + i]) < md->min_value)
				{
					max = fabs(arr[i * (count_unknowns + 1) + i]);

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-4) %20.20e %d %20.20e %20.20e\n", 
																			max, 
																			i * (count_unknowns + 1) + i,
																			arr[i * (count_unknowns + 1) + i],
																			md->min_value);
#endif				
				}
      }

      if (max == 0)
      {
				arr[i * (count_unknowns + 1) + i] = (LDBLE)1e-5 * x->moles;
				max = fabs ((LDBLE)1e-5 * x->moles);

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-5) %20.20e %d %20.20e\n", 
																		max, 
																		i * (count_unknowns + 1) + i,
																		arr[i * (count_unknowns + 1) + i]);
#endif
      }
    }

    if (x->type == MH)
    {    
      min = (LDBLE)1e-12;
      min = (LDBLE)MIN_TOTAL;
      arr[x->number * (count_unknowns + 1) + x->number] += min;
      
			if (fabs (arr[x->number * (count_unknowns + 1) + x->number]) < min)
				arr[x->number * (count_unknowns + 1) + x->number] = min;
      
			max = 0.0;

      for (j = 0; j < count_unknowns; j++)
      {
				if ((*unknown_list)[j]->type != MB &&
						(*unknown_list)[j]->type != SURFACE &&
						(*unknown_list)[j]->type != SURFACE_CB &&
						(*unknown_list)[j]->type != SURFACE_CB1 &&
						(*unknown_list)[j]->type != SURFACE_CB2 &&
						(*unknown_list)[j]->type != EXCH && 
						(*unknown_list)[j]->type != MH && 
						(*unknown_list)[j]->type != MH2O)
					continue;

				if (fabs (arr[j * (count_unknowns + 1) + i]) > max)
				{
					max = fabs (arr[j * (count_unknowns + 1) + i]);

					if (max > md->min_value)
						break;
				}
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-6) %20.20e %d %20.20e\n", 
																	max, 
																	x->number * (count_unknowns + 1) + x->number,
																	arr[x->number * (count_unknowns + 1) + x->number]);
#endif
    }

		if (max > 0.0 && max < md->min_value)
    {
      for (j = 0; j < count_unknowns; j++)
			{
				arr[j * (count_unknowns + 1) + i] *= md->min_value / max;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-7) %d %20.20e\n", 
																		j * (count_unknowns + 1) + i,
																		arr[j * (count_unknowns + 1) + i]);
#endif			
			}

			normal[i] = md->min_value / max;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-8) %d %20.20e\n", 
																	i,
																	normal[i]);
#endif
    }
  }

	// Allocate arrays for inequality solver
  max_row_count = 2 * count_unknowns + 2;
  max_column_count = count_unknowns + 2;

	IneqArrayNewMax(max_row_count * max_column_count);
	Fill(ineq_array, 0.0, ineq_array_max());

	BackEqNewMax(max_row_count);
	Fill(back_eq, 0, back_eq_max());

	ZeroNewMax(max_row_count);
	Fill(zero, 0.0, zero_max());

	ResNewMax(max_row_count);
	Fill(res, 0.0, res_max());

	Delta1NewMax(max_column_count);
	Fill(delta1, 0.0, delta1_max());

	// Copy equations to optimize into ineq_array
  count_rows = 0;
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];

    if (md->iterations < md->aqueous_only)
      continue;

#ifdef DEBUG_MOHID
		if (d.debug_status == 1)
			fprintf(d.ineq_f, "Ineq-9) %d %20.20e %20.20e\n", 
																x->type,
																x->moles,
																x->f);
#endif

			
		if (x->type == PP)
    {   
      if (!x->pure_phase->phase->in)
				continue;
      if (x->pure_phase->force_equality)
				continue;
      
			if (x->f > 0e-8 && x->moles <= 0 && x->pure_phase->add_formula.IsEmpty())
				continue;
      else if (x->f < 0e-8 && x->dissolve_only && (x->moles - x->pure_phase->initial_moles >= 0))
				continue;
      else
      {
				memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
				back_eq[count_rows] = i;
	
				if (x->pure_phase->add_formula.IsEmpty() && !x->dissolve_only)
				{
				  res[count_rows] = 1.0;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-10) %20.20e %d\n", res[count_rows], count_rows);
#endif
				}

				if (md->pp_scale != 1)
					for (j = 0; j < count_unknowns + 1; j++)
					{
						ineq_array[count_rows * max_column_count + j] *= md->pp_scale;

#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.ineq_f, "Ineq-11) %20.20e %d\n", ineq_array[count_rows * max_column_count + j], count_rows * max_column_count + j);
#endif					
					}

				if (in_kode != 1)
				{
					res[count_rows] = 0.0;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-12) %20.20e\n", res[count_rows], count_rows);
#endif
				}

				count_rows++;
      }
    }
    else if (x->type == ALK || x->type == SOLUTION_PHASE_BOUNDARY)
    {
      memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
      back_eq[count_rows] = i;
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-13) %d %d\n", back_eq[count_rows], count_rows);
#endif
      
			count_rows++;
    }
    else if (x->type == GAS_MOLES && md->gas_in)
    {
      memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
      back_eq[count_rows] = i;
      res[count_rows] = 1.0;

			if (in_kode != 1)
				res[count_rows] = 0.0;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-14) %d %20.20e %d\n", back_eq[count_rows], res[count_rows], count_rows);
#endif
				
			count_rows++;
    }
    else if (x->type == S_S_MOLES && x->s_s_in)
    {
      memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
      back_eq[count_rows] = i;
      res[count_rows] = 1.0;
      
			if (in_kode != 1)
				res[count_rows] = 0.0;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-15) %d %20.20e %d\n", back_eq[count_rows], res[count_rows], count_rows);
#endif
				
      count_rows++;
    }
  }
  
	count_optimize = count_rows;

	// Copy equality equations into ineq_array
	Unknown *x_b;
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];
		
		if (i > 0) 
			x_b = (*unknown_list)[i - 1];

		if (x->type != SOLUTION_PHASE_BOUNDARY && x->type != ALK && x->type != GAS_MOLES && x->type != S_S_MOLES)
    {
      if (x->type == PP && !x->pure_phase->force_equality)
				continue;

      if (md->mass_water_switch && x == md->mass_oxygen_unknown)
				continue;

      if (x->type == EXCH && x->moles <= MIN_RELATED_SURFACE)
				continue;

      if (x->type == SURFACE && x->phase_unknown == NULL && x->moles <= MIN_RELATED_SURFACE)
				continue;

      if ((x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2) && x->surface_charge->grams <= MIN_RELATED_SURFACE)
				continue;

			if (x->type == SURFACE && x->phase_unknown != NULL && x->phase_unknown->moles <= MIN_RELATED_SURFACE && x->phase_unknown->pure_phase->add_formula.IsEmpty())
				continue;

      if ((x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2) && x_b->phase_unknown != NULL && x->surface_charge->grams <= MIN_RELATED_SURFACE && x_b->phase_unknown->pure_phase->add_formula.IsEmpty())
				continue;


#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-16) %d %d %d %20.20e %20.20e %s\n", i, x->number, x->type, x->moles, x->f, x->name.CharPtr());
#endif		

      memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
      back_eq[count_rows] = i;
      
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-17) %d %d\n", count_rows, back_eq[count_rows]);
#endif	
			
			if (md->mass_water_switch && x == md->mass_hydrogen_unknown)
      {
				k = md->mass_oxygen_unknown->number;
				
				for (j = 0; j < count_unknowns; j++)
				{
					ineq_array[count_rows * max_column_count + j] -= 2 * arr[k * (count_unknowns + 1) + j];
					
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-18) %d %d %20.20e %20.20e\n", count_rows * max_column_count + j, k * (count_unknowns + 1) + j, ineq_array[count_rows * max_column_count + j], arr[k * (count_unknowns + 1) + j]);
#endif					
				}
      }

			count_rows++;
    }
    else if (x->type == PITZER_GAMMA)
    {
      memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
      back_eq[count_rows] = i;
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-19) %d %d\n", count_rows, back_eq[count_rows]);
#endif				
    }
  }

  count_equal = count_rows - count_optimize;

	// Copy inequality constraints into ineq
  if (md->pure_phase_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == PP)
      {
				if (x->pure_phase->phase->in == FALSE)
					continue;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-20) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif	

				if (x->moles <= 0.0 && x->f > 0e-8 && x->pure_phase->add_formula.IsEmpty())
					continue;
				else if (x->moles <= 0.0)
				{
					delta1[i] = -1.0;
					
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-21) %d %20.20e\n", i, delta1[i]);
#endif			
				}
				else if (x->f < 0e-8 && x->dissolve_only && (x->moles - x->pure_phase->initial_moles >= 0))
					continue;
				else
				{
					memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &zero[0], (size_t) (count_unknowns + 1) * sizeof (LDBLE));
					
					ineq_array[count_rows * max_column_count + i] = 1.0;
					ineq_array[count_rows * max_column_count + count_unknowns] = x->moles;
					back_eq[count_rows] = i;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-22) %d %20.20e %d %20.20e %d %d\n", 
																										count_rows * max_column_count + i, 
																										ineq_array[count_rows * max_column_count + i],
																										count_rows * max_column_count + count_unknowns,
																										ineq_array[count_rows * max_column_count + count_unknowns],
																										back_eq[count_rows],
																										count_rows);
#endif							
					count_rows++;
				}
	
				if (x->dissolve_only)
				{
					memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &zero[0], (size_t) (count_unknowns + 1) * sizeof (LDBLE));

					ineq_array[count_rows * max_column_count + i] = -1.0;
					ineq_array[count_rows * max_column_count + count_unknowns] = x->pure_phase->initial_moles - x->moles;
					back_eq[count_rows] = i;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-23) %d %20.20e %d %20.20e %d %d\n", 
																										count_rows * max_column_count + i, 
																										ineq_array[count_rows * max_column_count + i],
																										count_rows * max_column_count + count_unknowns,
																										ineq_array[count_rows * max_column_count + count_unknowns],
																										back_eq[count_rows],
																										count_rows);
#endif					
					count_rows++;								
				}
      }
    }
  }

	/*
	 *   Hydrogen mass balance is good
	 */
	/*
	 *   No moles and undersaturated, mass transfer must be zero
	 */
	if (md->pure_phase_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == PP)
      {
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-24) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif			

				if ((x->moles <= 0.0 && x->f > 0e-8 && x->pure_phase->add_formula.IsEmpty()) || !x->pure_phase->phase->in)
					for (j = 0; j < count_rows; j++)
					{
						ineq_array[j * max_column_count + i] = 0.0;
						
#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.ineq_f, "Ineq-25) %d %20.20e\n", j * max_column_count + i, ineq_array[j * max_column_count + i]);
#endif									
					}

				if (x->dissolve_only && x->f < 0e-8 && (x->moles - x->pure_phase->initial_moles >= 0))				
					for (j = 0; j < count_rows; j++)
					{
						ineq_array[j * max_column_count + i] = 0.0;
						
#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.ineq_f, "Ineq-26) %d %20.20e\n", j * max_column_count + i, ineq_array[j * max_column_count + i]);
#endif									
					}
      }
    }
  }

	// No moles of exchanger
	if (md->use.exc_p != NULL && (md->use.exc_p->related_phases || md->use.exc_p->related_rate))
    for (i = 0; i < count_unknowns; i++)
		{
			x = (*unknown_list)[i];

      if (x->type == EXCH && x->moles <= 0)
			{
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-27) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif			

				for (j = 0; j < count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
					
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-28) %d %20.20e\n", j * max_column_count + i, ineq_array[j * max_column_count + i]);
#endif									
					
				}
			}
		}

	// No moles of surface
  if (md->use.sur_p != NULL && (md->use.sur_p->related_phases || md->use.sur_p->related_rate))
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];
			
			if (i < 0)
				x_b = (*unknown_list)[i - 1];
			else
				x_b = NULL;

      if ((x->type == SURFACE && x->phase_unknown != NULL &&  x->phase_unknown->moles <= MIN_RELATED_SURFACE && x->phase_unknown->pure_phase->add_formula.IsEmpty()) || ((x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2) && x_b->phase_unknown != NULL && x->surface_charge->grams <= MIN_RELATED_SURFACE && x_b->phase_unknown->pure_phase->add_formula.IsEmpty()))
			{
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-29) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif		
				for (j = 0; j < count_rows; j++)
				{				
				  ineq_array[j * max_column_count + i] = 0.0;
					
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-30) %d %20.20e\n", j * max_column_count + i, ineq_array[j * max_column_count + i]);
#endif						
				}
			}
    }

	// No moles of surface
  if (md->use.sur_p != NULL)
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];
			
			if (i > 0)
				x_b = (*unknown_list)[i - 1];
			else
				x_b = NULL;		
				
      if ((x->type == SURFACE && x->phase_unknown == NULL && x->moles <= MIN_RELATED_SURFACE) || ((x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2) && x_b->phase_unknown == NULL && x->surface_charge->grams <= MIN_RELATED_SURFACE))
			{
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-31) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif	

				for (j = 0; j < count_rows; j++)
				{
				  ineq_array[j * max_column_count + i] = 0.0;   
		
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-32) %d %20.20e\n", j * max_column_count + i, ineq_array[j * max_column_count + i]);
#endif				
				}
			}
		}

	// Moles of gas must be >= zero
  if (md->gas_in)
  {
    i = md->gas_unknown->number;

    memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &zero[0], (size_t) (count_unknowns + 1) * sizeof (LDBLE));

		ineq_array[count_rows * max_column_count + i] = -1.0;
    ineq_array[count_rows * max_column_count + count_unknowns] = x->moles;
    back_eq[count_rows] = i;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.ineq_f, "Ineq-33) %d %20.20e %d %20.20e %d %d\n", 
																							count_rows * max_column_count + i, 
																							ineq_array[count_rows * max_column_count + i],
																							count_rows * max_column_count + count_unknowns,
																							ineq_array[count_rows * max_column_count + count_unknowns],
																							back_eq[count_rows],
																							count_rows);
#endif		
    count_rows++;
  }
  else if (md->use.gas_p != NULL && !md->gas_in)
  {
    i = md->gas_unknown->number;
    
		for (j = 0; j < count_rows; j++)
		{
      ineq_array[j * max_column_count + i] = 0.0;
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-34) %d %20.20e\n", 
																								j * max_column_count + i, 
																								ineq_array[j * max_column_count + i]);
#endif	
		}
  }

	// Phase must be "in" and moles of solid solution must be >= zero
  if (md->s_s_unknown != NULL)
    for (i = md->s_s_unknown->number; i < count_unknowns; i++)
    {
			//ToDo: this code was copied without the definition of x, so, x would be anything...
			
			x = (*unknown_list)[i];
			
      if (x->type != S_S_MOLES)
				break;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-35) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif	
				
      if (x->p->in && x->s_s_in)
      {
				memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &zero[0], (size_t) (count_unknowns + 1) * sizeof (LDBLE));

				ineq_array[count_rows * max_column_count + i] = 1.0;
				ineq_array[count_rows * max_column_count + count_unknowns] = (LDBLE)0.99 * x->moles - (LDBLE)MIN_TOTAL_SS;
				back_eq[count_rows] = i;
				
				
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-36) %d %20.20e %d %20.20e %d %d\n", 
																			count_rows * max_column_count + i, 
																			ineq_array[count_rows * max_column_count + i],
																			count_rows * max_column_count + count_unknowns,
																			ineq_array[count_rows * max_column_count + count_unknowns],
																			back_eq[count_rows],
																			count_rows);
#endif					
				count_rows++;
      }
      else
				for (j = 0; j < count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
					
#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.ineq_f, "Ineq-37) %d %20.20e\n", 
																										j * max_column_count + i, 
																										ineq_array[j * max_column_count + i]);
#endif						
				}
    }

	// Add inequality if total moles of element is less than zero
  if (md->negative_concentrations)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == MB && x->moles < 0.0)
      {

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-38) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif	
			
				memcpy ((void *) &ineq_array[count_rows * max_column_count], (void *) &arr[i * (count_unknowns + 1)], (size_t) (count_unknowns + 1) * sizeof (LDBLE));				
				back_eq[count_rows] = i;

				for (j = 0; j < count_unknowns; j++)
					if ((*unknown_list)[j]->type < PP)
					{
						ineq_array[count_rows * max_column_count + j] = 0.0;
						
#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.ineq_f, "Ineq-39) %d %20.20e\n", 
																					count_rows * max_column_count + j, 
																					ineq_array[count_rows * max_column_count + j]);
#endif									
					}

				count_rows++;
      }
    }
  }

	// Zero column for mass of water
  if (md->mass_oxygen_unknown != NULL && md->mass_water_switch)
  {
    k = md->mass_oxygen_unknown->number;

    for (j = 0; j < count_rows + 1; j++)
		{
      ineq_array[j * max_column_count + k] = 0;
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-40) %d %20.20e\n", 
																		j * max_column_count + k, 
																		ineq_array[j * max_column_count + k]);
#endif									
			
		}
  }

	// Scale column for pure phases
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];

    if ((x->type == PP || x->type == S_S_MOLES)	&& md->pp_column_scale != 1.0)
    {
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-41) %d %d %d %20.20e %20.20e\n", i, x->number, x->type, x->moles, x->f);
#endif

      for (j = 0; j < count_rows; j++)
			{
				ineq_array[j * max_column_count + i] *= md->pp_column_scale;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-42) %d %20.20e\n", 
																			j * max_column_count + i, 
																			ineq_array[j * max_column_count + i]);
#endif			
			}

			normal[i] = md->pp_column_scale;
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.ineq_f, "Ineq-43) %d %20.20e\n", 
																								i, 
																								normal[i]);
#endif			
    }
  }

	// Calculate dimensions
  k = count_optimize;		
  l = count_equal;		
  m = count_rows - l - k;	
	
  if (m < 0)
    m = 0;

  n = count_unknowns;		
  klmd = max_row_count - 2;
  nklmd = n + klmd;
  n2d = n + 2;

  kode = 1;

  if (in_kode == 2)
    kode = 1;
				
  iter = 200;

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.ineq_f, "Ineq-44) %d %d %d %d %d %d\n", 
																k, l, m, klmd, nklmd, n2d);
#endif
	
	// Allocate space for arrays
	CuNewMax(2 * nklmd);
	IuNewMax(2 * nklmd);
	IsNewMax(klmd);

  CL1(k, l, m, n, nklmd, n2d, &ineq_array[0], &kode, md->ineq_tol, &iter, &delta1[0], &res[0], &error, &cu[0], &iu[0], &is[0], false);

  if (kode == 1)
    return_code = ERROR;
  else if (kode == 2)
    return_code = 2;
  else
    return_code = OK;

	// Copy delta1 into delta and scale
  memcpy ((void *) &delta[0], (void *) &delta1[0], (size_t) count_unknowns * sizeof (LDBLE));
  
	for (i = 0; i < count_unknowns; i++)
    delta[i] *= normal[i];

#ifdef DEBUG_MOHID
		if (d.debug_status)
			for (int debug_i = 0; debug_i < count_unknowns; debug_i++)
				fprintf(d.ineq_f, "Ineq-45) %20.20e\n", delta[debug_i]);
#endif
		
	// Rescale columns of array
  for (i = 0; i < count_unknowns; i++)
    if (normal[i] != 1.0)
      for (j = 0; j < count_unknowns; j++)
			{
				arr[j * (count_unknowns + 1) + i] /= normal[i];
				
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.ineq_f, "Ineq-46) %d %20.20e %20.20e %d\n", 
																			j * (count_unknowns + 1) + i,
																			arr[j * (count_unknowns + 1) + i],
																			normal[i],
																			i);
#endif				
			}

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.ineq_f, "Ineq-FIM)====================\n");
#endif

  return (return_code);
}

bool ModelEngine::Reset()
{
	//
	// Checks deltas (changes to unknowns) to make sure they are reasonable
	// Scales deltas if necessary 
	// Updates unknowns with deltas
	//

  int i, j;
  bool converge;
  LDBLE up, down;
  LDBLE d_;
  LDBLE factor, f0;
  LDBLE sum_deltas;
  LDBLE step_up;
  LDBLE mu_calc;
  LDBLE old_moles;
	Master *m;

	/*
	if (d.si_count == 132)
		getch();
		*/

	// Calculate interphase mass transfers
	step_up = log(md->step_size_now);
  factor = 1.0;

	Unknown *x;
	STCoef *stc_p;

  if ((md->pure_phase_unknown != NULL || md->s_s_unknown != NULL) && !md->calculating_deriv)
  {
		// Don't take out more mineral than is present
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == PP || x->type == S_S_MOLES)
      {
				if (delta[i] < -1e8) 
					delta[i] = -10.0;
				else if (delta[i] > 1e8)
				  delta[i] = 10.0;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reset_f, "reset-1) %s %d %d %20.20e %20.20e\n", x->name.CharPtr(), x->number, x->type, x->moles, delta[i]); 
#endif

				if (x->dissolve_only)
				{
					if (delta[i] < 0.0 && (-delta[i] > (x->pure_phase->initial_moles - x->moles)))
					{
#ifdef DEBUG_MOHID
						if (d.debug_status)
							fprintf(d.reset_f, "reset-2) %s %20.20e\n", x->pure_phase->name.CharPtr(), x->pure_phase->initial_moles); 
#endif

						if ((x->pure_phase->initial_moles - x->moles) != 0.0)
						{
							f0 = fabs (delta[i] / (x->pure_phase->initial_moles - x->moles));
	      
							if (f0 > factor)
								factor = f0;

#ifdef DEBUG_MOHID
							if (d.debug_status)
								fprintf(d.reset_f, "reset-3) %20.20e %20.20e\n", f0, factor); 
#endif
						}
						else
						{
							delta[i] = 0;
						
#ifdef DEBUG_MOHID
							if (d.debug_status)
								fprintf(d.reset_f, "reset-4) %20.20e\n", delta[i]); 
#endif
						}
					}
				}

				if (x->moles > 0.0 && delta[i] > x->moles)
				{
					f0 = delta[i] / x->moles;

					if (f0 > factor)
						factor = f0;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-5) %20.20e %20.20e\n", f0, factor); 
#endif
				}
				else if (delta[i] > 0.0 && x->moles <= 0.0)
				{
					delta[i] = 0.0;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-6) %20.20e\n", delta[i]); 
#endif
				}
				else if (delta[i] < (LDBLE)-100.0)
				{
					f0 = -delta[i] / (LDBLE)100.0;

					if (f0 > factor)
						factor = f0;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-7) %20.20e %20.20e\n", f0, factor); 
#endif
				}
			}
    }
  }

	// Calculate change in element concentrations due to pure phases and gases
	for (i = 0; i < count_unknowns; i++)
    if (_isnan(delta[i]))
		{
      delta[i] = 0;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-8) %d %20.20e\n", i, delta[i]); 
#endif		
		}

	if (md->pure_phase_unknown != NULL || md->gas_unknown != NULL || md->s_s_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
      (*unknown_list)[i]->delta = 0.0;

    for (i = 0; i < md->sum_delta->Count(); i++)
		{
			stc_p = (*md->sum_delta)[i];
      *stc_p->target += *stc_p->source * stc_p->coef;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-9) %d %20.20e %20.20e %20.20e\n", i, *stc_p->target, *stc_p->source, stc_p->coef); 
#endif	
		}

		// Apply factor from minerals to deltas
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      x->delta /= factor;

      if (x->type == PP || x->type == S_S_MOLES) 
				delta[i] /= factor;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reset_f, "reset-10) %s %d %d %20.20e %20.20e\n", 
																				x->name.CharPtr(), 
																				x->number, x->type, 
																				x->delta, delta[i]); 
#endif
    }
  }

	// Calc factor for mass balance equations for aqueous unknowns
  factor = 1.0;
  sum_deltas = 0.0;
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];

		// fixes underflow problem on Windows
    if (delta[i] > 0)
      sum_deltas += delta[i];
    else
      sum_deltas -= delta[i];

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.reset_f, "reset-11) %d %d %20.20e %20.20e %20.20e %20.20e\n", i, x->type, x->moles, sum_deltas, step_up, md->mu_x); 
#endif	

		if (!md->calculating_deriv)
    {
      up = step_up;
      down = up;

      if (x->type <= SOLUTION_PHASE_BOUNDARY)
      {
				up = step_up;
				down = (LDBLE)1.3 * up;
      }
      else if (x->type == MU)
      {
				up = 100 * md->mu_x;
				down = md->mu_x;
      }
      else if (x->type == AH2O)
				down = up;
      else if (x->type == MH)
      {
				up = log (md->pe_step_size_now);
				down = (LDBLE)1.3 * up;
      }
      else if (x->type == MH2O)
      {
				up = log ((LDBLE)1.3);
				down = log ((LDBLE)1.2);
      }
      else if (x->type == PP)
				continue;
      else if (x->type == GAS_MOLES)
      {
				up = (LDBLE)1000. * x->moles;

				if (up <= 0.0)
					up = (LDBLE)1e-1;

				if (up >= 1.0)
					up = 1.0;

				down = x->moles;
      }
      else if (x->type == S_S_MOLES)
				continue;
      else if (x->type == EXCH)
      {
				up = step_up;
				down = (LDBLE)1.3 * up;
      }
      else if (x->type == SURFACE)
      {
				up = step_up;
				down = (LDBLE)1.3 * up;
      }
      else if (x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2)
      {
				up = step_up;
				down = (LDBLE)1.3 * up;
      }

      if (delta[i] > 0.0)
      {
				f0 = delta[i] / up;

				if (f0 > factor)
					factor = f0;
      }
      else
      {
				f0 = delta[i] / (-down);

				if (f0 > factor)
					factor = f0;
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-12) %20.20e %20.20e\n", f0, factor); 
#endif	
    }
  }

  factor = (LDBLE)1.0 / factor;

  for (i = 0; i < count_unknowns; i++)
	{
    if ((*unknown_list)[i]->type != PP && (*unknown_list)[i]->type != S_S_MOLES)
      delta[i] *= factor;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-13) %20.20e %20.20e\n", delta[i], factor); 
#endif	
	}

	SurfaceDiffLayer *sdl_p;
	// Solution mass balances: MB, ALK, CB, SOLUTION_PHASE_BOUNDARY
  for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];		

		if (x->master->Count() > 0)
			m = (*x->master)[0];
		else
			m = NULL;


#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.reset_f, "reset-14) %d %d %20.20e %20.20e\n", 
																		i, x->type, x->moles, delta[i]); 
#endif

    if (x->type == MB || x->type == ALK || x->type == EXCH || x->type == SURFACE)
    {			
      d_ = delta[i] / md->LOG_10;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.reset_f, "reset-15) %20.20e\n", 
																		d_); 
#endif
			if (x->type == SURFACE)
      {
				old_moles = x->moles;

				if (x->phase_unknown != NULL)
				{
					x->moles = x->surface_comp->phase_proportion * (x->phase_unknown->moles - delta[x->phase_unknown->number]);
					
					if (x->phase_unknown->moles - delta[x->phase_unknown->number] <= MIN_RELATED_SURFACE)
					{						
						x->moles = 0.0;

						if (fabs (x->f) > MIN_RELATED_SURFACE)
							m->s->la -= 5.;
					}

					if (old_moles <= 0 && x->moles > 0)
						m->s->la = log10 (x->moles) - (LDBLE)5.;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-16) %20.20e %s %20.20e %d %20.20e %20.20e %20.20e %20.20e\n", 
																					x->surface_comp->phase_proportion, 
																					x->phase_unknown->name.CharPtr(),
																					x->phase_unknown->moles,
																					x->phase_unknown->number,
																					delta[x->phase_unknown->number],
																					x->moles,
																					x->f,
																					m->s->la); 
#endif
				}
      }

      if (x->type == EXCH && x->moles <= MIN_RELATED_SURFACE)
      {
				x->moles = 0.0;

				if (fabs (x->f) > MIN_RELATED_SURFACE)
					m->s->la -= 5.;


#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reset_f, "reset-17) %20.15 %20.20e\n", 
																				x->moles,
																				m->s->la); 
#endif
      }

			m->s->la += d_;

      if (m->s->la < (LDBLE) (DBL_MIN_10_EXP + 10))
				m->s->la = (LDBLE) (DBL_MIN_10_EXP + 10);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-18) %20.20e\n", 
																			m->s->la); 
#endif
    }
    else if (x->type == SURFACE_CB || x->type == SURFACE_CB1 || x->type == SURFACE_CB2)
    {
      d_ = delta[i] / md->LOG_10;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-19) %20.20e\n", 
																			d_); 
#endif

      if (x->phase_unknown != NULL)
      {
				x->surface_charge->grams = (x->phase_unknown->moles - delta[x->phase_unknown->number]);

				if (x->surface_charge->grams <= MIN_RELATED_SURFACE)
					x->surface_charge->grams = 0.0;
      }

			if (x->surface_charge->grams <= MIN_RELATED_SURFACE)
				x->surface_charge->grams = 0.0;

			x->related_moles = x->surface_charge->grams;

			m->s->la += d_;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-20) %20.20e %s %20.20e %d %20.20e %20.20e %20.20e\n", 
																			x->surface_charge->grams, 
																			x->phase_unknown->name.CharPtr(),
																			x->phase_unknown->moles,
																			x->phase_unknown->number,
																			delta[x->phase_unknown->number],
																			x->related_moles,
																			m->s->la); 
#endif

			// recalculate g's for component
      if (md->dl_type_x != NO_DL && (md->use.sur_p->type == DDL || (md->use.sur_p->type == CD_MUSIC && x->type == SURFACE_CB2)))
      {
				for (j = 0; j < x->surface_charge->g->Count(); j++)
				{
					sdl_p = (*x->surface_charge->g)[j];

					if (md->use.sur_p->dl_type != DONNAN_DL)
						sdl_p->g += sdl_p->dg * delta[i];

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-21) %20.20e %20.20e\n", 
																					sdl_p->g, 
																					delta[i]); 
#endif
				}

				//ToDo: Put in "calcalldonnan file a "mark" saying that was called from here
				if (md->use.sur_p->dl_type == DONNAN_DL)
					CalcAllDonnan ();
      }
    }
    else if (x->type == SOLUTION_PHASE_BOUNDARY)
    {
      d_ = delta[i] / md->LOG_10;
      m->s->la += d_;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-22) %20.20e %20.20e\n", 
																					d_, 
																					m->s->la); 
#endif
    }
    else if (x->type == CB)
    {
      d_ = delta[i] / md->LOG_10;
      m->s->la += d_;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reset_f, "reset-23) %20.20e %20.20e\n", 
																					d_, 
																					m->s->la); 
#endif
    }
    else if (x->type == MU)
    {
      mu_calc = (LDBLE)0.5 * md->mu_unknown->f / md->mass_water_aq_x;
      d_ = md->mu_x + delta[i];

			if (d_ < 1e-7)
      {
				delta[i] = sqrt(mu_calc * md->mu_x) - md->mu_x;
				md->mu_x = sqrt(mu_calc * md->mu_x);
				//d.PrintSpeciesInfoToFile("Reset-1", md->species_info_list, md, gd);
      }
      else
			{
				md->mu_x += delta[i];
				//d.PrintSpeciesInfoToFile("Reset-2", md->species_info_list, md, gd);
			}

			if (md->mu_x <= 1e-8)
			{
				md->mu_x = (LDBLE)1e-8;
				//d.PrintSpeciesInfoToFile("Reset-3", md->species_info_list, md, gd);
			}

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-24) %20.20e %20.20e %20.20e %20.20e\n", 
																			d_, 
																			mu_calc,
																			delta[i],
																			md->mu_x); 
#endif
    }
    else if (x->type == AH2O)
    {
      d_ = delta[i] / md->LOG_10;
      gd->s_h2o->la += d_;

			if (gd->s_h2o->la < -1.0)
			{
				d_ = (LDBLE)-1.0 - gd->s_h2o->la;
				delta[i] = d_ * md->LOG_10;
				gd->s_h2o->la = -1.0;
			}

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-25) %20.20e %20.20e %20.20e\n", 
																			d_, 
																			gd->s_h2o->la,
																			delta[i]); 
#endif
    }
    else if (x->type == MH)
    {
      d_ = delta[i] / md->LOG_10;
      gd->s_eminus->la += d_;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-26) %20.20e %20.20e %20.20e\n", 
																			d_, 
																			gd->s_eminus->la,
																			delta[i]); 
#endif
    }
    else if (x->type == MH2O)
    {
      if (md->mass_water_switch)
				continue;

			d_ = exp(delta[i]);
      md->mass_water_aq_x *= d_;
      md->mass_water_bulk_x = md->mass_water_aq_x + md->mass_water_surfaces_x;
      m->s->moles = md->mass_water_aq_x / md->gfw_water;

			if (md->use.sur_p != NULL && md->use.sur_p->debye_lengths > 0)
				m->s->moles = md->mass_water_bulk_x / md->gfw_water;

      if (md->mass_water_aq_x < 1e-10)
				return false;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-27) %20.20e %20.20e %20.20e %20.20e\n", 
																			d_, 
																			md->mass_water_aq_x,
																			md->mass_water_bulk_x,
																			m->s->moles); 
#endif
    }
    else if (x->type == PP)
    {
      if (Equal(x->moles, delta[i], md->ineq_tol))
				x->moles = 0.0;
      else
				x->moles -= delta[i];

      if (x->dissolve_only && Equal(x->moles, x->pure_phase->initial_moles, md->ineq_tol))
				x->moles = x->pure_phase->initial_moles;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-28) %20.20e %20.20e %20.20e %20.20e\n", 
																			x->moles, 
																			x->pure_phase->initial_moles,
																			md->ineq_tol,
																			delta[i]); 
#endif
    }
    else if (x->type == GAS_MOLES)
    {
      x->moles += delta[i];
      
			if (x->moles < MIN_TOTAL)
				x->moles = (LDBLE)MIN_TOTAL;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-29) %20.20e %20.20e\n", 
																			x->moles, 
																			delta[i]); 
#endif
    }
    else if (x->type == S_S_MOLES)
    {
      x->moles -= delta[i];

			if (x->moles < MIN_TOTAL_SS && !md->calculating_deriv)
				x->moles = (LDBLE)MIN_TOTAL_SS;

      x->s_s_comp->moles = x->moles;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-30) %20.20e %20.20e\n", 
																			x->moles, 
																			delta[i]); 
#endif
    }
    else if (x->type == PITZER_GAMMA)
    {
      if (!md->full_pitzer)
				continue;

      d_ = delta[i];
      x->s->lg += d_;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-31) %20.20e %20.20e\n", 
																			x->s->lg, 
																			delta[i]); 
#endif
    }
  }

	// Reset total molalities in mass balance equations
	if (md->pure_phase_unknown != NULL || md->gas_unknown != NULL || md->s_s_unknown != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == MB || x->type == MH || x->type == MH2O || x->type == CB || x->type == EXCH || x->type == SURFACE)
      {
				if (x->type == SURFACE)
					x->delta = 0.0;

				x->moles += x->delta;
      }

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reset_f, "reset-32) %d %20.20e %20.20e\n", 
																			i,
																			x->delta, 
																			x->moles); 
#endif
    }
  }

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.reset_f, "FIM===)\n"); 
#endif
  converge = false;
  return (converge);
}

bool ModelEngine::Gammas(LDBLE mu)
{
	//
	// Calculates gammas and [moles * d(ln gamma)/d mu] for all aqueous
	// species.
	//
	//ToDo: Retirar a parte do debug

	Species *s;
  int i, j;
  LDBLE a_llnl, b_llnl, bdot_llnl, log_g_co2, dln_g_co2, c2_llnl;
  LDBLE s1, s2, s3;
  LDBLE c1, c2, a, b;
  LDBLE muhalf, equiv;

	a_llnl = b_llnl = bdot_llnl = log_g_co2 = dln_g_co2 = c2_llnl = 0;

	// compute temperature dependence of a and b for debye-huckel
  s1 = (LDBLE)374.11 - md->tc_x;
  s2 = pow ((LDBLE)s1, (LDBLE)(1.0 / 3.0));
  s3 = (LDBLE)1.0 + (LDBLE)0.1342489 * s2 - (LDBLE)3.946263e-03 * s1;
  s3 = s3 / ((LDBLE)3.1975 - (LDBLE)0.3151548 * s2 - (LDBLE)1.203374e-03 * s1 + (LDBLE)7.48908e-13 * (s1 * s1 * s1 * s1));
  s3 = sqrt(s3);

  if (md->tk_x >= 373.15)
    c1 = (LDBLE)5321.0 / md->tk_x + (LDBLE)233.76 - md->tk_x * (md->tk_x * ((LDBLE)8.292e-07 * md->tk_x - (LDBLE)1.417e-03) + (LDBLE)0.9297);
  else
    c1 = (LDBLE)2727.586 + (LDBLE)0.6224107 * md->tk_x - (LDBLE)466.9151 * log (md->tk_x) - (LDBLE)52000.87 / md->tk_x;

	c1 = sqrt (c1 * md->tk_x);

	a = (LDBLE)1824827.7 * s3 / (c1 * c1 * c1);
  b = (LDBLE)50.2905 * s3 / c1;

	// constants for equations
  muhalf = sqrt(mu);
  c1 = (-a) * md->LOG_10 * ((LDBLE)1.0 / ((LDBLE)2 * muhalf * (muhalf + (LDBLE)1.0) * (muhalf + (LDBLE)1.0)) - (LDBLE)0.3);
  c2 = -a / (2 * muhalf);

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.gammas_f, "1) %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e \n", 
														s1, s2, s3, c1, c2, a, b, mu, muhalf);
#endif 

	// Calculate activity coefficients
	Species *sx;
  for (i = 0; i < md->s_x->Count(); i++)
  {
		sx = (*md->s_x)[i];

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.gammas_f, "2) number:%i gflag:%d %20.20e %20.20e %20.20e %20.20e %20.20e %20.20e\n", 
														sx->number, sx->gflag, 
														sx->dhb, 
														sx->moles, 
														sx->z, 
														sx->dha, 
														gd->s_h2o->la, 
														md->gfw_water);
#endif
    switch (sx->gflag)
    {
    case 0: // uncharged
      sx->lg = sx->dhb * mu;
      sx->dg = sx->dhb * md->LOG_10 * sx->moles;
      break;
    case 1: // Davies			
      sx->lg = -sx->z * sx->z * a *	(muhalf / ((LDBLE)1.0 + muhalf) - (LDBLE)0.3 * mu);
      sx->dg = c1 * sx->z * sx->z * sx->moles;
      break;
    case 2: // Extended D-H, WATEQ D-H
      sx->lg = -a * muhalf * sx->z * sx->z / ((LDBLE)1.0 + sx->dha * b * muhalf) + sx->dhb * mu;
      sx->dg = (c2 * sx->z * sx->z / (((LDBLE)1.0 + sx->dha * b * muhalf) * ((LDBLE)1.0 + sx->dha * b * muhalf)) + sx->dhb) * md->LOG_10 * sx->moles;
      break;
    case 3: // Always 1.0
      sx->lg = 0.0;
      sx->dg = 0.0;
      break;
    case 4: // Exchange
			// Find CEC
			// z contains valence of cation for exchange species, alk contains cec
      for (j = 1; j < sx->rxn_x->token_list->Count(); j++)
			{
				s = sx->rxn_x->token_list->Element(j)->s;

				if (s->type == EX)
				{
					sx->alk = s->primary->u->moles;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.gammas_f, "3) j:%d s_number:%d type:%d primary_number:%d u_number:%d %20.20e %20.20e\n", 
																	i, 
																	(*sx->rxn_x->token_list)[j]->s->number, 
																	(*sx->rxn_x->token_list)[j]->s->type, 
																	(*sx->rxn_x->token_list)[j]->s->primary->number,
																	(*sx->rxn_x->token_list)[j]->s->primary->u->number,
																	(*sx->rxn_x->token_list)[j]->s->primary->u->moles,
																	sx->alk);
#endif

					break;
				}
			}
			
#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.gammas_f, "4) exch_flag:%d %20.20e\n", sx->exch_gflag, sx->equiv);
#endif

			if (sx->exch_gflag == 1 && sx->alk > 0) // Davies
      {
				sx->lg = -sx->equiv * sx->equiv * a * (muhalf / ((LDBLE)1.0 + muhalf) - (LDBLE)0.3 * mu) + log10 (fabs (sx->equiv) / sx->alk);
				sx->dg = c1 * sx->equiv * sx->equiv * sx->moles;
      }
      else if (sx->exch_gflag == 2 && sx->alk > 0) // Extended D-H, WATEQ D-H
      {
				sx->lg = -a * muhalf * sx->equiv * sx->equiv / ((LDBLE)1.0 + sx->dha * b * muhalf) + sx->dhb * mu + log10 (fabs (sx->equiv) / sx->alk);
				sx->dg = (c2 * sx->equiv * sx->equiv / (((LDBLE)1.0 + sx->dha * b * muhalf) * ((LDBLE)1.0 + sx->dha * b * muhalf)) + sx->dhb) * md->LOG_10 * sx->moles;
      }
      else if (sx->exch_gflag == 7 && sx->alk > 0)
				return false;
      else
      {
				// Master species is a dummy variable with meaningless activity and mass
				if (sx->primary != NULL)
				{
					sx->lg = 0.0;
					sx->dg = 0.0;
				}
				else
				{
					if (sx->alk <= 0)
						sx->lg = 0.0;
					else
						sx->lg = log10 (fabs (sx->equiv) / sx->alk);

					sx->dg = 0.0;
				}
      }
      break;
		case 5: // Always 1.0
      sx->lg = 0.0;
      sx->dg = 0.0;
      break;
    case 6: // Surface
			// Find moles of sites.
			// s_x[i]->equiv is stoichiometric coefficient of sites in species

      for (j = 1; j < sx->rxn_x->token_list->Count(); j++)
			{
				s = sx->rxn_x->token_list->Element(j)->s;

				if (s->type == SURF)
				{
					sx->alk = s->primary->u->moles;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.gammas_f, "5) j:%d s_number:%d type:%d primary_number:%d u_number:%d %20.20e %20.20e\n", 
																	i, 
																	(*sx->rxn_x->token_list)[j]->s->number, 
																	(*sx->rxn_x->token_list)[j]->s->type, 
																	(*sx->rxn_x->token_list)[j]->s->primary->number,
																	(*sx->rxn_x->token_list)[j]->s->primary->u->number,
																	(*sx->rxn_x->token_list)[j]->s->primary->u->moles,
																	sx->alk);
#endif

					break;
				}
			}

      if (md->use.sur_p->type == CD_MUSIC)
				equiv = 1.0; // mole fraction
      else
				equiv = sx->equiv;

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.gammas_f, "6) %20.20e\n", equiv);
#endif

			if (sx->alk > 0)
      {
				sx->lg = log10 (equiv / sx->alk);
				sx->dg = 0.0;
      }
      else
      {
				sx->lg = 0.0;
				sx->dg = 0.0;
      }
      break;
    case 7:
    case 8:			
			return false;
    case 9: // activity water

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.gammas_f, "7) %20.20e %20.20e\n", gd->s_h2o->la, md->gfw_water);
#endif

      sx->lg = log10 (exp (gd->s_h2o->la * md->LOG_10) * md->gfw_water);
      sx->dg = 0.0;
      break;
    }

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.gammas_f, "8) %20.20e %20.20e\n", sx->lg, sx->dg);
#endif
  }

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.gammas_f, "9) FIM------------\n");
#endif
  return true;
}

bool ModelEngine::Molalities(bool allow_overflow)
{
	//
	// Calculates la for master species
	// Calculates lm and moles from lk, lg, and la's of master species
	// Adjusts lm of h2 and o2.
	//

  int i, j, k, count_g;
  LDBLE total_g;
	ReactionToken *rxn_token;

	// la for master species
	Master *m;

	for (i = 0; i < gd->master_list.Count(); i++)
	{
		m = gd->master_list[i];
		if (m->rewrite)
		{
      m->s->la = m->s->lm + m->s->lg;

#ifdef DEBUG_MOHID
			if(d.debug_status)
				fprintf(d.molalities_f, "1) m:%d s:%d %20.20e %20.20e %20.20e\n", 
															    m->number,
																	m->s->number,
																	m->s->lg,
																	m->s->lm,
																	m->s->la);
#endif

		}
	}

  if (md->dl_type_x != NO_DL)
  {
    gd->s_h2o->tot_g_moles = gd->s_h2o->moles;
    gd->s_h2o->tot_dh2o_moles = 0.0;
  }

#ifdef DEBUG_MOHID
	if(d.debug_status)
		fprintf(d.molalities_f, "2) %20.20e %20.20e %20.20e\n", 
															gd->s_h2o->tot_g_moles,
															gd->s_h2o->tot_dh2o_moles);
#endif

	Species *sx;
  for (i = 0; i < md->s_x->Count(); i++)
  {
		sx = (*md->s_x)[i];

    if (sx->type > HPLUS && sx->type != EX && sx->type != SURF)
      continue;

		// lm and moles for all aqueous species
    sx->lm = sx->lk - sx->lg;

#ifdef DEBUG_MOHID
		if(d.debug_status)
			fprintf(d.molalities_f, "3) number:%d type:%d %20.20e %20.20e %20.20e\n", 															    
																sx->number,
																sx->type,
																sx->lg,
																sx->lk,
																sx->lm);
#endif
		
		for (k = 1; k < sx->rxn_x->token_list->Count(); k++)
		{
			rxn_token = sx->rxn_x->token_list->Element(k);
      sx->lm += rxn_token->s->la * rxn_token->coef;

#ifdef DEBUG_MOHID
			if(d.debug_status)
				fprintf(d.molalities_f, "4) number:%d %20.20e %20.20e\n", 															    
																	rxn_token->s->number,
																	rxn_token->s->la,
																	rxn_token->coef);
#endif
		}

		if (sx->type == EX)
      sx->moles = exp (sx->lm * md->LOG_10);
    else if (sx->type == SURF)
      sx->moles = exp (sx->lm * md->LOG_10); // formerly * mass water
    else
    {
      sx->moles = Under(sx->lm) * md->mass_water_aq_x;

      if (sx->moles / md->mass_water_aq_x > 30)
				if (md->iterations >= 0 && !allow_overflow)
 					return false;
    }

#ifdef DEBUG_MOHID
		if(d.debug_status)
			fprintf(d.molalities_f, "5) %20.20e %20.20e\n", 															    
																sx->lm,
																md->mass_water_aq_x,
																sx->moles);
#endif
  }

	// other terms for diffuse layer model
  if (md->use.sur_p != NULL && md->use.sur_p->type == CD_MUSIC && md->dl_type_x != NO_DL)
    CalcAllDonnan();

  for (i = 0; i < md->s_x->Count(); i++)
  {
		sx = (*md->s_x)[i];

    if (sx->type > HPLUS && sx->type != EX && sx->type != SURF)
      continue;

    if (md->use.sur_p != NULL && md->dl_type_x != NO_DL	&& sx->type <= HPLUS)
    {
      total_g = 0.0;
      sx->tot_dh2o_moles = 0.0;

			SpeciesDiffLayer *sdl;
			SurfaceDiffLayer *surdl_p;
			SurfaceCharge *sc_p;
			for (j = 0; j < md->use.sur_p->charge->Count(); j++)
      {
				sc_p = (*md->use.sur_p->charge)[j];
				sdl = (*sx->diff_layer)[j];

				count_g = sdl->count_g;
				surdl_p = (*sdl->charge->g)[count_g];

				sdl->g_moles = sx->moles * (surdl_p->g + sdl->charge->mass_water / md->mass_water_aq_x);

				if (sx->moles > 1e-30)
					sdl->dg_g_moles = sx->dg * sdl->g_moles / sx->moles;

				total_g += surdl_p->g + sdl->charge->mass_water / md->mass_water_aq_x;

				sdl->dx_moles = sx->moles * surdl_p->dg;

				sdl->dh2o_moles = -sx->moles * sdl->charge->mass_water / md->mass_water_aq_x;
				sx->tot_dh2o_moles += sdl->dh2o_moles;

				sdl->drelated_moles = sx->moles * sc_p->specific_area * md->use.sur_p->thickness / md->mass_water_aq_x;
      }

			sx->tot_g_moles = sx->moles * (1 + total_g);

      if (sx->moles > 1e-30)
				sx->dg_total_g = sx->dg * sx->tot_g_moles / sx->moles;
      else
				sx->dg_total_g = 0.0;
    }
  }

  if (md->use.gas_p != NULL)
    CalcGasPressures();

  if (md->use.ssa_p != NULL)
    CalcSSFractions();

#ifdef DEBUG_MOHID
	if(d.debug_status)
		fprintf(d.molalities_f, "FIM)----------------\n");
#endif

  return true;
}

bool ModelEngine::ReviseGuesses()
{
	// Revise la's of master species

  int i;
  int iter, max_iter;
	bool repeat, fail;
  LDBLE weight, f;
	String filename;
	Unknown *u;
	Master *m;

  max_iter = 10;
  
#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.gammas_f, "ReviseGuesses-1) called_by_reviseguesses:%i\n", d.gammas_count++);
#endif 	

	Gammas(md->mu_x);

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.reviseguesses_f, "2) %20.20e\n", md->mu_x);
#endif 	

  iter = 0;
  repeat = true;
  fail = false;

  while (repeat)
  {
    iter++;

		if (iter == max_iter + 1)
      fail = true;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.reviseguesses_f, "3) %i %i\n", iter, fail);
#endif 			

		if (iter > 2 * max_iter)
      return false;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.molalities_f, "called_by_revise_guesses) %i\n", d.molalities_count++);
#endif 	

		Molalities(true);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.mbsums_f, "called_by_revise_guesses) %i\n", d.mbsums_count++);
#endif 	

    MBSums();

		if (md->state < REACTION)
      SumSpecies();
    else
		{
			for (i = 0; i < count_unknowns; i++)
			{
				u = (*unknown_list)[i];
				u->sum = u->f;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reviseguesses_f, "4) %i %i %20.20e\n", i, u->number, u->f);
#endif 			
			}
		}

    repeat = false;
    for (i = 0; i < count_unknowns; i++)
    {
			u = (*unknown_list)[i];

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.reviseguesses_f, "5) i:%i number:%i type:%i\n", i, u->number, u->type);
#endif 	
				
			if (u == md->ph_unknown || u == md->pe_unknown)
				continue;

      if (u->type == MB || u->type == CB || u->type == SOLUTION_PHASE_BOUNDARY || u->type == EXCH || u->type == SURFACE)
      {
				if (fabs (u->moles) < 1e-30)
					u->moles = 0;
	
				f = fabs (u->sum);
	
#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reviseguesses_f, "6) %20.20e %20.20e\n", u->moles, f);
#endif 

				if (f == 0 && u->moles == 0)
				{
					m = (*u->master)[0];
					m->s->la = MIN_RELATED_LOG_ACTIVITY;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reviseguesses_f, "7) m:%i s:%i %20.20e\n", m->number, m->s->number, m->s->la);
#endif 

					continue;
				}
				else if (f == 0)
				{
					m = (*u->master)[0];

					repeat = true;
					m->s->la += 5;

					if (m->s->la < -999.0)
						m->s->la = MIN_RELATED_LOG_ACTIVITY;

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reviseguesses_f, "8) repeat:%i m:%i s:%i %20.20e\n", repeat, m->number, m->s->number, m->s->la);
#endif 
				}
				else if (fail && f < 1.5 * fabs (u->moles))
				{
					continue;
				}
				else if (f > 1.5 * fabs (u->moles) || f < 1e-5 * fabs (u->moles))
				{
					m = (*u->master)[0];
					weight = (f < (LDBLE)1e-5 * fabs (u->moles)) ? (LDBLE)0.3 : (LDBLE)1.0;
					
					if (u->moles <= 0)
						m->s->la = MIN_RELATED_LOG_ACTIVITY;
					else
					{
						repeat = true;
						m->s->la += weight * log10 (fabs (u->moles / u->sum));
				  }

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reviseguesses_f, "9) repeat:%i m:%i s:%i weight:%20.15 %20.20e\n", repeat, m->number, m->s->number, weight, m->s->la);
#endif 
				}
      }
      else if (u->type == ALK)
      {
				f = md->total_co2;

#ifdef DEBUG_MOHID
				if (d.debug_status)
					fprintf(d.reviseguesses_f, "10) %20.15 %20.20e\n", f, u->moles);
#endif 

				if (fail && f < 1.5 * fabs (u->moles))
				  continue;

				if (f > 1.5 * fabs (u->moles) || f < 1e-5 * fabs (u->moles))
				{
					m = (*u->master)[0];
					repeat = true;
					weight = (f < (LDBLE)1e-5 * fabs (u->moles)) ? (LDBLE)0.3 : (LDBLE)1.0;
					m->s->la += weight * log10 (fabs (u->moles / u->sum));

#ifdef DEBUG_MOHID
					if (d.debug_status)
						fprintf(d.reviseguesses_f, "11) repeat:%i m:%i s:%i weight:%20.15 %20.15 %20.15 %20.20e\n", repeat, m->number, m->s->number, weight, m->s->la, u->moles, u->sum);
#endif 
				}
      }
    }
  }


  md->mu_x = md->mu_unknown->f * (LDBLE)0.5 / md->mass_water_aq_x;

  if (md->mu_x <= 1e-8)
    md->mu_x = (LDBLE)1e-8;

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.reviseguesses_f, "12) %20.20e %20.20e %20.20e\n", md->mu_x, md->mu_unknown->f, md->mass_water_aq_x);
#endif 	

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.gammas_f, "ReviseGuesses-13) called_by_reviseguesses:%i\n", d.gammas_count++);
#endif 	
	Gammas(md->mu_x);

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.reviseguesses_f, "14) %20.20e\n", md->mu_x);
#endif 	

  return true;
}

bool ModelEngine::InitialSurfaceWater()
{
	//
	// In initial surface calculation, need to calculate
	// mass of water in diffuse layer.
	// diffuse layer water + aqueous solution water = bulk water.
	// Ionic strength is fixed, so diffuse-layer water will not change
	//

  int i;
  LDBLE debye_length, b, r, rd, ddl_limit, rd_limit, fraction, sum_surfs, s;
  LDBLE damp_aq;
	Unknown *u;

  // Debye  length = 1/k = sqrt[eta*eta_zero*R*T/(2*F**2*mu_x*1000)], Dzombak and Morel, p 36
  // 1000 converts kJ to J; 1000 converts Liters to meter**3; debye_length is in meters.
  debye_length = (EPSILON * EPSILON_ZERO * R_KJ_DEG_MOL * (LDBLE)1000.0 * md->tk_x) / ((LDBLE)2. * F_C_MOL * F_C_MOL * md->mu_x * (LDBLE)1000.);
  debye_length = sqrt (debye_length);

	// ddl is at most the fraction ddl_limit of bulk water
	ddl_limit = md->use.sur_p->DDL_limit;

	// Loop through all surface components, calculate each H2O surface (diffuse layer),
	// H2O aq, and H2O bulk (diffuse layers plus aqueous).
  if (md->use.sur_p->debye_lengths > 0)
  {
    sum_surfs = 0.0;

    for (i = 0; i < count_unknowns; i++)
    {
			u = (*unknown_list)[i];
      if (u->type != SURFACE_CB)
        continue;

      sum_surfs += u->surface_charge->specific_area * u->surface_charge->grams;
    }

    rd = debye_length * md->use.sur_p->debye_lengths;
    md->use.sur_p->thickness = rd;

    if (md->state == INITIAL_SURFACE)
    {
      // distribute water over DDL (rd) and free pore (r - rd)
      // find r: free pore (m3) = pi * (r - rd)^2 * L, where L = A / (2*pi*r), A = sum_surfs = sum of the surface areas
			b = (LDBLE)-2 * (rd + md->use.sol_p->mass_water / ((LDBLE)1000 * sum_surfs));
      r = (LDBLE)0.5 * (-b + sqrt (b * b - (LDBLE)4 * rd * rd));
      rd_limit = (1 - sqrt (1 - ddl_limit)) * r;

			// rd should be smaller than r and the ddl limit
			if (rd > rd_limit)
      {
				md->mass_water_surfaces_x = md->use.sol_p->mass_water * ddl_limit / (1 - ddl_limit);
				r = (LDBLE)0.002 * (md->mass_water_surfaces_x + md->use.sol_p->mass_water) / sum_surfs;
				rd_limit = (1 - sqrt (1 - ddl_limit)) * r;
				md->use.sur_p->thickness = rd = rd_limit;
      }
      else
				md->mass_water_surfaces_x = (r * r / pow (r - rd, 2) - 1) * md->use.sol_p->mass_water;
      
			for (i = 0; i < count_unknowns; i++)
      {
				u = (*unknown_list)[i];

				if (u->type != SURFACE_CB)
					continue;

        s = u->surface_charge->specific_area * u->surface_charge->grams;
				u->surface_charge->mass_water = md->mass_water_surfaces_x * s / sum_surfs;
      }
    }
    else
    {
      r = (LDBLE)0.002 * md->mass_water_bulk_x / sum_surfs;
      rd_limit = (1 - sqrt (1 - ddl_limit)) * r;
      
			if (rd > rd_limit)
      {
				md->use.sur_p->thickness = rd = rd_limit;
				fraction = ddl_limit;
      }
      else
				fraction = 1 - pow (r - rd, 2) / (r * r);      
      
			if (md->g_iterations > 10)
				damp_aq = (LDBLE)0.2;
      else if (md->g_iterations > 5)
				damp_aq = (LDBLE)0.5;
			else
				damp_aq = (LDBLE)1.0;

      md->mass_water_surfaces_x = damp_aq * fraction * md->mass_water_bulk_x + (1 - damp_aq) * md->mass_water_surfaces_x;

			for (i = 0; i < count_unknowns; i++)
      {
				u = (*unknown_list)[i];
				if (u->type != SURFACE_CB)
					continue;

        s = u->surface_charge->specific_area * u->surface_charge->grams;
				u->surface_charge->mass_water = md->mass_water_surfaces_x * s / sum_surfs;
      }
    }
  }
  else
  {
		// take constant thickness of, default 1e-8 m (100 Angstroms)
    md->mass_water_surfaces_x = 0.0;

		for (i = 0; i < count_unknowns; i++)
    {
			u = (*unknown_list)[i];

      if (u->type != SURFACE_CB)
				continue;

      u->surface_charge->mass_water = u->surface_charge->specific_area * u->surface_charge->grams * md->use.sur_p->thickness * 1000;
      md->mass_water_surfaces_x += u->surface_charge->mass_water;
    }
  }

  if (md->use.sur_p->type == CD_MUSIC)
    md->mass_water_bulk_x = md->mass_water_aq_x + md->mass_water_surfaces_x;
  else
  {
		// for variable distribution of water over DDL and bulk...
    if (md->state > INITIAL_SURFACE)
      md->mass_water_aq_x = md->mass_water_bulk_x - md->mass_water_surfaces_x;
    else
      md->mass_water_bulk_x = md->mass_water_aq_x + md->mass_water_surfaces_x;
  }

	return true;
}

bool ModelEngine::MBSums()
{
	//
	// Calculates sums of species for calculation of mass balances, charge
	// balance. Also calculates saturation indices for solution_phase_boundaries
	// and pure_phases. After this routine total calcium calculated from all
	// calcium species in solution is stored in x[i]->f.  Also calculates
	// x[i]->sum for some types of unknowns. Uses arrays sum_mb1 and
	// sum_mb1, which are generated in prep and reprep.
	//

  int k;
	Unknown *u;

	// Clear functions in unknowns
	for (k = 0; k < count_unknowns; k++)
  {
		u = (*unknown_list)[k];

    u->f = 0.0;
    u->sum = 0.0;
  }

	// Add terms with coefficients of 1.0
	STCoef *sc_p;
	for (k = 0; k < md->sum_mb1->Count(); k++)
	{
		sc_p = (*md->sum_mb1)[k];
		*sc_p->target += *sc_p->source;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.mbsums_f, "%20.20e %20.20e\n", *sc_p->target, *sc_p->source);
#endif 
	}

	// Add terms with coefficients != 1.0
  for (k = 0; k < md->sum_mb2->Count(); k++)
	{
		sc_p = (*md->sum_mb2)[k];
		*sc_p->target += (*sc_p->source * sc_p->coef);

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.mbsums_f, "%20.20e %20.20e %20.20e\n", *sc_p->target, *sc_p->source, sc_p->coef);
#endif 
	}

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.mbsums_f, "FIM----------\n");
#endif 	
	return true;
}

bool ModelEngine::SwitchBases()
{
	//
	// Check if activity of first master species is predominant among activities of
	// secondary master species included in mass balance.
	//

  int i, j;
  int first;
  bool return_value;
  LDBLE la, la1;

  return_value = false;

	Unknown *u;
	Master *m;
  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];
		
    if (u->type != MB)
      continue;

    first = 0;
    la = (*u->master)[0]->s->la;

    for (j = 1; j < u->master->Count() != NULL; j++)
    {
			m = (*u->master)[j];

      la1 = m->s->lm + m->s->lg;

      if (first == 0 && la1 > la + 10.)
      {
				la = la1;
				first = j;
      }
      else if (first != 0 && la1 > la)
      {
				la = la1;
				first = j;
      }
    }

		if (first != 0)
    {
			u->master->Swap(0, first);

      (*u->master)[0]->in = true;
			(*u->master)[0]->rewrite = false;
      (*u->master)[first]->rewrite = true;
			(*u->master)[first]->in = false;

			return_value = true;
    }
  }

  return (return_value);
}

bool ModelEngine::Reprep()
{
	//
	// If a basis species has been switched, makes new model.
	// Unknowns are not changed, but mass-action equations are
	// rewritten and lists for mass balance and jacobian are regenerated
	//

  int i;
	Master *m;

	// Initialize s, master, and unknown pointers
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];
		if (!m->in)
      continue;

		m->rxn_primary->CopyTo(m->rxn_secondary);
  }

  if (!ResetupMaster()) 
		return false;

	// Set unknown pointers, unknown types, validity checks
  if (!TidyRedox())
		return false;

	// Free arrays built in build_model
	md->s_x->Clear();
	md->sum_mb1->Clear();
	md->sum_mb2->Clear();
	md->sum_jacob0->Clear();
	md->sum_jacob1->Clear();
	md->sum_jacob2->Clear();
	md->sum_delta->Clear();

	// Build model again
  if (!BuildModel ())
		return false;

  md->same_model = false;

  KTemp(md->tc_x);

  return true;
}

Master *ModelEngine::SurfaceGetPSIMaster(String name, int plane)
{
  Master *master_ptr;
  String token;

	if (name.IsEmpty())
    return NULL;

	token = name;
	token += "_psi";

	switch (plane)
  {
  case SURF_PSI:
    break;
  case SURF_PSI1:
    token += "b";
    break;
  case SURF_PSI2:
    token += "d";
    break;
  default:
    throw exception();
  }

	master_ptr = gd->master_list.Search(&token);

  return (master_ptr);
}

LDBLE ModelEngine::Under(LDBLE xval)
{
	//
	// Exponentiate a number, x, but censor large and small numbers
	// log values less than MIN_LM are set to 0
	// log values greater than MAX_LM are set to 10**MAX_LM
	//

  if (xval < -40.)
    return (0.0);

	if (xval > 3.)
    return (1.0e3);

	return (exp (xval * md->LOG_10));
}

bool ModelEngine::SSPrep(LDBLE t, SS *s_s_ptr)
{
/*
  int i, j, k, converged, divisions;
  LDBLE r, rt, ag0, ag1, crit_pt;
  LDBLE xc, tc;
  LDBLE x, x0, x1, xsm1, xsm2, xb1, xb2;
  LDBLE xc1, xc2;
  LDBLE facb1, faca1, spim1, xblm1, acrae, acrael, xliapt, xliapm;
  LDBLE xaly, xaly1, xaly2;
  LDBLE faca, facb, spialy, facal, facbl;
  LDBLE tol;

  if (pr.s_s_assemblage == FALSE)
    print = FALSE;
  tol = 1e-6;
  r = R_KJ_DEG_MOL;
  rt = r * t;
  a0 = s_s_ptr->ag0 / rt;
  a1 = s_s_ptr->ag1 / rt;
  s_s_ptr->a0 = a0;
  s_s_ptr->a1 = a1;
  ag0 = a0 * rt;
  ag1 = a1 * rt;
  kc = exp (k_calc (ss0->phase->rxn->logk, t) * md->LOG_10);
  kb = exp (k_calc (ss1->phase->rxn->logk, t) * md->LOG_10);
  crit_pt = fabs (a0) + fabs (a1);

  s_s_ptr->miscibility = FALSE;
  s_s_ptr->spinodal = FALSE;
  xsm1 = 0.5;
  xsm2 = 0.5;
  xb1 = 0.5;
  xb2 = 0.5;
  xc1 = 0;
  xc2 = 0;

  if (crit_pt >= tol)
  {
    if (fabs (a1) < tol)
    {
      xc = 0.5;
      tc = ag0 / (2 * r);
    }
    else
    {
      xc = 0.5 + (pow ((ag0 * ag0 + 27 * ag1 * ag1), 0.5) - ag0) / (18 * ag1);
      tc = (12 * ag1 * xc - 6 * ag1 + 2 * ag0) * (xc - xc * xc) / r;
    }
    if (print == TRUE)
    {
      sprintf (error_string, "Description of Solid Solution %s",
	       s_s_ptr->name);
      dup_print (error_string, TRUE);
    }
    if (print == TRUE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t                              Temperature: %g kelvin\n",
		  (double) t);
      output_msg (OUTPUT_MESSAGE,
		  "\t                       A0 (dimensionless): %g\n",
		  (double) a0);
      output_msg (OUTPUT_MESSAGE,
		  "\t                       A1 (dimensionless): %g\n",
		  (double) a1);
      output_msg (OUTPUT_MESSAGE,
		  "\t                              A0 (kJ/mol): %g\n",
		  (double) ag0);
      output_msg (OUTPUT_MESSAGE,
		  "\t                              A1 (kJ/mol): %g\n\n",
		  (double) ag1);
    }
    if (xc < 0 || xc > 1)
    {
      if (print == TRUE)
	output_msg (OUTPUT_MESSAGE,
		    "No miscibility gap above 0 degrees kelvin.\n");
    }
    else
    {
      if (print == TRUE)
      {
	output_msg (OUTPUT_MESSAGE,
		    "\t    Critical mole-fraction of component 2: %g\n",
		    (double) xc);
	output_msg (OUTPUT_MESSAGE,
		    "\t                     Critical temperature: %g kelvin\n",
		    (double) tc);
	output_msg (OUTPUT_MESSAGE,
		    "\n(The critical temperature calculation assumes that the Guggenheim model\ndefined at %g kelvin md->is valid at the critical temperature.)\n\n\n",
		    (double) t);
      }
    }

    if (tc >= t)
    {
      x0 = 0;
      x1 = 1;
      if (scan (f_spinodal, &x0, &x1) == TRUE)
      {

	xsm1 = halve (f_spinodal, x0, x1, tol);
	s_s_ptr->spinodal = TRUE;

	x0 = x1;
	x1 = 1;
	if (scan (f_spinodal, &x0, &x1) == TRUE)
	{
	  xsm2 = halve (f_spinodal, x0, x1, tol);
	}
	else
	{
	  error_msg ("Failed to find second spinodal point.", STOP);
	}
      }
    }
  }
  if (s_s_ptr->spinodal == TRUE)
  {
    if (print == TRUE)
      output_msg (OUTPUT_MESSAGE,
		  "\t Spinodal-gap mole fractions, component 2: %g\t%g\n",
		  (double) xsm1, (double) xsm2);
    converged = FALSE;
    if (converged == FALSE)
    {
      for (i = 1; i < 3; i++)
      {
	divisions = (int) pow (10., i);
	for (j = 0; j < divisions; j++)
	{
	  for (k = divisions; k > 0; k--)
	  {
	    xc1 = (LDBLE) j / divisions + 0.001;
	    xc2 = (LDBLE) k / divisions;
	    converged = solve_misc (&xc1, &xc2, tol);
	    if (converged == TRUE)
	      break;
	  }
	  if (converged == TRUE)
	    break;
	}
	if (converged == TRUE)
	  break;
      }
    }
    if (converged == FALSE)
    {
      error_msg ("Failed to find miscibility gap.", STOP);
    }
    s_s_ptr->miscibility = TRUE;
    if (xc1 < xc2)
    {
      xb1 = 1 - xc2;
      xb2 = 1 - xc1;
      xc1 = 1 - xb1;
      xc2 = 1 - xb2;
    }
    else
    {
      xb1 = 1 - xc1;
      xb2 = 1 - xc2;
    }
    facb1 = kb * xb1 * exp (xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
    faca1 = kc * xc1 * exp (xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
    spim1 = log10 (faca1 + facb1);
    xblm1 = 1. / (1. + faca1 / facb1);
    acrae = facb1 / faca1;
    acrael = log10 (acrae);
    xliapt = log10 (facb1);
    xliapm = log10 (faca1);

    if (print == TRUE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t   Miscibility-gap fractions, component 2: %g\t%g\n",
		  (double) xb1, (double) xb2);
      fprintf(f, "\n\t\t\tEutectic Point Calculations\n\n");
      output_msg (OUTPUT_MESSAGE,
		  "\t     Aqueous activity ratio (comp2/comp1): %g\n",
		  (double) acrae);
      output_msg (OUTPUT_MESSAGE,
		  "\t Log aqueous activity ratio (comp2/comp1): %g\n",
		  (double) acrael);
      output_msg (OUTPUT_MESSAGE,
		  "\t Aqueous activity fraction of component 2: %g\n",
		  (double) xblm1);
      output_msg (OUTPUT_MESSAGE,
		  "\t                    Log IAP (component 2): %g\n",
		  (double) xliapt);
      output_msg (OUTPUT_MESSAGE,
		  "\t                    Log IAP (component 1): %g\n",
		  (double) xliapm);
      output_msg (OUTPUT_MESSAGE,
		  "\t                               Log Sum Pi: %g\n",
		  (double) spim1);
    }
    s_s_ptr->tk = t;
    s_s_ptr->xb1 = xb1;
    s_s_ptr->xb2 = xb2;
  }
  xaly = -1.0;
  x = a0 * a0 + 3 * a1 * a1 + 6 * a1 * log (kb / kc);
  if (x > 0)
  {
    if (fabs (x - a0 * a0) >= tol)
    {
      xaly1 = (-(a0 - 3 * a1) + pow (x, 0.5)) / (6 * a1);
      xaly2 = (-(a0 - 3 * a1) - pow (x, 0.5)) / (6 * a1);
      if (xaly1 >= 0 && xaly1 <= 1)
      {
	xaly = xaly1;
      }
      if (xaly2 >= 0 && xaly2 <= 1)
      {
	xaly = xaly2;
      }
    }
    else
    {
      xaly = 0.5 + log (kb / kc) / (2 * a0);
    }
    if (xaly > 0 && xaly < 1)
    {
      faca = kc * (1 - xaly) * exp (xaly * xaly * (a0 - a1 * (3 - 4 * xaly)));
      facb =
	kb * xaly * exp ((1 - xaly) * (1 - xaly) *
			 (a0 + a1 * (4 * xaly - 1.0)));
      spialy = log10 (faca + facb);
      facal = log10 (faca);
      facbl = log10 (facb);
      if (xaly > xb1 && xaly < xb2)
      {
	if (print == TRUE)
	  output_msg (OUTPUT_MESSAGE,
		      "\nLocal minimum in the solidus curve coresponding to a maximum\nin the minimum stoichiometric saturation curve.\n\n");
      }
      else
      {
	if (print == TRUE)
	  fprintf(f, "\n\t\t\tAlyotropic Point\n\n");
      }
      if (print == TRUE)
      {
	output_msg (OUTPUT_MESSAGE,
		    "\t       Solid mole fraction of component 2: %g\n",
		    (double) xaly);
	output_msg (OUTPUT_MESSAGE,
		    "\t                    Log IAP (component 2): %g\n",
		    (double) facbl);
	output_msg (OUTPUT_MESSAGE,
		    "\t                    Log IAP (component 1): %g\n",
		    (double) facal);
	output_msg (OUTPUT_MESSAGE,
		    "\t                               Log Sum Pi: %g\n",
		    (double) spialy);
      }
    }
  }
	*/
  return true;
}

bool ModelEngine::InitialGuesses()
{
	//
	// Make initial guesses for activities of master species and ionic strength
	//

  int i;

	Solution *sol_p = md->use.sol_p;

  md->mu_x = gd->s_hplus->moles + exp((sol_p->ph - (LDBLE)14.) * md->LOG_10) * md->mass_water_aq_x;
  md->mu_x /= md->mass_water_aq_x;
	//d.PrintSpeciesInfoToFile("Initialguesses-1", md->species_info_list, md, gd);
  gd->s_h2o->la = 0.0;

	Unknown *u;
	Master *m;
	for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];		

    if (u == md->ph_unknown || u == md->pe_unknown)
      continue;

    if (u->type < CB)
    {
			m = (*u->master)[0];

      md->mu_x += u->moles / md->mass_water_aq_x * (LDBLE)0.5 * m->s->z * m->s->z;
			//d.PrintSpeciesInfoToFile("InitialGuesses-2", md->species_info_list, md, gd);
      m->s->la = log10 (u->moles / md->mass_water_aq_x);
    }
    else if (u->type == CB || u->type == SOLUTION_PHASE_BOUNDARY)
		{
			m = (*u->master)[0];

      m->s->la = log10 ((LDBLE)0.001 * u->moles / md->mass_water_aq_x);
		}
    else if (u->type == EXCH)
    {
			m = (*u->master)[0];

      if (u->moles <= 0)
				m->s->la = MIN_RELATED_LOG_ACTIVITY;
      else
				m->s->la = log10 (u->moles);
    }
    else if (u->type == SURFACE)
    {
			m = (*u->master)[0];

      if (u->moles <= 0)
				m->s->la = MIN_RELATED_LOG_ACTIVITY;
      else
				m->s->la = log10 ((LDBLE)0.1 * u->moles);
    }
    else if (u->type == SURFACE_CB)
		{
			m = (*u->master)[0];

      m->s->la = 0.0;
		}
  }

  return true;
}

bool ModelEngine::CalculateValues()
{
	/*
  int j;
  struct calculate_value *calculate_value_ptr;
  struct isotope_ratio *isotope_ratio_ptr;
  struct isotope_alpha *isotope_alpha_ptr;
  struct master_isotope *master_isotope_ptr;
  char command[] = "run";

  for (j = 0; j < count_calculate_value; j++)
  {
    calculate_value[j]->calculated = FALSE;
    calculate_value[j]->value = MISSING;
  }

  for (j = 0; j < count_calculate_value; j++)
  {
    calculate_value_ptr = calculate_value[j];
    rate_moles = NAN;
    if (calculate_value_ptr->new_def == TRUE)
    {
      if (basic_compile
	  (calculate_value[j]->commands, &calculate_value[j]->linebase,
	   &calculate_value[j]->varbase, &calculate_value[j]->loopbase) != 0)
      {
	sprintf (error_string, "Fatal Basic error in CALCULATE_VALUES %s.",
		 calculate_value[j]->name);
	error_msg (error_string, STOP);
      }
      calculate_value_ptr->new_def = FALSE;
    }
    if (basic_run
	(command, calculate_value[j]->linebase, calculate_value[j]->varbase,
	 calculate_value[j]->loopbase) != 0)
    {
      sprintf (error_string, "Fatal Basic error in calculate_value %s.",
	       calculate_value[j]->name);
      error_msg (error_string, STOP);
    }
    if (rate_moles == NAN)
    {
      sprintf (error_string, "Calculated value not SAVE'd for %s.",
	       calculate_value[j]->name);
      error_msg (error_string, STOP);
    }
    else
    {
      calculate_value[j]->calculated = TRUE;
      calculate_value[j]->value = rate_moles;
    }
  }
  for (j = 0; j < count_isotope_ratio; j++)
  {
    isotope_ratio_ptr = isotope_ratio[j];
    master_isotope_ptr =
      master_isotope_search (isotope_ratio_ptr->isotope_name);
    calculate_value_ptr = calculate_value_search (isotope_ratio_ptr->name);
    
    if (calculate_value_ptr->value == MISSING)
    {
      isotope_ratio_ptr->ratio = MISSING;
      isotope_ratio_ptr->converted_ratio = MISSING;
    }
    else
    {
      isotope_ratio_ptr->ratio = calculate_value_ptr->value;
      isotope_ratio_ptr->converted_ratio =
	convert_isotope (master_isotope_ptr, calculate_value_ptr->value);
    }
  }
  for (j = 0; j < count_isotope_alpha; j++)
  {
    isotope_alpha_ptr = isotope_alpha[j];
    calculate_value_ptr = calculate_value_search (isotope_alpha_ptr->name);
   
    if (calculate_value_ptr->value == MISSING)
    {
      isotope_alpha_ptr->value = MISSING;
    }
    else
    {
      isotope_alpha_ptr->value = calculate_value_ptr->value;
    }
  }
	*/
  return true;
}

bool ModelEngine::IneqInit(int max_row_count, int max_column_count)
{
	if (normal_max() <= 0)
		NormalNewMax(count_unknowns);

	if (ineq_array_max() <= 0)
		IneqArrayNewMax(max_row_count * max_column_count);

	if (back_eq_max() <= 0)
		BackEqNewMax(max_row_count);

	if (zero_max() <= 0)
		ZeroNewMax(max_row_count);

	if (res_max() <= 0)
		ResNewMax(max_row_count);

	if (delta1_max() <= 0)
		Delta1NewMax(max_column_count);

	if (cu_max() <= 0)
		CuNewMax(3 * max_row_count);

	if (iu_max() <= 0)
		IuNewMax(3 * max_row_count);

	if (is_max() <= 0)
		IsNewMax(3 * max_row_count);

	return true;
}

bool ModelEngine::CL1(int k, int l, int m, int n,
										  int nklmd, int n2d,
										  LDBLE * q,
										  int *kode, LDBLE toler,
										  int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
										  LDBLE *cu, int *iu, int *s, int check)
{
  union double_or_int
  {
    int ival;
    LDBLE dval;
  } *q2;

  static int nklm;
  static LDBLE xmin, xmax;
  static int iout, i, j;
  static LDBLE z;
  static int maxit, n1, n2;
  static LDBLE pivot;
  static int ia, ii, kk, in, nk, js;
  static LDBLE sn;
  static int iphase, kforce;
  static LDBLE zu, zv;
  static LDBLE tpivot;
  static int klm, jmn, nkl, jpn;
  static LDBLE cuv, sum;
  static int klm1;
  int q_dim, cu_dim;
  int kode_arg;
  LDBLE check_toler;

  zv = 0;
  kode_arg = *kode;
  CL1Space(check, n2d, k+l+m, nklmd);

  q_dim = n2d;
  q2 = (union double_or_int *) q;
  cu_dim = nklmd;

  maxit = *iter;
  n1 = n + 1;
  n2 = n + 2;
  nk = n + k;
  nkl = nk + l;
  klm = k + l + m;
  klm1 = klm + 1;
  nklm = n + klm;
  kforce = 1;
  *iter = 0;
  js = 0;
  ia = -1;

#ifdef DEBUG_MOHID
	/*
	if (d.debug_status)
		fprintf(cl1_f, "CL1-1) %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
										kode_arg,
										q_dim,
										cu_dim,
										max_it,
										n1,
										n2,
										nk,
										nlk,
										klm,
										klm1,
										nklm,
										kforce,
										*iter,
										js,
										ia);	
										*/
#endif	
	
	
  for (j = 0; j < n; ++j)
	{
    q2[klm1 * q_dim + j].ival = j + 1;
		
#ifdef DEBUG_MOHID
		/*
		if (d.debug_status)
			fprintf(cl1_f, "CL1-2) %d %d", 
											klm1 * q_dim + j,
											q2[klm1 * q_dim + j].ival);										
											*/
#endif			
	}

  for (i = 0; i < klm; ++i)
  {
    q2[i * q_dim + n1].ival = n + i + 1;
		
		if (q2[i * q_dim + n].dval < 0.)
    {
      for (j = 0; j < n1; ++j)
				q2[i * q_dim + j].dval = -q2[i * q_dim + j].dval;

			q2[i * q_dim + n1].ival = -q2[i * q_dim + n1].ival;
    }
  }

	
  iphase = 2;

  memcpy ((void *) &cu[0], (void *) &scratch_e[0], (size_t) nklm * sizeof (LDBLE));

  for (j = 0; j < nklm; ++j)
    iu[j] = 0;

  if (l != 0)
  {
    for (j = nk; j < nkl; ++j)
    {
      cu[j] = 1.;
      iu[j] = 1;
    }

    iphase = 1;
  }

  memcpy ((void *) &cu[cu_dim], (void *) &cu[0], (size_t) nklm * sizeof (LDBLE));
  memcpy ((void *) &iu[cu_dim], (void *) &iu[0], (size_t) nklm * sizeof (int));

  if (m != 0)
  {
    for (j = nkl; j < nklm; ++j)
    {
      cu[cu_dim + j] = 1.;
      iu[cu_dim + j] = 1;

      jmn = j - n;

      if (q2[jmn * q_dim + n1].ival < 0)
				iphase = 1;
    }
  }

	if (*kode != 0)
  {
    for (j = 0; j < n; ++j)
    {
      if (x[j] < 0.)
      {
				cu[j] = 1.;
				iu[j] = 1;
      }
      else if (x[j] > 0.)
      {
				cu[cu_dim + j] = 1.;
				iu[cu_dim + j] = 1;
      }
    }

		for (j = 0; j < k; ++j)
    {
      jpn = j + n;

      if (res[j] < 0.)
      {
				cu[jpn] = 1.;
				iu[jpn] = 1;

				if (q2[j * q_dim + n1].ival > 0)
					iphase = 1;
      }
      else if (res[j] > 0.)
      {
				cu[cu_dim + jpn] = 1.;
				iu[cu_dim + jpn] = 1;

				if (q2[j * q_dim + n1].ival < 0)
					iphase = 1;
      }
    }
  }

	if (iphase == 2)
  {
    goto L500;
  }

L160:

	for (j = js; j < n1; ++j)
  {
    sum = 0.;

		for (i = 0; i < klm; ++i)
    {
      ii = q2[i * q_dim + n1].ival;

			if (ii < 0)
				z = cu[cu_dim - ii - 1];
      else
				z = cu[ii - 1];

			sum += q2[i * q_dim + j].dval * z;
    }

    q2[klm * q_dim + j].dval = sum;
  }

  for (j = js; j < n; ++j)
  {
    ii = q2[klm1 * q_dim + j].ival;
  
		if (ii < 0)
      z = cu[cu_dim - ii - 1];
    else
      z = cu[ii - 1];

    q2[klm * q_dim + j].dval -= z;
  }

L240:

	xmax = 0.;

	if (js >= n)
    goto L490;			

  for (j = js; j < n; ++j)
  {
    zu = q2[klm * q_dim + j].dval;
    ii = q2[klm1 * q_dim + j].ival;
  
		if (ii > 0)
      zv = -zu - cu[ii - 1] - cu[cu_dim + ii - 1];
    else
		{
      ii = -ii;
      zv = zu;
			zu = -zu - cu[ii - 1] - cu[cu_dim + ii - 1];
    }

		if (kforce == 1 && ii > n)
      continue;

		if (iu[ii - 1] != 1 && zu > xmax)
    {
      xmax = zu;
      in = j;
    }

		if (iu[cu_dim + ii - 1] != 1 && zv > xmax)
    {
      xmax = zv;
      in = j;
    }
  }

	if (xmax <= toler)
    goto L490;

	if (q2[klm * q_dim + in].dval != xmax)
  {
    for (i = 0; i < klm1; ++i)
      q2[i * q_dim + in].dval = -q2[i * q_dim + in].dval;

		q2[klm1 * q_dim + in].ival = -q2[klm1 * q_dim + in].ival;
    q2[klm * q_dim + in].dval = xmax;
  }

	if (iphase != 1 && ia != -1)
  {
    xmax = 0.;

		for (i = 0; i <= ia; ++i)
    {
      z = fabs (q2[i * q_dim + in].dval);

			if (z > xmax)
      {
				xmax = z;
				iout = i;
      }
    }

    if (xmax > toler)
    {
      memcpy ((void *) &scratch_e[0], (void *) &(q2[ia * q_dim]), (size_t) n2 * sizeof (LDBLE));
      memcpy ((void *) &(q2[ia * q_dim]), (void *) &(q2[iout * q_dim]), (size_t) n2 * sizeof (LDBLE));
      memcpy ((void *) &(q2[iout * q_dim]), (void *) &scratch_e[0], (size_t) n2 * sizeof (LDBLE));

			iout = ia;
      --ia;
      pivot = q2[iout * q_dim + in].dval;
      goto L420;
    }
  }

  kk = -1;

  for (i = 0; i < klm; ++i)
  {
    z = q2[i * q_dim + in].dval;
    if (z > toler)
    {
      ++kk;
      res[kk] = q2[i * q_dim + n].dval / z;
      s[kk] = i;
    }
  }

L350:

	if (kk < 0)
  {
    *kode = 2;
    goto L590;
  }

  xmin = res[0];
  iout = s[0];
  j = 0;

  if (kk != 0)
  {
    for (i = 1; i <= kk; ++i)
    {
      if (res[i] < xmin)
      {
				j = i;
				xmin = res[i];
				iout = s[i];
      }
    }

    res[j] = res[kk];
    s[j] = s[kk];
  }

  --kk;
  pivot = q2[iout * q_dim + in].dval;
  ii = q2[iout * q_dim + n1].ival;

	if (iphase != 1)
  {
    if (ii < 0 && iu[-ii - 1] == 1)
			goto L420;
    else if (iu[cu_dim + ii - 1] == 1)
			goto L420;
  }

  ii = abs (ii);
  cuv = cu[ii - 1] + cu[cu_dim + ii - 1];
  
	if (q2[klm * q_dim + in].dval - pivot * cuv > toler)
  {
    for (j = js; j < n1; ++j)
    {
      z = q2[iout * q_dim + j].dval;
      q2[klm * q_dim + j].dval -= z * cuv;
      q2[iout * q_dim + j].dval = -z;
    }

		q2[iout * q_dim + n1].ival = -q2[iout * q_dim + n1].ival;

		goto L350;
  }

L420:

	if (*iter >= maxit)
  {
    *kode = 3;

		goto L590;
  }

	++(*iter);
  for (j = js; j < n1; ++j)
    if (j != in)
      q2[iout * q_dim + j].dval /= pivot;

  for (j = js; j < n1; ++j)
  {
    if (j != in)
    {
      z = -q2[iout * q_dim + j].dval;

			for (i = 0; i < klm1; ++i)
				if (i != iout)
					q2[i * q_dim + j].dval += z * q2[i * q_dim + in].dval;
    }
  }

	tpivot = -pivot;

	for (i = 0; i < klm1; ++i)
    if (i != iout)
      q2[i * q_dim + in].dval /= tpivot;

  q2[iout * q_dim + in].dval = (LDBLE)1. / pivot;
  ii = q2[iout * q_dim + n1].ival;
  q2[iout * q_dim + n1].ival = q2[klm1 * q_dim + in].ival;
  q2[klm1 * q_dim + in].ival = ii;
  ii = abs (ii);

	if (iu[ii - 1] == 0 || iu[cu_dim + ii - 1] == 0)
    goto L240;

  for (i = 0; i < klm1; ++i)
  {
    z = q2[i * q_dim + in].dval;
    q2[i * q_dim + in].dval = q2[i * q_dim + js].dval;
    q2[i * q_dim + js].dval = z;
  }

	i = q2[klm1 * q_dim + in].ival;
  q2[klm1 * q_dim + in].ival = q2[klm1 * q_dim + js].ival;
  q2[klm1 * q_dim + js].ival = i;
  ++js;

	goto L240;

L490:

	if (kforce == 0)
  {
    if (iphase == 1)
    {
      if (q2[klm * q_dim + n].dval <= toler)
				goto L500;

			*kode = 1;

			goto L590;
    }

		*kode = 0;

		goto L590;
  }

	if (iphase != 1 || q2[klm * q_dim + n].dval > toler)
  {
    kforce = 0;

		goto L240;
  }

L500:

  iphase = 2;

	for (j = 0; j < nklm; ++j)
    cu[j] = 0.;

  for (j = n; j < nk; ++j)
    cu[j] = 1.;

	memcpy ((void *) &cu[cu_dim], (void *) &cu[0], (size_t) nklm * sizeof (LDBLE));

  for (i = 0; i < klm; ++i)
  {
    ii = q2[i * q_dim + n1].ival;
    if (ii <= 0)
    {
      if (iu[cu_dim - ii - 1] == 0)
				continue;
      
      cu[cu_dim - ii - 1] = 0.;
    }
    else
    {
      if (iu[ii - 1] == 0)
				continue;
      
      cu[ii - 1] = 0.;
    }

    ++ia;

    memcpy ((void *) &scratch_e[0], (void *) &(q2[ia * q_dim]), (size_t) n2 * sizeof (LDBLE));
    memcpy ((void *) &(q2[ia * q_dim]), (void *) &(q2[i * q_dim]), (size_t) n2 * sizeof (LDBLE));
    memcpy ((void *) &(q2[i * q_dim]), (void *) &scratch_e[0], (size_t) n2 * sizeof (LDBLE));
  }

	goto L160;

L590:

	sum = 0.;

	for (j = 0; j < n; ++j)
    x[j] = 0.;

	for (i = 0; i < klm; ++i)
    res[i] = 0.;

	for (i = 0; i < klm; ++i)
  {
    ii = q2[i * q_dim + n1].ival;
    sn = 1.;

		if (ii < 0)
    {
      ii = -ii;
      sn = -1.;
    }

		if (ii <= n)
      x[ii - 1] = sn * q2[i * q_dim + n].dval;
    else
    {
      res[ii - n - 1] = sn * q2[i * q_dim + n].dval;

			if (ii >= n1 && ii <= nk)
				sum += q2[i * q_dim + n].dval;
    }
  }

  *error = sum;

  if ((check == 1) && (*kode == 0))
  {
    check_toler = (LDBLE)10. * toler;

		if (kode_arg == 1)
    {
      for (i = 0; i < k; i++)
      {
				if (res_arg_e[i] < 0.0 && res[i] > check_toler)
					*kode = 1;
				else if (res_arg_e[i] > 0.0 && res[i] < -check_toler)
					*kode = 1;
      }
    }

    for (i = k; i < k + l; i++)
      if (fabs (res[i]) > check_toler)
				*kode = 1;

    for (i = k + l; i < k + l + m; i++)
      if (res[i] < -check_toler)
				*kode = 1;

    if (kode_arg == 1)
    {
      for (i = 0; i < n; i++)
      {
				if (x_arg_e[i] < 0.0 && x[i] > check_toler)
					*kode = 1;
				else if (x_arg_e[i] > 0.0 && x[i] < -check_toler)
					*kode = 1;
      }
    }

		if (*kode == 1)
			return false;
  }

	//arr.PrintToFile("array_cl1.txt");

	return true;
}

bool ModelEngine::CL1Space(int check, int n2d, int klm, int nklmd)
{
	if (check == 1)
  {
		XArgENewMax(max(n2d, x_arg_e_max()));
		Fill(x_arg_e, 0.0, x_arg_e_max());

		ResArgENewMax(max(klm, res_arg_e_max()));
		Fill(res_arg_e, 0.0, res_arg_e_max());
	}

	ScratchENewMax(max(nklmd, scratch_e_max()));
	Fill(scratch_e, 0.0, scratch_e_max());

	return true;
}

bool ModelEngine::CalcAllDonnan()
{
  int i, j, k;
  int count_g, count_charge;
	bool converge;
  //char name[MAX_LENGTH];
  LDBLE new_g, f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;
  LDBLE new_g2, f_psi2, surf_chrg_eq2, psi_avg2, dif;
  LDBLE cz, cm, cp;

	if (md->use.sur_p == NULL)
    return true;

	if (md->use.sur_p->type == CD_MUSIC)
    return (CalcAllDonnanMusic ());

  f_sinh = sqrt ((LDBLE)8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * md->tk_x * md->mu_x);
  cz = cm = 1.0;
  cp = 1.0;

  converge = true;
  count_charge = 0;

	Unknown *x;
	ChargeGroup *cg_p;
	Species *s_x;
	Master *master;
	SurfaceDiffLayer *sdl_p;
	SpeciesDiffLayer *spdl_p;
  for (j = 0; j < count_unknowns; j++)
  {
		x = (*unknown_list)[j];
		master = (*x->master)[0];

    if (x->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = x->surface_charge;

		count_g = x->surface_charge->g->Count();
		
    for (i = 0; i < count_g; i++)
		{
      cg_p = gd->charge_group[i];
			cg_p->eq = 0.0;
		}

    for (i = 0; i < md->s_x->Count(); i++)
    {
			s_x = (*md->s_x)[i];

      if (s_x->type > HPLUS)
				continue;

      for (k = 0; k < count_g; k++)
      {
				cg_p = gd->charge_group[k];

				if (Equal(cg_p->z, s_x->z, md->G_TOL))
				{
					cg_p->eq += s_x->z * s_x->moles;
					break;
				}
      }
    }
   
    A_surf = x->surface_charge->specific_area * x->surface_charge->grams;
    f_psi = master->s->la * md->LOG_10;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;
   
    dif = (LDBLE)1e-5;
    f_psi2 = f_psi + dif;
    surf_chrg_eq2 = A_surf * f_sinh * sinh (f_psi2) / F_C_MOL;

    psi_avg = CalcPSIAvg(surf_chrg_eq);
    psi_avg2 = CalcPSIAvg(surf_chrg_eq2);

    ratio_aq = md->surface_charge_ptr->mass_water / md->mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
			sdl_p = (*x->surface_charge->g)[k];
			cg_p = gd->charge_group[k];
     	sdl_p->charge = cg_p->z;
      new_g = ratio_aq * (exp (-cg_p->z * psi_avg) - 1);

      if (md->use.sur_p->only_counter_ions && ((surf_chrg_eq < 0 && cg_p->z < 0) || (surf_chrg_eq > 0 && cg_p->z > 0)))
				new_g = -ratio_aq;

      if (new_g <= -ratio_aq)
				new_g = -ratio_aq + md->G_TOL * (LDBLE)1e-3;

      new_g2 = ratio_aq * (exp (-cg_p->z * psi_avg2) - 1);

      if (md->use.sur_p->only_counter_ions && ((surf_chrg_eq < 0 && cg_p->z < 0) || (surf_chrg_eq > 0 && cg_p->z > 0)))
				new_g2 = -ratio_aq;

      if (new_g2 <= -ratio_aq)
				new_g2 = -ratio_aq + md->G_TOL * (LDBLE)1e-3;

      if (fabs (new_g) >= 1)
      {
				if (fabs ((new_g - sdl_p->g) / new_g) > md->convergence_tolerance)
					converge = false;
      }
      else
      {
				if (fabs (new_g - sdl_p->g) > md->convergence_tolerance)
				  converge = false;
      }

      sdl_p->g = new_g;

      if (new_g != 0)
				sdl_p->dg = (new_g2 - new_g) / dif;
      else
				sdl_p->dg = -cg_p->z;

			for (i = 0; i < md->s_x->Count(); i++)
      {
				s_x = (*md->s_x)[i];
				spdl_p = (*s_x->diff_layer)[count_charge];

				if (Equal (cg_p->z, s_x->z, md->G_TOL))
				{
					spdl_p->charge = x->surface_charge;
					spdl_p->count_g = k;
				}
      }
    }

		count_charge++;
  }
  
	return (converge);
}

bool ModelEngine::CalcGasPressures()
{
  int i, j;
  LDBLE lp;
  ReactionToken *rxn_ptr;
  GasComp *gas_comp_ptr;
  Phase *phase_ptr;
 
	// moles and partial pressures for gases
  if (md->use.gas_p == NULL)
    return true;

  if (md->use.gas_p->type == VOLUME)
  {
    md->use.gas_p->total_p = 0;
    md->use.gas_p->total_moles = 0;
  }

	for (i = 0; i < md->use.gas_p->comps->Count(); i++)
  {
    gas_comp_ptr = (*md->use.gas_p->comps)[i];
    phase_ptr = (*md->use.gas_p->comps)[i]->phase;

    if (phase_ptr->in)
    {
      lp = -phase_ptr->lk;

			for(j = 1; j < phase_ptr->rxn_x->token_list->Count(); j++)
      {
				rxn_ptr = phase_ptr->rxn_x->token_list->Element(j);
				
				lp += rxn_ptr->s->la * rxn_ptr->coef;
      }

      gas_comp_ptr->phase->p_soln_x = exp (lp * md->LOG_10);

      if (md->use.gas_p->type == PRESSURE)
      {
				gas_comp_ptr->phase->moles_x = gas_comp_ptr->phase->p_soln_x * md->gas_unknown->moles / md->gas_unknown->gas_phase->total_p;
				gas_comp_ptr->phase->fraction_x = gas_comp_ptr->phase->moles_x / md->gas_unknown->moles;
      }
      else
      {
				gas_comp_ptr->phase->moles_x = gas_comp_ptr->phase->p_soln_x * md->use.gas_p->volume / (R_LITER_ATM * md->tk_x);
				md->use.gas_p->total_p += gas_comp_ptr->phase->p_soln_x;
				md->use.gas_p->total_moles += gas_comp_ptr->phase->moles_x;
      }
    }
    else
    {
      gas_comp_ptr->phase->moles_x = 0;
      gas_comp_ptr->phase->fraction_x = 0;
    }
  }

  return true;
}

bool ModelEngine::CalcSSFractions()
{
  int i, k;
  LDBLE moles, n_tot;
  SS *s_s_ptr;
	SSComp *ssc_p;

	// moles and lambdas for solid solutions
  if (md->s_s_unknown == NULL)
    return true;
		
	// Calculate mole fractions and log lambda and derivative factors
	if (md->use.ssa_p == NULL)
    return true;
		
	for (i = 0; i < md->use.ssa_p->ss_list->Count(); i++)
  {
    s_s_ptr = (*md->use.ssa_p->ss_list)[i];
    n_tot = 0;
		
		for (k = 0; k < s_s_ptr->comps_list->Count(); k++)
    {
			ssc_p = (*s_s_ptr->comps_list)[k];
      moles = ssc_p->moles;
			
      if (ssc_p->moles < 0)
      {
				moles = (LDBLE)MIN_TOTAL_SS;
				ssc_p->initial_moles = moles;
      }
			
      n_tot += moles;
    }
		
    s_s_ptr->total_moles = n_tot;
		
		for (k = 0; k < s_s_ptr->comps_list->Count(); k++)
    {
 			ssc_p = (*s_s_ptr->comps_list)[k];
			moles = ssc_p->moles;
			
      if (ssc_p->moles < 0)
				moles = (LDBLE)MIN_TOTAL_SS;
			
      ssc_p->fraction_x = moles / n_tot;
      ssc_p->log10_fraction_x = log10 (moles / n_tot);
			
      // all mb and jacobian items must be in x or phase to be static between models
      ssc_p->phase->log10_fraction_x = ssc_p->log10_fraction_x;
    }
		
    if (s_s_ptr->a0 != 0.0 || s_s_ptr->a1 != 0)
      SSBinary (s_s_ptr);	
    else
      SSIdeal (s_s_ptr);
  }
	
  return true;
}

bool ModelEngine::ResetupMaster()
{
  int i, j;
  Master *master_ptr, *master_ptr0;

  for (i = 0; i < count_unknowns; i++)
  {
    if ((*unknown_list)[i]->type != MB)
      continue;

    master_ptr0 = (*(*unknown_list)[i]->master)[0];

    for (j = 0; j < (*unknown_list)[i]->master->Count(); j++)
    {
			master_ptr = (*(*unknown_list)[i]->master)[j];

      if (j == 0 && master_ptr->s->primary == NULL)
				master_ptr->s->rxn_s->CopyTo(master_ptr->rxn_secondary);
      else if (master_ptr0->s->primary == NULL)
			{
				r_temp->Reset();
				RewriteMasterToSecondary(master_ptr, master_ptr0);
				master_ptr->rxn_secondary->Reset();
				r_temp->CopyReactions(master_ptr->rxn_secondary);
			}
    }
  }

  return true;
}

bool ModelEngine::SetupMasterReaction(ListOfPointers<Master> *m, Reaction **r)
{
  int j;
  Master *master_ptr, *master_ptr0;

	master_ptr0 = (*m)[0];
  for (j = 0; j < m->Count(); j++)
  {
		master_ptr = (*m)[j];

    if (master_ptr->s == gd->s_h2o)
			return false;

		if ((master_ptr->in || master_ptr->rewrite) && (master_ptr->s != gd->s_eminus && master_ptr->s != gd->s_hplus))
			return false;

    if (j == 0)
    {
      master_ptr->in = true;
			master_ptr->rewrite = false;

      if (master_ptr->s->primary == NULL)
				master_ptr->s->rxn_s->CopyTo(master_ptr->rxn_secondary);
    }
    else
    {
      master_ptr->rewrite = true;

      if (master_ptr0->s->primary == NULL)
      {
				r_temp->Reset();
				RewriteMasterToSecondary(master_ptr, master_ptr0);
				master_ptr->rxn_secondary->Reset();
				r_temp->CopyReactions(master_ptr->rxn_secondary);
      }
    }

    *master_ptr->pe_rxn = *r;
  }

  return true;
}

void ModelEngine::PrintTotals(FILE *f)
{
  int i;
	bool pure_water;
  LDBLE EC;

	pure_water = true;

  PrintCentered(f, "Solution composition");
  fprintf(f, "\t%-15s%12s%12s\n\n", "Elements", "Molality", "Moles");
  
	for (i = 0; i < count_unknowns; i++)
  {
    if ((*unknown_list)[i] == md->alkalinity_unknown)
    {
			fprintf(f, "\t%-15s%12.3e%12.3e\n", (*unknown_list)[i]->total->name.CharPtr(), (double) ((*unknown_list)[i]->f / md->mass_water_aq_x), (double) (*unknown_list)[i]->f);
      pure_water = false;
    }
    if ((*unknown_list)[i] == md->ph_unknown)
      continue;
    if ((*unknown_list)[i] == md->pe_unknown)
      continue;
    if ((*unknown_list)[i] == md->charge_balance_unknown)
    {
      fprintf(f, "\t%-15s%12.3e%12.3e", (*unknown_list)[i]->name.CharPtr(), (double) ((*unknown_list)[i]->sum / md->mass_water_aq_x), (double) (*unknown_list)[i]->sum);
      fprintf(f, "  Charge balance\n");
      pure_water = false;
      continue;
    }
    if ((*unknown_list)[i]->type == SOLUTION_PHASE_BOUNDARY)
    {
      fprintf(f, "\t%-15s%12.3e%12.3e", (*unknown_list)[i]->name.CharPtr(),
		  (double) ((*unknown_list)[i]->sum / md->mass_water_aq_x), (double) (*unknown_list)[i]->sum);
      fprintf(f, "  Equilibrium with %s\n",
		  (*unknown_list)[i]->p->name);
      pure_water = FALSE;
      continue;
    }
    if ((*unknown_list)[i]->type == MB)
    {
      fprintf(f, "\t%-15s%12.3e%12.3e\n", (*unknown_list)[i]->name.CharPtr(),
		  (double) ((*unknown_list)[i]->sum / md->mass_water_aq_x), (double) (*unknown_list)[i]->sum);
      pure_water = FALSE;
    }
  }

  if (pure_water == true)
  {
    fprintf(f, "\t%-15s\n", "Pure water");
  }
/*
 *   Description of solution
 */
  fprintf(f, "\n");
  PrintCentered(f, "Description of solution");
/*
 *   pH
 */
  fprintf(f, "%45s%7.3f    ", "pH  = ", (double) (-(gd->s_hplus->la)));
  if (md->ph_unknown == NULL)
  {
    fprintf(f, "\n");
  }
  else if (md->ph_unknown == md->charge_balance_unknown)
  {
    fprintf(f, "  Charge balance\n");
  }
  else if (md->ph_unknown->type == SOLUTION_PHASE_BOUNDARY)
  {
    fprintf(f, "  Equilibrium with %s\n",
		md->ph_unknown->p->name);
  }
  else if (md->ph_unknown->type == ALK)
  {
    fprintf(f, "  Adjust alkalinity\n");
  }

  fprintf(f, "%45s%7.3f    ", "pe  = ", (double) (-(gd->s_eminus->la)));
  if (md->pe_unknown == NULL)
  {
    fprintf(f, "\n");
  }
  else if (md->pe_unknown == md->charge_balance_unknown)
  {
    fprintf(f, "  Charge balance\n");
  }
  else if (md->pe_unknown->type == SOLUTION_PHASE_BOUNDARY)
  {
    fprintf(f, "  Equilibrium with %s\n",
		md->pe_unknown->p->name);
  }
  else if (md->pe_unknown->type == MH)
  {
    fprintf(f, "  Adjusted to redox equilibrium\n");
  }

  EC = CalcSC ();
  if (EC > 0) {
  fprintf(f, "%45s%i\n", "Specific Conductance (uS/cm, 25 oC) = ",
	      (int) EC);
  }
  fprintf(f, "%45s%7.3f\n", "Activity of water  = ",
	      exp (gd->s_h2o->la * md->LOG_10));
  fprintf(f, "%45s%11.3e\n", "Ionic strength  = ",
	      (double) md->mu_x);
  fprintf(f, "%45s%11.3e\n", "Mass of water (kg)  = ",
	      (double) md->mass_water_aq_x);
  if (md->alkalinity_unknown == NULL)
  {
    fprintf(f, "%45s%11.3e\n",
		"Total alkalinity (eq/kg)  = ",
		(double) (md->total_alkalinity / md->mass_water_aq_x));
  }
  if (md->carbon_unknown == NULL)
  {
    fprintf(f, "%45s%11.3e\n", "Total carbon (mol/kg)  = ",
		(double) (md->total_carbon / md->mass_water_aq_x));
  }
  fprintf(f, "%45s%11.3e\n", "Total CO2 (mol/kg)  = ",
	      (double) (md->total_co2 / md->mass_water_aq_x));
  fprintf(f, "%45s%7.3f\n", "Temperature (deg C)  = ",
	      (double) md->tc_x);
  fprintf(f, "%45s%11.3e\n", "Electrical balance (eq)  = ",
	      (double) md->cb_x);
  fprintf(f, "%45s%6.2f\n",
	      "Percent error, 100*(Cat-|An|)/(Cat+|An|)  = ",
	      (double) (100 * md->cb_x / md->total_ions_x));
  fprintf(f, "%45s%3d\n", "Iterations  = ", md->iterations);
  fprintf(f, "%45s%e\n", "Total H  = ", (double) md->total_h_x);
  fprintf(f, "%45s%e\n", "Total O  = ", (double) md->total_o_x);
  fprintf(f, "\n");
}

LDBLE ModelEngine::CalcSC()
{
  int i;
  LDBLE lm, a, z, Dw, SC, ff;
	SpeciesInfo *s, *s_b;

  SC = 0;
  for (i = 0; i < md->species_info_list->Count(); i++)
  {
		s = (*md->species_info_list)[i];

		if (s->s->type == EX)
      continue;
    if (s->s->type == SURF)
      continue;
    if (s->s == gd->s_h2o)
      continue;
    if ((Dw = s->s->dw) == 0)
      continue;
    if ((z = fabs(s->s->z)) == 0)
      continue;
    if (i > 0)
		{
			s_b = (*md->species_info_list)[i - 1];
			
			if(s->s->name == s_b->s->name)
				continue;
		}

    lm = s->s->lm;
    if (lm > -9)
    {
      ff = (md->mu_x < (LDBLE).36 * z ? (LDBLE)0.6 / sqrt(z) : sqrt(md->mu_x) / z);
	  
      a = Under(lm + ff * s->s->lg);
      SC += a * z * z * Dw;
    }
  }
  SC *= (LDBLE)1e7 * F_C_MOL * F_C_MOL / (R_KJ_DEG_MOL * (LDBLE)298160.0);
  return SC;
}

void ModelEngine::PrintCentered(FILE *f, String str)
{
  int i, l, l1, l2;
  String token;

  l = (int) str.Length();
  l1 = (79 - l) / 2;
  l2 = 79 - l - l1;
  for (i = 0; i < l1; i++)
    token += '-';
  token += str;
  for (i = 0; i < l2; i++)
    token += '-';
  fprintf (f, "%s\n\n", token);
}

void ModelEngine::PrintSpecies(FILE *f)
{
  int i;
  String name, name1;
  Master *master_ptr;
  LDBLE min;
  LDBLE lm;

  min = -1000;
  PrintCentered(f, "Distribution of species");

  fprintf(f, "   %-20s%12s%12s%10s%10s%10s\n", " ", " ", " ", "Log   ", "Log   ", "Log ");
  fprintf(f, "   %-20s%12s%12s%10s%10s%10s\n\n", "Species", "Molality", "Activity", "Molality", "Activity", "Gamma");

  gd->s_h2o->lm = gd->s_h2o->la;
  name = gd->s_hplus->secondary->e->name;

	SpeciesInfo *s;
  for (i = 0; i < md->species_info_list->Count(); i++)
  {
		s = (*md->species_info_list)[i];

    if (s->s->type == EX)
      continue;
    if (s->s->type == SURF)
      continue;
    if (s->master_s->secondary != NULL)
    {
      master_ptr = s->master_s->secondary;
      name1 = s->master_s->secondary->e->name;
    }
    else
    {
      master_ptr = s->master_s->primary;
      name1 = s->master_s->primary->e->name;
    }

    if (name1 != name)
    {
      name = name1;
			fprintf(f, "%-14s%12.3e\n", name.CharPtr(), (double) (master_ptr->total / md->mass_water_aq_x));
      min = md->censor * master_ptr->total / md->mass_water_aq_x;
      if (min > 0)
      {
				min = log10 (min);
      }
      else
      {
				min = -1000.;
      }
    }

    if (s->s->lm > min)
    {
      if (s->s == gd->s_h2o)
      {
				lm = log10 (gd->s_h2o->moles / md->mass_water_aq_x);
      }
      else
      {
				lm = s->s->lm;
      }

			fprintf(f, "   %-20s%12.3e%12.3e%10.3f%10.3f%10.3f\n", s->s->name.CharPtr(),
						  (double) ((s->s->moles) / md->mass_water_aq_x), 
							(double) Under(s->s->lm + s->s->lg), 
							(double) lm, 
							(double) (s->s->lm + s->s->lg), 
							(double) s->s->lg);
    }
  }
  fprintf(f, "\n");
}

void ModelEngine::PrintResults(String filename)
{
	FILE *f = fopen(filename.CharPtr(), "wt");

	PrintTotals(f);
	PrintSpecies(f);

	fclose(f);
}


bool ModelEngine::Prepare()
{
  if (md->use.sol_p == NULL)
		return false;

	Solution	*sol_p = md->use.sol_p;

  if (!Clear())
		return false;

	SetupUnknowns();

  if (md->state == INITIAL_SOLUTION)
		if (!ConvertUnits(sol_p))
			return false;

	if (!SetupSolution())
		return false;

	if (!SetupExchange())
		return false;

	if (!SetupSurface())
		return false;

	if (!SetupPurePhases())
		return false;

	if (!SetupGasPhases())
		return false;

	if (!SetupSSAssemblage())
		return false;

	if (!SetupRelatedSurface())
		return false;

	if (!TidyRedox())
		return false;

	ArrNewMax((md->max_unknowns + 1) * md->max_unknowns);
	ResidualNewMax(md->max_unknowns);
	DeltaNewMax(md->max_unknowns);

	/*
	arr_max = (md->max_unknowns + 1) * md->max_unknowns;
	if (arr_max > arr_capacity)
	{
		delete [] arr;
		arr = NULL;
		arr_capacity = arr_max;
		arr = new LDBLE [arr_capacity];
	}

	residual_max = md->max_unknowns;
	if (residual_max > residual_capacity)
	{
		delete [] residual;
		residual = NULL;
		residual_capacity = residual_max;
		residual = new LDBLE [residual_capacity];
	}

	delta_max = md->max_unknowns;
	if (delta_max > delta_capacity)
	{
		delete [] delta;
		delta = NULL;	
		delta_capacity = delta_max;
		delta = new LDBLE [delta_capacity];
	}
	*/

  if (!BuildModel())
		return false;
  
  return true;
}

bool ModelEngine::BuildPurePhases()
{
	//
	// Includes calculation of inverse saturation index in sum_mb.
	// Puts coefficients in iap and mass balance equations for each phase.
	//

	if (md->pure_phase_unknown == NULL) //ToDo: CHECK THIS!!!!
    return true;
	
	int i;
  bool stop;
	int j, k, l;
	char token[DEFAULT_STRING_LENGTH];
  char *ptr;
  Master *master_ptr, *m2, *m2_0;
  ReactionToken *rxn_ptr;
	Unknown *x2;

	// Build into sums the logic to calculate inverse saturation indices for pure phases

	// Calculate inverse saturation index
	Unknown *u;
  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

		if (u->type != PP || u->p->rxn_x->token_list->Count() <= 0)      
			continue;

		// the code below is "strange", because if md->pure_phase_unknown was NULL, this function do not "execute"
    if (md->pure_phase_unknown == NULL)
      md->pure_phase_unknown = u;

    StoreMB(&u->p->lk, &u->f, 1.0);
    StoreMB(&u->si, &u->f, 1.0);

		for (j = 1; j < u->p->rxn_x->token_list->Count(); j++)
    {
			rxn_ptr = u->p->rxn_x->token_list->Element(j);
			StoreMB(&rxn_ptr->s->la, &u->f, -rxn_ptr->coef);
    }
  }

  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

		// rxn_x is null if an element in phase is not in solution
    if (u->type != PP || u->p->rxn_x == NULL)
      continue;

		// Put coefficients into IAP equations
		for (j = 1; j < u->p->rxn_x->token_list->Count(); j++)
    {
			rxn_ptr = u->p->rxn_x->token_list->Element(j);

			if (rxn_ptr->s->secondary != NULL && rxn_ptr->s->secondary->in)
				master_ptr = rxn_ptr->s->secondary;
      else
				master_ptr = rxn_ptr->s->primary;

			if (master_ptr == NULL || master_ptr->u == NULL)
				continue;

      StoreJacob0(u->number, master_ptr->u->number, rxn_ptr->coef);
    }

		//Put coefficients into mass balance equations
		eos_list->Clear();
    parent_count = 0;

    if (!u->pure_phase->add_formula.IsEmpty())
      u->pure_phase->add_formula.Copy(token);
    else
      u->p->formula.Copy(token);

    ptr = &token[0];
		GetElementsInSpecies(&ptr, 1.0);

		// Go through elements in phase
		ChangeHydrogenInEOSList(0);

		ElementOfSpecies *eos;
		for (j = 0; j < eos_list->Count(); j++)
    {
			eos = (*eos_list)[j];

      if (eos->e->name == "H" && md->mass_hydrogen_unknown != NULL)
      {
				StoreJacob0 (md->mass_hydrogen_unknown->number, u->number, -eos->coef);
				StoreSumDeltas (&delta[i], &md->mass_hydrogen_unknown->delta, eos->coef);

      }
      else if (eos->e->name == "O" && md->mass_oxygen_unknown != NULL)
      {
				StoreJacob0 (md->mass_oxygen_unknown->number, u->number, -eos->coef);
				StoreSumDeltas (&delta[i], &md->mass_oxygen_unknown->delta, eos->coef);
      }
      else
      {
				master_ptr = eos->e->primary;

				if (!master_ptr->in)
					master_ptr = master_ptr->s->secondary;

				if (master_ptr == NULL || !master_ptr->in)
					return false;
				else if (master_ptr->in)
				{
					StoreJacob0 (master_ptr->u->number, u->number, -eos->coef);
				  StoreSumDeltas (&delta[i], &master_ptr->u->delta,	eos->coef);
				}
				else if (master_ptr->rewrite)
				{
					stop = false;

					for (k = 0; k < count_unknowns; k++)
					{
						x2 = (*unknown_list)[k];

						if (x2->type != MB)
							continue;

						for (l = 0; l < x2->master->Count(); l++)
						{
							m2 = (*x2->master)[l];							

							if (m2 == master_ptr)
							{
								m2_0 = (*x2->master)[0];
								StoreJacob0 (m2_0->u->number, u->number, -eos->coef);
								StoreSumDeltas (&delta[i], &m2_0->u->delta, eos->coef);

								stop = true;
								break;
							}
						}

						if (stop)
							break;
					}
				}
      }
    }
  }

  return true;
}

bool ModelEngine::SetInitialMoles(int i)
{
  int j, k;
	PurePhase *pp;
	GasComp *gc_p;
	SS *ss_p;
	SSComp *ssc_p;
	PPAssemblage *ppa;
	GasPhase *gas;
	SSAssemblage  *ssa;

	// Pure phase assemblage
	ppa = md->use.ppa_list->Element(i);
	if (ppa != NULL)
	{
    for (j = 0; j < ppa->pure_phases->Count(); j++)
    {
      pp = (*ppa->pure_phases)[j];
			pp->initial_moles = pp->moles;

      if (pp->initial_moles < 0)
				pp->initial_moles = 0;
    }
	}

	// Gas phase  
	gas = md->use.gas_list->Element(i);
  if (gas != NULL)
  {
		for (j = 0; j < gas->comps->Count(); j++)
		{
			gc_p = (*gas->comps)[j];
      gc_p->initial_moles = gc_p->moles;
		}
  }

	// Solid solutions	
	ssa = md->use.ssa_list->Element(i);
	if (ssa != NULL)
  {
		for (k = 0; k < ssa->ss_list->Count(); k++)
    {
			ss_p = (*ssa->ss_list)[k];

			for (j = 0; j < ss_p->comps_list->Count(); j++)
			{
				ssc_p = (*ss_p->comps_list)[j];
				ssc_p->init_moles = ssc_p->moles;
			}
    }
  }

  return true;
}

bool ModelEngine::RunReactions(int i, LDBLE step_fraction)
{
	CONVERGE_RESULT converge;

  md->run_reactions_iterations = 0;
  
  converge = SetAndRunWrapper (i, step_fraction, i);

  if (converge == MASS_BALANCE_CR)
		throw EModelEngine(1, RUN_REACTIONS_MEF);

  md->run_reactions_iterations += md->iterations;
	md->iterations = md->run_reactions_iterations;

  return true;
}

CONVERGE_RESULT ModelEngine::SetAndRunWrapper(int i, LDBLE step_fraction, int nsaver)
{
  int j, max_try;
  int old_itmax;
	bool old_diag;
  LDBLE old_tol, old_min_value, old_step, old_pe, old_pp_column_scale;
  LDBLE small_pe_step, small_step;

	bool ppa_saved, ssa_saved;

	CONVERGE_RESULT converge;

	/*
  struct pp_assemblage *pp_assemblage_save = NULL;
  struct s_s_assemblage *s_s_assemblage_save = NULL;
  struct kinetics *kinetics_save = NULL;
	*/

  /* 0 -- normal */
  /* 1 -- try smaller step size, more iterations */
  /* 2 -- try diagonal scaling */
  /* 3 -- try smaller tolerance */
  /* 4 -- try alternate scaling */
  small_pe_step = 5.;
  small_step = 10.;
	converge = NOT_CONVERGED_CR;

  old_diag = md->diagonal_scale;
  old_itmax = md->itmax;
  old_tol = md->ineq_tol;
  old_step = md->step_size;
  old_pe = md->pe_step_size;
  old_min_value = md->min_value;
  old_pp_column_scale = md->pp_column_scale;

  if (md->state == REACTION)
  {
    if (!SetReaction(i)) 
			throw exception(); 
  }

  if (md->use.ppa_p != NULL)
	{
		md->use.ppa_p->CopyTo(&ppa_save_1);
		ppa_saved = true;
	}
	else
		ppa_saved = false;

  if (md->use.ssa_p != NULL)
  {
		md->use.ssa_p->CopyTo(&ssa_save_1);
		ssa_saved = true;
  }
	else
		ssa_saved = false;

	max_try = 11;

  for (j = 0; j < max_try; j++)
  {
    if (j == 1)
    {
      if (md->pe_step_size <= small_pe_step && md->step_size <= small_step)
				continue;

      md->itmax *= 2;
      md->step_size = small_step;
      md->pe_step_size = small_pe_step;
    }
    else if (j == 2)
    {
      md->itmax *= 2;
			md->diagonal_scale = !md->diagonal_scale;
    }
    else if (j == 3)
    {
      md->itmax *= 2;
      md->ineq_tol /= 10.;
    }
    else if (j == 4)
    {
      md->itmax *= 2;
      md->ineq_tol *= 10.;
    }
    else if (j == 5)
    {
      md->itmax *= 2;
			md->diagonal_scale = !md->diagonal_scale;
      md->ineq_tol /= 10.;
    }
    else if (j == 6)
    {
      md->itmax *= 2;
      md->pp_column_scale = (LDBLE)1e-10;
    }
    else if (j == 7)
    {
      md->itmax *= 2;
      md->pp_column_scale = (LDBLE)1e-10;
			md->diagonal_scale = !md->diagonal_scale;
    }
    else if (j == 8)
    {
      md->itmax *= 2;
      md->min_value *= 10;
    }
    else if (j == 9)
    {
      md->aqueous_only = 5;
    }
    else if (j == 10)
    {
      md->negative_concentrations = true;
    }

    if (j > 0)
    {
			if (ppa_saved)
				ppa_save_1.CopyTo(md->use.ppa_p);

			if (ssa_saved)
				ssa_save_1.CopyTo(md->use.ssa_p);
    }

		converge = SetAndRun(i, nsaver, step_fraction);

    md->diagonal_scale = old_diag;
    md->itmax = old_itmax;
    md->ineq_tol = old_tol;
    md->step_size = old_step;
    md->pe_step_size = old_pe;
    md->min_value = old_min_value;
    md->pp_column_scale = old_pp_column_scale;
    md->aqueous_only = 0;
    md->negative_concentrations = false;
		
    if (converge == OK_CR || converge == CONVERGED_CR || converge == MASS_BALANCE_CR)
      break;
  }

  if (converge == NOT_CONVERGED_CR)
		throw EModelEngine(2, SET_AND_RUN_WRAPPER_MEF);

	return converge;
}

CONVERGE_RESULT ModelEngine::SetAndRun(int i, int nsaver, LDBLE step_fraction)
{
	//
	//  i			--user number for soln, reaction, etc.
	//  use_mix	--integer flag
	//	state == TRANSPORT: DISP, STAG, NOMIX
	//	state == REACTION: TRUE, FALSE
	//  use_kinetics --true or false flag to calculate kinetic reactions
	//  nsaver	   --user number to store solution
	//  step_fraction--fraction of irreversible reaction to add
	//

  //int n, n1, n2;
	
	CONVERGE_RESULT converge;

  if (md->state == REACTION)
  {
    if (!SetReaction (i))
			throw exception();
  }

	// Take step
	if (md->state >= REACTION)
  {
    if (Step(step_fraction) == MASS_BALANCE_CR)
      return MASS_BALANCE_CR;

		// Always use solution, exchange, and surface -1
		md->use.sol_p = md->use.sol_list->Element(1); //md->use.SearchSolution(-1);

    if (md->use.exc_p != NULL)
			md->use.exc_p = md->use.exc_list->Element(1); //md->use.SearchExchange(-1);

		if (md->use.sur_p != NULL)
      md->use.sur_p = md->use.sur_list->Element(1); //md->use.SearchSurface(-1);
  }

  if (md->use.sur_p != NULL)
    md->dl_type_x = md->use.sur_p->dl_type;

  if (md->use.sur_p != NULL && md->dl_type_x != NO_DL)
    converge = ExecuteSurfaceModel();
  else
	{

#ifdef DEBUG_MOHID
		//d.debug_status = true;
#endif

		if (!Prepare())
			throw exception();

		if (!KTemp(md->use.sol_p->tc))
			throw exception();

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.set_f, "called_by_set_and_run: %d\n", d.set_count++);
#endif

		if (!Set(false))
			throw exception();

		converge = ExecuteModel();

#ifdef DEBUG_MOHID
		//d.debug_status = true;
#endif
	}

	SumSpecies();

	return (converge);
}

CONVERGE_RESULT ModelEngine::Step(LDBLE step_fraction)
{
	//
  //   zero global solution, add solution or mixture, add exchange,
  //   add surface, add gas phase, add solid solutions,
  //   set temperature, and add reaction. 
  //   Ensure all elements
  //   included in any of these are present in small amounts.
  //   Save result as n_user -1.
	//

  //LDBLE difftemp;
  int step_number;

	// Zero out global solution data
	md->ResetXVariables(); //xsolution_zero ();

	// Set reaction to zero
  md->step_x = 0.0;
  step_number = md->reaction_step;	

  if (!AddSolution(md->use.sol_p, 1.0, 1.0))
		throw exception();

	// Exchange
	if (md->use.exc_p != NULL && !AddExchange(md->use.exc_p))
		throw exception();

	// Surface
  if (md->use.sur_p != NULL && !AddSurface(md->use.sur_p))
		throw exception();

	// Gases
	if (md->use.gas_p != NULL && !AddGasPhase(md->use.gas_p))
		throw exception();

	// Pure phases and solid solutions are added to avoid zero or negative concentrations
	// Pure phases
	if (md->use.ppa_p != NULL)
	{
		md->use.ppa_p->CopyTo(&ppa_save_2);

		if (!AddPPAssemblage(md->use.ppa_p))
			throw exception();
	}

 // Solid solutions
  if (md->use.ssa_p != NULL)
  {
		md->use.ssa_p->CopyTo(&ssa_save_2);

    if (!AddSSAssemblage(md->use.ssa_p))
			throw exception();
  }

 // Check that elements are available for gas components, pure phases, and solid solutions
  if (md->use.gas_p != NULL)
    GasPhaseCheck(md->use.gas_p);

	if (md->use.ppa_p != NULL)
		PPAssemblageCheck(md->use.ppa_p);

	if (md->use.ssa_p != NULL)
		SSAssemblageCheck(md->use.ssa_p);

	// Check that element moles are >= zero
  if (SolutionCheck() == MASS_BALANCE_CR)
  {
		// reset moles and deltas
		if (md->use.ppa_p != NULL)
			ppa_save_2.CopyTo(md->use.ppa_p);
		
    if (md->use.ssa_p != NULL)
			ssa_save_2.CopyTo(md->use.ssa_p);

    return (MASS_BALANCE_CR);
  }

	// Copy global into solution n_user = -1
	SaveSolution(md->use.sol_list->Element(1)); //ToDo: CHECK how will be done this
  StepSaveSurf(md->use.sur_list->Element(1)); //ToDo: Criar estas funções
  StepSaveExch(md->use.exc_list->Element(1)); //ToDo: Criar estas funções

  return (OK_CR);
}


void ModelEngine::SaveSolution(Solution *dest)
{
  int i; //, j, n;
  //int count_mass_balance, count_master_activity;
  int max_mass_balance, max_master_activity;
  
  Master *m; //, *mi;
	MasterActivity *ma;
	Conc *c;

  max_mass_balance = MAX_MASS_BALANCE;
  max_master_activity = MAX_MASS_BALANCE;

  dest->tc = md->tc_x;
  dest->ph = md->ph_x;
  dest->solution_pe = md->solution_pe_x;
  dest->mu = md->mu_x;
	//d.PrintSpeciesInfoToFile("SaveSolution-1", md->species_info_list, md, gd);
  dest->ah2o = md->ah2o_x;
  dest->density = md->density_x;
  dest->total_h = md->total_h_x;
  dest->total_o = md->total_o_x;
  dest->cb = md->cb_x;	
  dest->mass_water = md->mass_water_aq_x;
  dest->total_alk = md->total_alkalinity;
	dest->units = _mol_kgw_;

	dest->InitPE();
  dest->default_pe = 0;

	dest->totals->Clear();
	dest->ma->Clear();

  for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];

    if (m->s->type == EX || m->s->type == SURF || m->s->type == SURF_PSI)
      continue;
    if (m->s == gd->s_hplus || m->s == gd->s_h2o) 
			continue;

		// Save list of log activities
		if (m->in || m->rewrite)
    {
			ma = dest->ma->AddNew();

			ma->description = m->name;
      ma->la = m->s->la;
    }

    if (m->total <= MIN_TOTAL)
    {
      m->total = 0.0;
      m->total_primary = 0.0;
      continue;
    }

		// Save list of concentrations
		c = dest->totals->AddNew();

    c->name = m->name;
    c->input_conc = m->total;
    c->moles = m->total;
    c->units = dest->units;
    c->equation_name = "";
    c->n_pe = 0;
    c->p = NULL;
    c->phase_si = 0.0;
    c->as = "";
    c->gfw = 0.0;
  }

	dest->species_gamma->Clear();
}

bool ModelEngine::CheckPPAssemblage(PPAssemblage *ppa)
{
	//
	// Check list of all elements in pure_phase assemblage to see
	// if all are in model. Return true if all are present,
	// Return false if one or more is missing.
	//

 int j;

  Master *m;
  for (j = 0; j < ppa->eos_list->Count(); j++)
  {
		m = (*ppa->eos_list)[j]->e->primary;

    if (m->s == gd->s_h2o || m->s == gd->s_hplus)
      continue;

    if (m->total > MIN_TOTAL)
      continue;

    return false;
  }

  return true;
}

CONVERGE_RESULT ModelEngine::SolutionCheck()
{
	//
	// Check for missing elements
	//

	int i;
  Master *m;

	// Check that all elements are in solution for phases with zero mass
  for (i = 0; i < gd->master_list.Count(); i++)
  {
    m = gd->master_list[i];

    if (m->total >= 0.0)
      continue;

    if (m->total > -MIN_TOTAL)
    {
      m->total = 0;
      continue;
    }

    if (m->s == gd->s_eminus || m->s == gd->s_h2o	|| m->s == gd->s_hplus || m->s == gd->s_h3oplus)
    {
      m->total = 0;
      continue;
    }

		return (MASS_BALANCE_CR);
  }

  return (OK_CR);
}




//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupSolution()
{
	// 
	// Fills in data in unknown structure for the solution
	//

  int x_index, j;
  Master *m;
	Conc *conc;
	String token, token2, token_list;
	Unknown *x;

	Solution *sol_p = md->use.sol_p;
	count_unknowns = 0;

  for (x_index = 0; x_index < sol_p->totals->Count(); x_index++)
  {
		x = (*unknown_list)[count_unknowns];

		conc = (*sol_p->totals)[x_index];
		token = conc->name(0, " ");

    if (conc->input_conc <= 0.0 && token != "H(1)" && token != "E")
			continue;

		m = gd->master_list.Search(&token, true);
    if (m == NULL)
			return false;

    if (m->type != AQ)
			return false;

		conc->name.Trim(_all_);
		int p = conc->name.Find(" ");
		if (p != -1)
			token_list = conc->name(p + 1);
		else
			token_list = "";

		// Store list of master species pointers, set master[i].in and master[i].rxn for list
		GetMasterList(token_list, m, x->master);
		SetupMasterReaction(x->master, &(*md->pe_x)[conc->n_pe]->rxn);

		// Set default unknown data
    x->type = MB;
    x->name = conc->name;
    x->total = conc;

    for (j = 0; j < x->master->Count(); j++)
		{
			m = (*x->master)[j];
      m->u = x;
		}

		x->moles = conc->moles;

#ifdef DEBUG_MOHID
		if (d.debug_status)
			fprintf(d.setup_solution_f, "SetupSolution-1) %d %20.20e\n", x_index, x->moles);

#endif

		token.ToLower();

		// Set pointers
		if (token.Find("alk", true) != -1)
    {
      if (md->alkalinity_unknown == NULL)
      {
				x->type = ALK;
				md->alkalinity_unknown = x;
      }
      else
				return false;
    }
    else if (token == "c" || token == "c(4)")
    {
      if (md->carbon_unknown == NULL)
				md->carbon_unknown = x;
      else
				return false;
    }
    else if (token == "h(1)")
    {
      if (md->ph_unknown == NULL)
				md->ph_unknown = x;
      else
				return false;
    }
    else if (token == "e")
    {
      if (md->pe_unknown == NULL)
				md->pe_unknown = x;
      else
				return false;
    }

		// Charge balance unknown
		if (!conc->equation_name.IsEmpty())
    {
			if (conc->charge)
      {
				if (md->charge_balance_unknown == NULL)
				{
					md->charge_balance_unknown = x;
					x->type = CB;

					if (md->charge_balance_unknown == md->ph_unknown)
						x->moles = sol_p->cb;
				}
				else
					return false;
      }
      else
      {
				conc->p = gd->phase_list.Search(&conc->equation_name, true);

				if (conc->p == NULL)
					return false;

				x->type = SOLUTION_PHASE_BOUNDARY;
				x->p = conc->p;
				x->si = conc->phase_si;

				if (md->solution_phase_boundary_unknown == NULL)
					md->solution_phase_boundary_unknown = x;
      }
    }

		//conc->x = x;
		count_unknowns++;
  }

	// Set mb_unknown
  if (count_unknowns > 0)
    md->mb_unknown = (*unknown_list)[0];

	// Special for alkalinity
  if (md->alkalinity_unknown != NULL)
  {
    if (md->carbon_unknown != NULL)
    {
			// pH adjusted to obtain given alkalinity
      if (md->ph_unknown == NULL)
      {
				md->ph_unknown = md->alkalinity_unknown;
				m = gd->master_list.Search("H(1)", true);

				if (!md->alkalinity_unknown->master->Store(0, m))
					return false;

				m->in = true;
				m->rewrite = false;
				m->u = md->ph_unknown;

				if (!md->ph_unknown->master->Store(0, m))
					return false;

				md->ph_unknown->name = "H(1)";
      }
      else
				return false;
    }
    else
    {
			m = (*md->alkalinity_unknown->master)[0]->s->secondary;

      if (m != NULL)
      {
				m->in = true;
				m->rewrite = false;
				m->u =	md->alkalinity_unknown;
      }
      else
				return false;
    }
  }

	// Ionic strength
  md->mu_unknown = (*unknown_list)[count_unknowns];
	md->mu_unknown->number = count_unknowns++;
  md->mu_unknown->name = "Mu";
  md->mu_unknown->type = MU;
  md->mu_unknown->moles = 0.0;

	// Activity of water
  md->ah2o_unknown = (*unknown_list)[count_unknowns];
	md->ah2o_unknown->number = count_unknowns++;
  md->ah2o_unknown->name = "A(H2O)";
  md->ah2o_unknown->type = AH2O;
	md->ah2o_unknown->master->Clear();
  md->ah2o_unknown->master->AddNew(gd->master_list.Search("O", true));
  (*md->ah2o_unknown->master)[0]->u = md->ah2o_unknown;
  md->ah2o_unknown->moles = 0.0;

  if (md->state >= REACTION)
  {
		// Reaction: pH for charge balance
    md->ph_unknown = (*unknown_list)[count_unknowns];
		md->ph_unknown->number = count_unknowns++;
    md->ph_unknown->name = "pH";
    md->ph_unknown->type = CB;
    md->ph_unknown->moles = sol_p->cb;
		md->ph_unknown->master->Clear();
    md->ph_unknown->master->AddNew(gd->s_hplus->primary);
    (*md->ph_unknown->master)[0]->u = md->ph_unknown;
    md->charge_balance_unknown = md->ph_unknown;
    
		//Reaction: pe for total hydrogen
		md->pe_unknown = (*unknown_list)[count_unknowns];
		md->pe_unknown->number = count_unknowns++;
    md->mass_hydrogen_unknown = md->pe_unknown;
    md->mass_hydrogen_unknown->name = "Hydrogen";
    md->mass_hydrogen_unknown->type = MH;
		md->mass_hydrogen_unknown->moles = sol_p->total_h - 2 * sol_p->total_o;
		md->mass_hydrogen_unknown->master->Clear();
    md->mass_hydrogen_unknown->master->AddNew(gd->s_eminus->primary);
    (*md->mass_hydrogen_unknown->master)[0]->u = md->mass_hydrogen_unknown;  

		// Reaction H2O for total oxygen
		md->mass_oxygen_unknown = (*unknown_list)[count_unknowns];
		md->mass_oxygen_unknown->number = count_unknowns++;
		md->mass_oxygen_unknown->name = "Oxygen";
    md->mass_oxygen_unknown->type = MH2O;
    md->mass_oxygen_unknown->moles = sol_p->total_o;
		md->mass_oxygen_unknown->master->Clear();
    md->mass_oxygen_unknown->master->AddNew(gd->s_h2o->primary);
  }

	// Validity tests
	if (md->ph_unknown != NULL && md->ph_unknown == md->charge_balance_unknown && md->alkalinity_unknown != NULL)
		return false;

  if (md->alkalinity_unknown != NULL && (md->alkalinity_unknown->type == CB ||  md->alkalinity_unknown->type == SOLUTION_PHASE_BOUNDARY))
		return false;

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupExchange()
{
	//
	// Fill in data for exchanger in unknowns structures
	//

  int i, j;
  Master *m_p;
	Use *use;
	ExchComp *excc_p;
	ElementOfSpecies *eos_p;
	Unknown *u;

	use = &md->use;
	if (use->exc_p == NULL)
    return true;

	for (j = 0; j < use->exc_p->comps->Count(); j++)
  {
		excc_p = (*use->exc_p->comps)[j];

		for (i = 0; i < excc_p->totals->Count(); i++)
    {
			eos_p = (*excc_p->totals)[i];

			// Find master species
      if ((m_p = eos_p->e->master) == false)
				return false;

      if (m_p->type != EX)
				continue;

			// Check for data already given
			if (m_p->in || m_p->rewrite) //ToDo: Check if the "rewrite" check is necessary or if is wrong
				(*unknown_list)[m_p->u->number]->moles += eos_p->coef;
      else
      {
				// Set flags
				m_p->in = true;
				m_p->rewrite = false;

				// Set unknown data
				u = (*unknown_list)[count_unknowns++];

				u->type = EXCH;
				u->exch_comp = (*use->exc_p->comps)[j];
				u->name = eos_p->e->name;
				u->moles = eos_p->coef;
				u->master->Clear();
				u->master->AddNew(m_p);
				(*u->master)[0]->u = u;
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupSurface()
{
	//
	// Fill in data for surface assemblage in unknown structure
	//
 
  int i, j, plane;
  int mb_unknown_number, type;

	Use *use;
	SurfaceComp *sc_p;
	ElementOfSpecies *eos_p;
	Master *m_p;
	Unknown *u, *u_prior, *u_p, **u_target, *u2;

	String token, 
				 mass_balance_name,
				 cb_suffix,
				 psi_suffix,
				 name1,
				 name2;

	use = &md->use;

	if (use->sur_p == NULL)
    return true;

	for (i = 0; i < use->sur_p->comps->Count(); i++)
  {
		sc_p = (*use->sur_p->comps)[i];

		// Find master species for each surface, setup unknown structure
    for (j = 0; j < sc_p->totals->Count(); j++)
    {
			eos_p = (*sc_p->totals)[j];

			if ((m_p = eos_p->e->master) == NULL)
				return false;

      if (m_p->type != SURF)
				continue;

  		// Check that data not already given
			if (m_p->in || m_p->rewrite)
				return false;

			// Set flags
      sc_p->master = m_p;
			sc_p->name = m_p->name;
			m_p->in = true;
			m_p->rewrite = false;

			// Setup mass balance unknown
			u = (*unknown_list)[count_unknowns];			
			u_prior = (*unknown_list)[count_unknowns - 1];
      
			u->type = SURFACE;
      u->name = eos_p->e->name;
      u->number = count_unknowns++;
      u->surface_comp = sc_p;
      u->master->Clear();
			u->master->AddNew(m_p);
      (*u->master)[0]->u = u;
      u->moles = eos_p->coef;

      if (md->surface_unknown == NULL)
				md->surface_unknown = u;

      u->potential_unknown = NULL;

      if (use->sur_p->type == DDL)
      {
				// Setup surface-potential unknown
				token = m_p->e->name;
				u_p = FindSurfaceChargeUnknown (&token, SURF_PSI);

				if (u_p != NULL)
					u_prior->potential_unknown = u_p;
				else
				{
					// Find master species
 					token.Replace("_CB", "_psi");
	  
					m_p = gd->master_list.Search(&token, true);
					m_p->in = true;
					m_p->rewrite = false;
	    
					// Find surface charge structure
					u = (*unknown_list)[count_unknowns];
					u_prior = (*unknown_list)[count_unknowns - 1];

					u->type = SURFACE_CB; 
					u->surface_charge = (*use->sur_p->charge)[sc_p->charge];
					u->related_moles = u->surface_charge->grams;
					u->mass_water = (*use->sur_p->charge)[sc_p->charge]->mass_water;

					token.Replace("_psi", "_CB");

					u->name = token;
					u->master->Clear();
					u->master->AddNew(m_p);
					(*u->master)[0]->u = u;
					u->moles = 0.0;

					u_prior->potential_unknown = u;
					u->surface_comp = u_prior->surface_comp;
					count_unknowns++;
				}
      }
			else if (use->sur_p->type == CD_MUSIC)
      {
				// Setup 3 surface-potential unknowns
 				mb_unknown_number = count_unknowns - 1;			
				mass_balance_name = token;

				for (plane = SURF_PSI; plane <= SURF_PSI2; plane++)
				{
					cb_suffix = "_CB";
					psi_suffix = "_psi";

				  u_target = NULL;
	  
					type = SURFACE_CB;

					switch (plane)
					{
					case SURF_PSI:
						type = SURFACE_CB;
						u_target = &(*unknown_list)[mb_unknown_number]->potential_unknown;
						break;
					case SURF_PSI1:
						cb_suffix += "b";
						psi_suffix += "b";
						type = SURFACE_CB1;
						u_target = &(*unknown_list)[mb_unknown_number]->potential_unknown1;
						break;
					case SURF_PSI2:
						cb_suffix += "d";
						psi_suffix += "d";
						type = SURFACE_CB2;
						u_target = &(*unknown_list)[mb_unknown_number]->potential_unknown2;
						break;
					}

					token = m_p->e->name;
					u_p = FindSurfaceChargeUnknown (token, plane);

					if (u_p != NULL)
						*u_target = u_p;
					else
					{
						// Find master species
						token.Replace(cb_suffix, psi_suffix);

						m_p = gd->master_list.Search(&token, true);
						m_p->in = true;
						m_p->rewrite = false;
	    
	      		// Find surface charge structure
						u = (*unknown_list)[count_unknowns];
						u->type = type;
						u->surface_charge = (*use->sur_p->charge)[sc_p->charge];
						u->related_moles = u->surface_charge->grams;
						u->mass_water = (*use->sur_p->charge)[sc_p->charge]->mass_water;
						
						token.Replace(psi_suffix, cb_suffix);

						u->name = token;
						u->master->Clear();
						u->master->AddNew(m_p);
						
						// Find surface charge structure
				    if (plane == SURF_PSI)
							(*unknown_list)[mb_unknown_number]->potential_unknown = u;
						else if (plane == SURF_PSI1)
							(*unknown_list)[mb_unknown_number]->potential_unknown1 = u;
						else if (plane == SURF_PSI2)
							(*unknown_list)[mb_unknown_number]->potential_unknown2 = u;

						(*u->master)[0]->u = u;
						u->moles = 0.0;
						u->surface_comp = (*unknown_list)[mb_unknown_number]->surface_comp;
						count_unknowns++;
					}
				}

				// Add SURFACE unknown to a list for SURF_PSI
				u_p = FindSurfaceChargeUnknown(token, SURF_PSI);
				u_p->comp_unknowns->AddNew((*unknown_list)[mb_unknown_number]);
      }
    }
  }

	// check related phases
	if (use->sur_p->related_phases)
  {
    for (i = 0; i < count_unknowns; i++)
    {
			u = (*unknown_list)[i];

      if (u->type != SURFACE_CB)
				continue;

      for (j = 0; j < count_unknowns; j++)
      {
				u2 = (*unknown_list)[j];

				if (u2->type != SURFACE)
					continue;

				if (u2->potential_unknown != u)
					continue;

				if (u2->surface_comp->phase_name != u->surface_comp->phase_name)
					return false;
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
Unknown *ModelEngine::FindSurfaceChargeUnknown (String *str, int plane)
{
	//
	// Makes name for the potential unknown and returns in str_ptr
	// Returns NULL if this unknown not in unknown list else
	// returns a pointer to the potential unknown
	//

  int i;

	str->Replace("_", " ");
	*str = (*str)(0, " ");

	if (plane == SURF_PSI)
    (*str) += "_CB";
  else if (plane == SURF_PSI1)
    (*str) += "_CBb";
  else if (plane == SURF_PSI2)
    (*str) += "_CBd";

  for (i = 0; i < count_unknowns; i++)
  {
    if ((*unknown_list)[i]->name == *str)
      return (*unknown_list)[i];
  }

  return NULL;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupPurePhases()
{
	//
	// Fills in data for pure_phase assemglage in unknown structure
	//

  int i;
	Use *use;
	Unknown *u;
	PurePhase *pp;

	use = &md->use;
	if (use->ppa_p == NULL)
    return true;

	for (i = 0; i < use->ppa_p->pure_phases->Count(); i++)
  {
		pp = (*use->ppa_p->pure_phases)[i];

		u = (*unknown_list)[count_unknowns++];

		u->type = PP;
		u->name = pp->name;
		u->moles = pp->moles;
    u->p = pp->phase;
    u->si = pp->si;
    u->delta = pp->delta;
    u->pure_phase = pp;
    u->dissolve_only = pp->dissolve_only;
    
		if (md->pure_phase_unknown == NULL)
      md->pure_phase_unknown = u;
  }
  
	return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupGasPhases()
{
	//
	// Fill in data for gas phase unknown (sum of partial pressures) in unknown structure
	//

	int i;
	Use *use;
	Unknown *u;
	GasComp *gc_p;

	use = &md->use;
  if (use->gas_p == NULL)
    return true;

	// One for total moles in gas
	u = (*unknown_list)[count_unknowns++];

  u->type = GAS_MOLES;
  u->name = "gas moles";
  u->moles = 0.0;

	for (i = 0; i < use->gas_p->comps->Count(); i++)
  {
		gc_p = (*use->gas_p->comps)[i];
    u->moles += gc_p->moles;
  }

  if (u->moles <= 0)
    u->moles = (LDBLE)MIN_TOTAL;

  u->ln_moles = log (u->moles);
  u->gas_phase = use->gas_p;
  md->gas_unknown = u;

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupSSAssemblage()
{
	//
	// Fill in data for solid solution unknowns (sum of partial pressures) in unknown structure
	//
  int i, j;
	Use *use;
	SS *ss_p;
	SSComp *ssc_p;
	Unknown *u;

	use = &md->use;
  if (use->ssa_p == NULL)
    return true;

	// One for each component in each solid solution
	md->s_s_unknown = NULL;

	for (j = 0; j < use->ssa_p->ss_list->Count(); j++)
  {
		ss_p = (*use->ssa_p->ss_list)[j];

		for (i = 0; i < ss_p->comps_list->Count(); i++)
    {
			ssc_p = (*ss_p->comps_list)[i];

			u = (*unknown_list)[count_unknowns];

      u->type = S_S_MOLES;
			u->name =	ssc_p->name;

      if (ssc_p->moles <= 0)
				ssc_p->moles = (LDBLE)MIN_TOTAL_SS;

      u->moles = ssc_p->moles;
      ssc_p->initial_moles = u->moles;
      u->ln_moles = log(u->moles);
			u->s_s = ss_p;
      u->s_s_comp =	ssc_p;
      u->s_s_comp_number = i;
      u->p = ssc_p->phase;
      u->number = count_unknowns++;
			u->p->dn = ssc_p->dn;
			u->p->dnb =	ssc_p->dnb;
			u->p->dnc = ssc_p->dnc;
	    u->p->log10_fraction_x = ssc_p->log10_fraction_x;
      u->p->log10_lambda = ssc_p->log10_lambda;

			if (md->s_s_unknown == NULL)
				md->s_s_unknown = u;      
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetupRelatedSurface()
{
	//
	// Fill in data for surface assemblage in unknown structure
	//

  int i, k;
	Use *use;
	Unknown *u, *u2, *u_prior;
	SurfaceComp *sc_p;
  //struct surface_comp *comp_ptr;

	use = &md->use;
	if (use->sur_p == NULL || !use->sur_p->related_phases)
    return true;

  for (i = 0; i < count_unknowns; i++)
  {
		u = (*unknown_list)[i];

		if (i > 0)
			u_prior = (*unknown_list)[i - 1];

		if (u->type == SURFACE && !u->surface_comp->phase_name.IsEmpty())
    {
      for (k = count_unknowns - 1; k >= 0; k--)
      {
				u2 = (*unknown_list)[k];

				if (u2->type != PP)
					continue;

				if (u2->p->name == u->surface_comp->phase_name)
					break;
      }

      if (k == -1)
				continue;

      sc_p = u->surface_comp;
      u->phase_unknown = u2;

      u->moles = u2->moles * sc_p->phase_proportion;
    }
		else if (u->type == SURFACE_CB && !u_prior->surface_comp->phase_name.IsEmpty())
    {
      for (k = count_unknowns - 1; k >= 0; k--)
      {
				u2 = (*unknown_list)[k];

				if (u2->type != PP)
					continue;
	
				if (u2->p->name == u->surface_comp->phase_name)
					break;
      }

      if (k == -1)
				continue;

      sc_p = u->surface_comp;
      u->phase_unknown = u2;

      u->related_moles = u2->moles * sc_p->phase_proportion;
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::TidyRedox()
{
	//
	// Write pe redox reactions (rxn in struct pe_data) in terms of master species
	// defined in analytical data
	//
 
  int i, j;
	String *tok1, *tok2;

	//Keep valences of oxygen and hydrogen in model, if not already in
	Master *m;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];

    if (m->primary == true && (m->s == gd->s_hplus || m->s == gd->s_h2o))
    {
      j = i + 1;

      while (j < gd->master_list.Count() && gd->master_list[j]->e->primary == m)
      {
				if (!gd->master_list[j]->in && gd->master_list[j]->s != m->s)
				{
					gd->master_list[j]->rewrite = true;
					*(gd->master_list[j]->pe_rxn) = *(m->pe_rxn);
				}

				j++;
      }
    }
  }

	//Writes equations for e- for each redox couple used in solution n
	PEData *pe_ptr;
	for (i = 0; i < md->pe_x->Count(); i++)
	{
		pe_ptr = (*md->pe_x)[i];

		if (pe_ptr->name.Compare("pe", true) == 0)
			gd->s_eminus->rxn->CopyTo(pe_ptr->rxn);
    else
    {
			String t = pe_ptr->name;
			t.Replace("/", " ");

			token.SetString(t);
			tok1 = token.NextToken();
      tok2 = token.NextToken();

			Master *m1 = gd->master_list.Search(tok1, true);
      Master *m2 = gd->master_list.Search(tok2, true);

      if (m1 != NULL && m2 != NULL)
      {
				r_temp->Reset();

				if (!RewriteMasterToSecondary(m1, m2))
					return false;

				r_temp->TRXNSwap("e-");
      }
      else
				return false;


      if (!CheckIn())
				return false;
      else
				r_temp->CopyTo(pe_ptr->rxn);
    }
  }

	//Rewrite equations to master species that are "in" the model
	for (i = 0; i < md->pe_x->Count(); i++)
	{
		pe_ptr = (*md->pe_x)[i];

		r_temp->Reset();
		r_temp->AddTRXN(pe_ptr->rxn, 1.0, false);

    if (!WriteMassActionEqnX())
			return false;
		else
			r_temp->CopyTo(pe_ptr->rxn);
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::BuildModel()
{
	//
	// Guts of prep. Determines species in model, rewrites equations,
	// builds lists for mass balance and jacobian sums.
	//

  int i;
  LDBLE coef_e;
	String h2o("H2O");

  if (gd->s_hplus == NULL || gd->s_eminus == NULL || gd->s_h2o == NULL)
		return false;

	// All these "lists" are "created" on initialization of ModelData instance
	// Insted of "allocate" the space, is used the functions Add and AddNew that take 
	// care of allocation if necessary
	// ToDo: verify if this is ok...
	// Originally: "Make space for lists of pointers to species in the model"
	md->s_x->Clear();
	md->sum_mb1->Clear();
	md->sum_mb2->Clear();
	md->sum_jacob0->Clear();
	md->sum_jacob1->Clear();
	md->sum_jacob2->Clear();
	md->sum_delta->Clear();
	md->species_info_list->Clear();

	// Pick species in the model, determine reaction for model, build jacobian
	ComputeGFW(h2o, md->gfw_water);
  md->gfw_water *= (LDBLE)0.001;

#ifdef DEBUG_MOHID
	if (d.debug_status)
		fprintf(d.buildmodel_f, "1) %20.20e\n", md->gfw_water);
#endif


	Species *s;
	for (i = 0; i < gd->species_list.Count(); i++)
  {
		s = gd->species_list[i];

    if (s->type > H2O && s->type != EX && s->type != SURF) 
			continue;

    s->in = false;
		//s->rewrite = false;

    r_temp->Reset();
    r_temp->AddTRXN(s->rxn_s, 1.0, false);

		// Check if species is in model
		s->in = CheckIn();

    if (s->in)
    {
			// for isotopes, activity of water is for 1H and 16O
			// ToDo: Probably this will not be used until "isotopes" are seted in program
      if (s->gflag == 9)
				md->gfw_water = (LDBLE)18.0 / (LDBLE)1000.0;

			s->lg = 0.0;

			md->s_x->AddNew(s);

			// Write mass action equation for current model
			if (!WriteMassActionEqnX()) 
				return false;

      if (s->type == SURF)
      {
				if (!AddPotentialFactor()) //add_potential_factor ();
					return false; 

				if (!AddCDMusicFactors(i))
					return false; //add_cd_music_factors (i);
      }

			r_temp->CopyTo(s->rxn_x);

			// Determine mass balance equations, build sums for mass balance, build sums for jacobian
			r_temp->Reset();
      r_temp->AddTRXN(s->rxn_s, 1.0, false);

			if (s->e_sec_list->Count() <= 0) 
				WriteMBEqnX();
      else
      {
				eos_list->Clear();
				CopyToTempEOSList(s->e_sec_list, 1.0);
      }

      if (s->type == SURF)
      {
				if (!AddPotentialFactor())
					return false;

				if (!AddCDMusicFactors(i))
					return false;

				if (!AddSurfaceChargeBalance())
					return false;

				if (!AddCDMusicChargeBalances(i))
					return false;
      }

      if (s->type < EMINUS)     
				MBForSpeciesAQ(i);
      else if (s->type == EX)	  
				MBForSpeciesEX(i);
      else if (s->type == SURF)	
				MBForSpeciesSURF(i);

			BuildMBSums ();
			BuildJacobianSums (i);

			// Build list of species for summing and printing
			if (s->e_sec_list->Count() <= 0)
				WriteMBForSpeciesList(i);
      else
      {
				eos_list->Clear();
				CopyToTempEOSList(s->e_sec_list, 1.0);
      }

      BuildSpeciesList(i);

#ifdef DEBUG_MOHID
			if (d.debug_status)
				fprintf(d.buildmodel_f, "2) number:%d %20.20e\n", s->number, s->la);
#endif
		}
  }

	// Sum diffuse layer water into hydrogen and oxygen mass balances
	Unknown *x;
	if (md->dl_type_x != NO_DL && md->state >= REACTION)
  {
		for (i = 0; i < count_unknowns; i++)
    {
			x = (*unknown_list)[i];

      if (x->type == SURFACE_CB)
      {
				if (md->mass_oxygen_unknown != NULL)
					StoreMB (&x->mass_water, &md->mass_oxygen_unknown->f, 1 / md->gfw_water);
      }
    }
  }

	// Rewrite phases to current master species
	Phase *phases;
	for (i = 0; i < gd->phase_list.Count(); i++)
  {
		phases = gd->phase_list[i];

		r_temp->Reset();
		r_temp->AddPhaseTRXN(phases->rxn_s, 1.0, false);
		r_temp->TRXNReverseK();

    phases->in = CheckIn();

    if (phases->in)
    {
			//Replace e- in original equation with default redox reaction
			if (!GetReactionCoef(r_temp, "e-", coef_e)) //trxn_find_coef ("e-", 1);
				return false; 

      if (!Equal(coef_e, 0.0, TOLERANCE))
				r_temp->AddTRXN((*md->pe_x)[md->default_pe_x]->rxn, coef_e, true);

			if (!WriteMassActionEqnX())
				return false;

			r_temp->TRXNReverseK();
			r_temp->CopyTo(phases->rxn_x);
			
      WritePhaseSysTotal(i);
    }
  }

  BuildSolutionPhaseBoundaries ();
  BuildPurePhases ();
  BuildMinExch ();
  BuildMinSurface ();
  BuildGasPhase ();
  BuildSSAssemblage ();

	md->species_info_list->SortByMasterOnly(gd->s_hplus);

#ifdef DEBUG_MOHID
	if (d.debug_status)
	{
		Species *d_s;
		for(int d_i = 0; d_i < md->s_x->Count(); d_i++)
		{
			d_s = (*md->s_x)[d_i];
			fprintf(d.buildmodel_f, "3) %d %20.20e\n", d_s->number, d_s->la);
		}
	}
#endif


	SaveModel (); //ToDo: it's doing nothing for now

	return true;
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::SSRoot(LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq)
{
  int i;
  LDBLE x0, y0, x1, y1, xb, miny;

	// Bracket answer
	x0 = 0.0;
  x1 = 0.0;
  y0 = SSF (x0, a0, a1, kc, kb, xcaq, xbaq);
  miny = fabs (y0);

  for (i = 1; i <= 10; i++)
  {
    x1 = (LDBLE) i / 10;
    y1 = SSF (x1, a0, a1, kc, kb, xcaq, xbaq);

    if (fabs (y1) < miny)
      miny = fabs (y1);

		if (y0 * y1 < 0)
      break;
    else
    {
      x0 = x1;
      y0 = y1;
    }
  }

	// Interval halve
  if (i > 10)
    xb = 0.0;
  else
    xb = SSHalve (a0, a1, x0, x1, kc, kb, xcaq, xbaq);

  return (xb);
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::SSF(LDBLE xb, LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq)
{
	//
	// Need root of this function to determine xb
	//
  LDBLE lb, lc, f, xc, r;

  xc = 1 - xb;

  if (xb == 0)
    xb = (LDBLE)1e-20;

  if (xc == 0)
    xc = (LDBLE)1e-20;

  lc = exp ((a0 - a1 * (-4 * xb + 3)) * xb * xb);
  lb = exp ((a0 + a1 * (4 * xb - 1)) * xc * xc);

  r = lc * kc / (lb * kb);
  
	f = xcaq * (xb / r + xc) + xbaq * (xb + r * xc) - 1;
  
	return (f);
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::SSHalve(LDBLE a0, LDBLE a1, LDBLE x0, LDBLE x1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq)
{
  int i;
  LDBLE x, y0, dx, y;

  y0 = SSF (x0, a0, a1, kc, kb, xcaq, xbaq);
  dx = (x1 - x0);

	// Loop for interval halving
  for (i = 0; i < 100; i++)
  {
    dx *= 0.5;
    x = x0 + dx;
    y = SSF(x, a0, a1, kc, kb, xcaq, xbaq);

    if (dx < 1e-8 || y == 0)
      break;

		if (y0 * y >= 0)
    {
      x0 = x;
      y0 = y;
    }
  }

  return (x0 + dx);
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::CalcPSIAvg(LDBLE surf_chrg_eq)
{
	//
	// calculate the average (F * Psi / RT) that lets the DL charge counter the surface charge
	//

  int i, iter, count_g;
  LDBLE fd, fd1, p, temp, ratio_aq;

	count_g = md->surface_charge_ptr->g->Count();
  ratio_aq = md->surface_charge_ptr->mass_water / md->mass_water_aq_x;
  p = 0;

  if (surf_chrg_eq == 0)
    return (0.0);
  else if (surf_chrg_eq < 0)
    p = (LDBLE)-0.5 * log (-surf_chrg_eq * ratio_aq / md->mu_x + (LDBLE)1);
  else if (surf_chrg_eq > 0)
    p = (LDBLE)0.5 * log (surf_chrg_eq * ratio_aq / md->mu_x + (LDBLE)1);

	// Optimize p in SS{s_x[i]->moles * z_i * g(p)} = -surf_chrg_eq
	// g(p) = exp(-p * z_i) * ratio_aq
	// Elsewhere in PHREEQC, g is the excess, after subtraction of conc's for p = 0:
	//                       g(p) = (exp(-p *z_i) - 1) * ratio_aq
  iter = 0;
	ChargeGroup *cg;

  do
  {
    fd = surf_chrg_eq;
    fd1 = 0.0;

    for (i = 1; i < count_g; i++)
    {
			cg = gd->charge_group[i];

			if (md->use.sur_p->type == CD_MUSIC)
				temp = exp (-cg->z * p);
      else // multiply with ratio_aq for multiplier options cp and cm in calc_all_donnan (not used now)...  */
				temp = exp (-cg->z * p) * ratio_aq;

      if (md->use.sur_p->only_counter_ions && ((surf_chrg_eq < 0 && cg->z < 0) || (surf_chrg_eq > 0 && cg->z > 0)))
				temp = 0.0;

      fd += cg->eq * temp;
      fd1 -= cg->z * cg->eq * temp;
    }

    fd /= -fd1;
    p += (fd > 1) ? 1 : ((fd < -1) ? -1 : fd);

    if (fabs (p) < md->G_TOL)
      p = 0.0;

    iter++;
    if (iter > 50)
			throw exception(); //ToDo: create exceptions manager
  }
  while (fabs (fd) > 1e-12 && p != 0.0);

  return (p);
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SSBinary (SS *s_s_ptr)
{
  LDBLE nb, nc, n_tot, xb, xc, dnb, dnc, a0, a1;
  LDBLE xb2, xb3, xb4, xc2, xc3;
  LDBLE xb1, xc1;
	SSComp *ss0, *ss1;

	// component 0 is major component
	// component 1 is minor component
	// xb is the mole fraction of second component (formerly trace)
	// xc is the mole fraction of first component (formerly major)

	// Calculate mole fractions and log lambda and derivative factors
  n_tot = s_s_ptr->total_moles;

	ss0 = (*s_s_ptr->comps_list)[0];
	ss1 = (*s_s_ptr->comps_list)[1];

  nc = ss0->moles;
  xc = nc / n_tot;
  nb = ss1->moles;
  xb = nb / n_tot;

  a0 = s_s_ptr->a0;
  a1 = s_s_ptr->a1;

  if (s_s_ptr->miscibility && xb > s_s_ptr->xb1 && xb < s_s_ptr->xb2)
  {
		// In miscibility gap
    xb1 = s_s_ptr->xb1;
    xc1 = (LDBLE)1.0 - xb1;

    ss0->fraction_x = xc1;
    ss0->log10_fraction_x = log10 (xc1);
    ss0->phase->log10_fraction_x = ss0->log10_fraction_x;

    ss1->fraction_x = xb1;
    ss1->log10_fraction_x = log10 (xb1);
    ss1->phase->log10_fraction_x = ss1->log10_fraction_x;

    ss0->log10_lambda = xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)) / md->LOG_10;
    ss0->phase->log10_lambda = ss0->log10_lambda;

    ss1->log10_lambda = xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)) / md->LOG_10;
    ss1->phase->log10_lambda = ss1->log10_lambda;

    ss0->dnb = 0;
    ss0->dnc = 0;
    ss1->dnb = 0;
    ss1->dnc = 0;
    ss0->phase->dnb = 0;
    ss0->phase->dnc = 0;
    ss1->phase->dnb = 0;
    ss1->phase->dnc = 0;
  }
  else
  {
		// Not in miscibility gap
    ss0->fraction_x = xc;
    ss0->log10_fraction_x = log10 (xc);
    ss0->phase->log10_fraction_x =
      ss0->log10_fraction_x;

    ss1->fraction_x = xb;
    ss1->log10_fraction_x = log10 (xb);
    ss1->phase->log10_fraction_x = ss1->log10_fraction_x;

    ss0->log10_lambda = xb * xb * (a0 - a1 * (3 - 4 * xb)) / md->LOG_10;
    ss0->phase->log10_lambda = ss0->log10_lambda;

    ss1->log10_lambda = xc * xc * (a0 + a1 * (4 * xb - 1)) / md->LOG_10;
    ss1->phase->log10_lambda = ss1->log10_lambda;

    xc2 = xc * xc;
    xc3 = xc2 * xc;
    xb2 = xb * xb;
    xb3 = xb2 * xb;
    xb4 = xb3 * xb;

		// used derivation that did not substitute x2 = 1-x1

    // first component, df1/dn1
    dnc = 2 * a0 * xb2 + 12 * a1 * xc * xb2 + 6 * a1 * xb2;
    ss0->phase->dnc = -xb / nc + dnc / n_tot;

    // first component, df1/dn2
    dnb = 1 - 2 * a0 * xb + 2 * a0 * xb2 + 8 * a1 * xc * xb - 12 * a1 * xc * xb2 - 2 * a1 * xb + 2 * a1 * xb2;
    ss0->phase->dnb = dnb / n_tot;

    // second component, df2/dn1
    dnc = 1 - 2 * a0 * xc + 2 * a0 * xc2 - 8 * a1 * xb * xc + 12 * a1 * xb * xc2 + 2 * a1 * xc - 2 * a1 * xc2;
    ss1->phase->dnc = dnc / n_tot;

    // second component, df2/dn2
    dnb = 2 * a0 * xc2 + 12 * a1 * xb * xc2 - 6 * a1 * xc2;
    ss1->phase->dnb = -xc / nb + dnb / n_tot;
  }
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SSIdeal(SS *s_s_ptr)
{
  int k, j;
  LDBLE n_tot, n_tot1;

	// component 0 is major component
	// component 1 is minor component
	// xb is the mole fraction of second component (formerly trace)
	// xc is the mole fraction of first component (formerly major)

	// Calculate mole fractions and log lambda and derivative factors
  n_tot = s_s_ptr->total_moles;

	// Ideal solid solution
  s_s_ptr->dn = (LDBLE)1.0 / n_tot;

	SSComp *ssc_k, *ssc_j;

	for (k = 0; k < s_s_ptr->comps_list->Count(); k++)
  {
		ssc_k = (*s_s_ptr->comps_list)[k];

    n_tot1 = 0;
		for (j = 0; j < s_s_ptr->comps_list->Count(); j++)
    {
			ssc_j = (*s_s_ptr->comps_list)[j];

			if (j != k)
				n_tot1 += ssc_j->moles;
    }

    ssc_k->log10_lambda = 0;
    ssc_k->phase->log10_lambda = 0;

    ssc_k->dnb = -(n_tot1) / ssc_k->moles * n_tot;
    ssc_k->phase->dnb = ssc_k->dnb;

    ssc_k->dn = s_s_ptr->dn;
    ssc_k->phase->dn = s_s_ptr->dn;
  }
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::CalcInitG()
{
  int i, j, k;
  int count_g, count_charge;
	LDBLE la;

  if (md->use.sur_p == NULL)
    return true;

	// calculate g for each surface
	Unknown *u;
	SurfaceDiffLayer *new_g, *g_p;
	SpeciesDiffLayer *dl_p;
	Species *s;

  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
		u = (*unknown_list)[j];

    if (u->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = u->surface_charge;

    count_g = 0;
    if (u->surface_charge->g != NULL)
			count_g = u->surface_charge->g->Count();

		if (count_g == 0)
    {
			new_g = u->surface_charge->g->AddNew();

      new_g->charge = 0.0;
      new_g->g = 0.0;
      new_g->dg = 0.0;

			la = (*u->master)[0]->s->la;
      md->xd = exp (-2 * la * md->LOG_10);

			md->alpha =	sqrt (EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * (LDBLE)1000.0 * md->tk_x * (LDBLE)0.5);
    }

		// calculate g for given surface for each species
    count_g = 1;
		for (i = 0; i < md->s_x->Count(); i++)
    {
			s = (*md->s_x)[i];

      if (s->type > HPLUS)
				continue;

      for (k = 0; k < count_g; k++)
      {
				g_p = (*u->surface_charge->g)[k];

				if (Equal(g_p->charge, s->z, md->G_TOL))
				{
					dl_p = (*s->diff_layer)[count_charge];

					dl_p->charge = u->surface_charge;
					dl_p->count_g = k;
					dl_p->g_moles = 0.0;
					dl_p->dg_g_moles = 0.0;
					break;
				}
      }

      if (k >= count_g)
      {
				g_p = u->surface_charge->g->AddNew();

				// save g for charge
				g_p->charge = s->z;
	
				if (u->surface_charge->grams > 0.0)
				{
					g_p->g = 2 * md->alpha * sqrt (md->mu_x) * (pow ((LDBLE)md->xd, (LDBLE)(s->z / 2.0)) - 1) * md->surface_charge_ptr->grams * md->surface_charge_ptr->specific_area / F_C_MOL;
					g_p->dg = -s->z;

					if ((md->use.sur_p->only_counter_ions) && g_p->g < 0)
					{
						g_p->g = 0;
						g_p->dg = 0;
					}
				}
				else
				{
					g_p->g = 0.0;
					g_p->dg = -s->z;
				}

				/* save g for species */
				dl_p = (*s->diff_layer)[count_g];

				dl_p->charge = u->surface_charge;
				dl_p->count_g = count_g;
				dl_p->g_moles = 0.0;
				dl_p->dg_g_moles = 0.0;
				count_g++;
      }
    }

		count_charge++;
		u->surface_charge->g->SetNewCapacity(count_g);
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::CalcAllG()
{
  int i, j, k;
  bool converge, converge1;
  int count_g, count_charge;
  LDBLE new_g, xd1;
  LDBLE epsilon, la;
	Unknown *u;
	SurfaceDiffLayer *sdl_p;
	SpeciesDiffLayer *spdl_p;
	Species *s;

  if (md->use.sur_p == NULL)
    return (OK);

	// calculate g for each surface
  epsilon = md->convergence_tolerance;

	if (md->convergence_tolerance >= 1e-8)
    md->G_TOL = (LDBLE)1e-9;
  else
    md->G_TOL = (LDBLE)1e-10;

  converge = true;
  count_charge = 0;

  for (j = 0; j < count_unknowns; j++)
  {
		u = (*unknown_list)[j];

    if (u->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = u->surface_charge;

		sdl_p = (*u->surface_charge->g)[0];
    sdl_p->charge = 0.0;
    sdl_p->g = 0.0;
    sdl_p->dg = 0.0;

		count_g = 1;

		la = (*u->master)[0]->s->la;

    md->xd = exp (-2 * la * md->LOG_10);
    md->alpha = sqrt (EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * (LDBLE)1000.0 * md->tk_x * (LDBLE)0.5);

		// calculate g for given surface for each species
    for (i = 0; i < md->s_x->Count(); i++)
    {
			s = (*md->s_x)[i];

      if (s->type > HPLUS)
				continue;

      for (k = 0; k < count_g; k++)
      {
				sdl_p = (*u->surface_charge->g)[k];

				if (Equal(sdl_p->charge, s->z, md->G_TOL))
				{
					spdl_p = (*s->diff_layer)[count_charge];

					spdl_p->charge = u->surface_charge;
					spdl_p->count_g = k;
					break;
				}
      }

      if (k < count_g)
				continue;

      if (u->surface_charge->grams > 0.0)
      {
				md->z = s->z;

				if (
						(md->use.sur_p->only_counter_ions) || 
						(((la > 0) && (md->z < 0)) || 
						((la < 0) && 
						(md->z > 0))))
				{
					if (md->xd > 0.1)
						new_g = QRombMidPnt (1.0, md->xd);
					else if (md->xd > 0.01)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, md->xd);
					else if (md->xd > 0.001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, md->xd);
					else if (md->xd > 0.0001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, md->xd);
					else if (md->xd > 0.00001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, (LDBLE).0001) + QRombMidPnt((LDBLE)0.0001, md->xd);
					else if (md->xd > 0.000001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, (LDBLE).0001) + QRombMidPnt((LDBLE)0.0001, (LDBLE).00001) + QRombMidPnt((LDBLE)0.00001, md->xd);
					else if (md->xd > 0.0000001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, (LDBLE).0001) + QRombMidPnt((LDBLE)0.0001, (LDBLE).00001) + QRombMidPnt ((LDBLE)0.00001, (LDBLE).000001) + QRombMidPnt((LDBLE)0.000001, md->xd);
					else if (md->xd > 0.00000001)
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, (LDBLE).0001) + QRombMidPnt((LDBLE)0.0001, (LDBLE).00001) + QRombMidPnt ((LDBLE)0.00001, (LDBLE).000001) + QRombMidPnt((LDBLE)0.000001, (LDBLE).0000001) + QRombMidPnt((LDBLE)0.0000001, md->xd);
					else
						new_g = QRombMidPnt((LDBLE)1.0, (LDBLE)0.1) + QRombMidPnt((LDBLE)0.1, (LDBLE)0.01) + QRombMidPnt((LDBLE)0.01, (LDBLE).001) + QRombMidPnt((LDBLE)0.001, (LDBLE).0001) + QRombMidPnt((LDBLE)0.0001, (LDBLE).00001) + QRombMidPnt ((LDBLE)0.00001, (LDBLE).000001) + QRombMidPnt((LDBLE)0.000001, (LDBLE).0000001) + QRombMidPnt((LDBLE)0.0000001, (LDBLE).00000001) + QRombMidPnt((LDBLE)0.00000001, md->xd);
				}
				else
					new_g = 0;
      }
      else
				new_g = 0.0;

      if (md->use.sur_p->only_counter_ions && new_g < 0)
				new_g = 0;

			sdl_p = (*u->surface_charge->g)[count_g];

      sdl_p->charge = s->z;
      converge1 = true;

      if ((fabs (new_g) >= 1.) && (fabs ((new_g - sdl_p->g) / new_g) > epsilon))
				converge1 = false;
      else if (fabs (new_g - sdl_p->g) > epsilon)
				converge1 = false;

			if (converge1 == false)
				converge = false;

      sdl_p->g = new_g;

			if (new_g == 0)
				sdl_p->dg = 0;
      else
      {
				if (u->surface_charge->grams > 0.0)
				{
					sdl_p->dg = md->surface_charge_ptr->grams * md->surface_charge_ptr->specific_area * md->alpha * GFunction(md->xd) / F_C_MOL;
					sdl_p->dg *= (LDBLE)-2. / (exp (la * md->LOG_10) * exp (la * md->LOG_10));

					if ((md->xd - 1) < 0.0)
						sdl_p->dg *= -1.0;

					if (fabs (sdl_p->dg) < 1e-8)
					{
						xd1 = exp ((LDBLE)-2 * (LDBLE)1e-3 * md->LOG_10);
						new_g = QRombMidPnt(1.0, xd1);
						sdl_p->dg = new_g / (LDBLE).001;
					}
				}
				else
					sdl_p->dg = 0.0;
      }

			spdl_p = (*s->diff_layer)[count_charge];

			spdl_p->charge = u->surface_charge;
      spdl_p->count_g = count_g;
      count_g++;
    }

    count_charge++;
  }

  return (converge);
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::QRombMidPnt(LDBLE x1, LDBLE x2)
{
  LDBLE ss, dss;
  LDBLE sv[MAX_QUAD + 2], h[MAX_QUAD + 2];
  int j;

  h[0] = 1.0;
  sv[0] = MidPnt (x1, x2, 1);
  for (j = 1; j < MAX_QUAD; j++)
  {
    sv[j] = MidPnt (x1, x2, j + 1);
    h[j] = h[j - 1] / (LDBLE)9.0;

    if (fabs (sv[j] - sv[j - 1]) <= md->G_TOL * fabs (sv[j]))
    {
      sv[j] *= md->surface_charge_ptr->grams * md->surface_charge_ptr->specific_area * md->alpha / F_C_MOL;

      if ((x2 - 1) < 0.0)
				sv[j] *= -1.0;

      return (sv[j]);
    }

    if (j >= K_POLY - 1)
    {
      Polint(&h[j - K_POLY], &sv[j - K_POLY], K_POLY, 0.0, &ss, &dss);

      if (fabs (dss) <= md->G_TOL * fabs (ss) || fabs (dss) < md->G_TOL)
      {
				ss *= md->surface_charge_ptr->grams * md->surface_charge_ptr->specific_area * md->alpha / F_C_MOL;	

				if ((x2 - 1) < 0.0)
				  ss *= -1.0;
	
				return (ss);
      }
    }
  }

	throw exception();
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::MidPnt(LDBLE x1, LDBLE x2, int n)
{
  LDBLE xv, tnm, sum, del, ddel;
  static LDBLE sv;
  int it, j;

  if (n == 1)
  {
    sv = (x2 - x1) * GFunction ((LDBLE)0.5 * (x1 + x2));
    return (sv);
  }
  else
  {
    for (it = 1, j = 1; j < n - 1; j++)
      it *= 3;

    tnm = (LDBLE) it;
    del = (x2 - x1) / (3 * tnm);
    ddel = del + del;
    xv = x1 + (LDBLE)0.5 * del;
    sum = 0.0;
    
		for (j = 1; j <= it; j++)
    {
      sum += GFunction (xv);
      xv += ddel;
      sum += GFunction (xv);
      xv += del;
    }

		sv = (sv + (x2 - x1) * sum / tnm) / (LDBLE)3.0;
    return sv;
  }
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::Polint(LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv, LDBLE * dy)
{
  int i, m, ns;
  LDBLE den, dif, dift, ho, hp, w;
  LDBLE *c, *d;

  ns = 1;
  dif = fabs (xv - xa[1]);

	c = NULL;
	d = NULL;

	try
	{
		c = new LDBLE [n + 1];
		d = new LDBLE [n + 1];
	}
	catch(...)
	{
		delete [] c;
		delete [] d;

		throw;
	}

  for (i = 1; i <= n; i++)
  {
    dift = fabs (xv - xa[i]);

    if (dift < dif)
    {
      ns = i;
      dif = dift;
    }
    
		c[i] = ya[i];
    d[i] = ya[i];
  }

  *yv = ya[ns--];

  for (m = 1; m < n; m++)
  {
    for (i = 1; i <= n - m; i++)
    {
      ho = xa[i] - xv;
      hp = xa[i + m] - xv;
      w = c[i + 1] - d[i];

      if ((den = ho - hp) == 0.0)
				throw exception();

			den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }

    if (2 * ns < (n - m))
      *dy = c[ns + 1];
    else
      *dy = d[ns--];

		*yv += *dy;
  }

	delete [] c;
	delete [] d;
}
//-----------------------------------------------------------------------------------------------------------
LDBLE ModelEngine::GFunction (LDBLE x_value)
{
  LDBLE sum, return_value;
  int i, j;
  LDBLE ln_x_value;
	SurfaceCharge *sc_p;
	SurfaceDiffLayer *sdl_p;

  if (Equal(x_value, 1.0, md->G_TOL * 100))
    return (0.0);

  sum = 0.0;
  ln_x_value = log(x_value);

	sc_p = (*md->use.sur_p->charge)[0];

	for (j = 0; j < sc_p->g->Count(); j++)
  {
		sdl_p = (*sc_p->g)[j];

    sdl_p->psi_to_z = exp (ln_x_value * sdl_p->charge) - (LDBLE)1.0;
  }

	Species *s;
	for (i = 0; i < md->s_x->Count(); i++)
  {
		s = (*md->s_x)[i];

    if (s->type < H2O && s->z != 0.0)
    {
      for (j = 0; j < sc_p->g->Count(); j++)
      {
				sdl_p = (*sc_p->g)[j];
	
				if (sdl_p->charge == s->z)
				{
					sum += s->moles * sdl_p->psi_to_z;
					break;
				}
      }
    }
  }

  if (sum < 0.0)
		throw exception();

  return_value = (exp (ln_x_value * md->z) - 1) / sqrt ((x_value * x_value * md->mass_water_aq_x * sum));
  return (return_value);
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::CalcInitDonnanMusic()
{
  int i, j, k;
  int count_g, count_charge;
  //char name[MAX_LENGTH];
  LDBLE psi_avg, f_sinh, ratio_aq;
	Species *s;
	Unknown *u;
	SurfaceDiffLayer *sdl_p;
	SpeciesDiffLayer *spdl_p;

  if (md->use.sur_p == NULL)
    return true;

  f_sinh = sqrt ((LDBLE)8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * md->tk_x);

  if (md->convergence_tolerance >= 1e-8)
  {
    md->G_TOL = (LDBLE)1e-9;
  }
  else
  {
    md->G_TOL = (LDBLE)1e-10;
  }

	// sum eq of each charge number in solution...
	gd->charge_group.Clear();

	ChargeGroup *cg_p = gd->charge_group.AddNew();

	cg_p->z = 0.0;
  cg_p->eq = 0.0;

  count_g = 1;
	for (i = 0; i < md->s_x->Count(); i++)
  {
		s = (*md->s_x)[i];

    if (s->type > HPLUS)
      continue;

    for (k = 0; k < count_g; k++)
    {
      if (Equal(cg_p->z, s->z, md->G_TOL))
      {
				cg_p->eq += s->z * s->moles;
				break;
      }
    }

    if (k >= count_g)
    {
			cg_p = gd->charge_group.AddNew();

			cg_p->z = s->z;
      cg_p->eq = s->z * s->moles;

      count_g++;
    }
  }

	// calculate g for each surface...
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
		u = (*unknown_list)[j];

    if (u->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = u->surface_charge;

    u->surface_charge->g->Clear();
		u->surface_charge->g->AddNew(count_g);

    // find psi_avg that matches surface charge...
		psi_avg = CalcPSIAvg(0);

    // fill in g's 
    ratio_aq = md->surface_charge_ptr->mass_water / md->mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
			sdl_p = (*u->surface_charge->g)[k];
			cg_p = gd->charge_group[k];

      sdl_p->charge = cg_p->z;
      sdl_p->g = exp (cg_p->z * psi_avg) - 1;

      // save g for species
			for (i = 0; i < md->s_x->Count(); i++)
      {
				s = (*md->s_x)[i];

				if (Equal(cg_p->z, s->z, md->G_TOL))
				{
					spdl_p = (*s->diff_layer)[count_charge];

					spdl_p->charge = u->surface_charge;
					spdl_p->count_g = k;
					spdl_p->g_moles = 0.0;
					spdl_p->dg_g_moles = 0.0;
				}
      }
    }

    count_charge++;
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::GasPhaseCheck(GasPhase *gas_phase_ptr)
{
	//
	//   Check for missing elements
	//

  int i, j;
	LDBLE coef;

  GasComp *gas_comp_ptr;
  Master *master_ptr;

  if (gas_phase_ptr == NULL)
    return true;

	// Check that all elements are in solution for phases with zero mass
	coef = 1.0;
	for (i = 0; i < gas_phase_ptr->comps->Count(); i++)
  {
		gas_comp_ptr = (*gas_phase_ptr->comps)[i];

    eos_list->Clear();
    parent_count = 0;

    if (gas_comp_ptr->moles <= 0.0)
    {
      CopyToTempEOSList(gas_comp_ptr->phase->eos_list, coef);

      for (j = 0; j < eos_list->Count(); j++)
      {
				master_ptr = (*eos_list)[j]->e->primary;

				if (master_ptr->s == gd->s_hplus)
					continue;
				else if (master_ptr->s == gd->s_h2o)
					continue;
				else if (master_ptr->total > MIN_TOTAL)
					continue;
				else
				{
					return false;
				}
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::PPAssemblageCheck (PPAssemblage *pp_assemblage_ptr)
{
	//
	// Check for missing elements
	//

  int i, j, k;
  char token[MAX_LENGTH];
  char *ptr;
  PurePhase *pure_phase_ptr;
  Master *master_ptr;

  if (CheckPPAssemblage(pp_assemblage_ptr))
    return true;

	// Check that all elements are in solution for phases with zero mass
  for (j = 0; j < pp_assemblage_ptr->pure_phases->Count(); j++)
  {
		pure_phase_ptr = (*pp_assemblage_ptr->pure_phases)[j];

		eos_list->Clear();
    parent_count = 0;

    if (pure_phase_ptr->moles <= 0.0)
    {
      pure_phase_ptr->delta = 0.0;

			if (!pure_phase_ptr->add_formula.IsEmpty())
      {
				pure_phase_ptr->add_formula.Copy(token);
				ptr = &(token[0]);
				
				if (!GetElementsInSpecies(&ptr, 1.0))
					return false;
      }
      else
				CopyToTempEOSList(pure_phase_ptr->phase->eos_list, 1.0);

			for (i = 0; i < eos_list->Count(); i++)
      {
				master_ptr = (*eos_list)[i]->e->primary;

				if (master_ptr->s == gd->s_hplus)
					continue;
				else if (master_ptr->s == gd->s_h2o)
					continue;
				else if (master_ptr->total > MIN_TOTAL)
					continue;
				else
				{
					sprintf (error_string,
									 "Element %s is contained in %s (which has 0.0 mass),"
									 "\t\nbut is not in solution or other phases.",
									 master_ptr->name.CharPtr(), pure_phase_ptr->phase->name.CharPtr());
					//printf ("\nPHREEQC WARNING: %s\n", error_string);
					//ToDo: Preparar para colocar a mensagem
					
					// Make la's of all master species for the element small, so SI will be small
					// and no mass transfer will be calculated
					for (k = 0; k < gd->master_list.Count(); k++)
					{
						Master *m = gd->master_list[k];
						if (m->e->primary == master_ptr)
						{
							m->s->la = (LDBLE)-9999.999;
						}
					}
				}
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SSAssemblageCheck(SSAssemblage *s_s_assemblage_ptr)
{
	//
	// Check for missing elements
	//

  int i, j, l;
  Master *master_ptr;
	SS *ss_p;
	SSComp *ssc_p;

  if (s_s_assemblage_ptr == NULL)
    return true;

	// Check that all elements are in solution for phases with zero mass
	for (i = 0; i < s_s_assemblage_ptr->ss_list->Count(); i++)
  {
		ss_p = (*s_s_assemblage_ptr->ss_list)[i];

		for (j = 0; j < ss_p->comps_list->Count(); j++)
    {
			ssc_p = (*ss_p->comps_list)[j];

      eos_list->Clear();
      parent_count = 0;

      if (ssc_p->moles <= 0.0)
      {
				CopyToTempEOSList(ssc_p->phase->eos_list, 1.0);

				for (l = 0; l < eos_list->Count(); l++)
				{
					master_ptr = (*eos_list)[l]->e->primary;

					if (master_ptr->s == gd->s_hplus)
						continue;
					else if (master_ptr->s == gd->s_h2o)
						continue;
					else if (master_ptr->total > MIN_TOTAL_SS)
						continue;
					else
						return false;
				}
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::CalcInitDonnan()
{
  int i, j, k;
  int count_g, count_charge;
  //char name[MAX_LENGTH];
  LDBLE f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq, la;
	SurfaceDiffLayer *sdl_p;
	SpeciesDiffLayer *spdl_p;

  if (md->use.sur_p == NULL)
    return true;

  if (md->use.sur_p->type == CD_MUSIC)
    return CalcInitDonnanMusic();

  f_sinh = sqrt ((LDBLE)8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * md->tk_x * md->mu_x);

  if (md->convergence_tolerance >= 1e-8)
    md->G_TOL = (LDBLE)1e-9;
  else
    md->G_TOL = (LDBLE)1e-13;

	// sum eq of each charge number in solution...
  gd->charge_group.Clear();

	ChargeGroup *cg_p = gd->charge_group.AddNew();

	cg_p->z = 0.0;
  cg_p->eq = 0.0;

  count_g = 1;

	Species *s;
  for (i = 0; i < md->s_x->Count(); i++)
  {
		s = (*md->s_x)[i];

    if (s->type > HPLUS)
      continue;

    for (k = 0; k < count_g; k++)
    {
      if (Equal(cg_p->z, s->z, md->G_TOL))
      {
				cg_p->eq += s->z * s->moles;
				break;
      }
    }

    if (k >= count_g)
    {
			cg_p = gd->charge_group.AddNew();

			cg_p->z = s->z;
      cg_p->eq = s->z * s->moles;

      count_g++;
    }
  }

	// calculate g for each surface...
  count_charge = 0;
	Unknown *u;
  for (j = 0; j < count_unknowns; j++)
  {
		u = (*unknown_list)[j];

    if (u->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = u->surface_charge;

		u->surface_charge->g->Clear();
		u->surface_charge->g->AddNew(count_g);

    // find surface charge from potential...
    A_surf = u->surface_charge->specific_area * u->surface_charge->grams;

		la = (*u->master)[0]->s->la; 

    f_psi = la * md->LOG_10;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;

    // find psi_avg that matches surface charge...
		psi_avg = CalcPSIAvg(0);

    // fill in g's
    ratio_aq = md->surface_charge_ptr->mass_water / md->mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
			sdl_p = (*u->surface_charge->g)[k];
			cg_p = gd->charge_group[k];
      
			sdl_p->charge = cg_p->z;
      sdl_p->g = ratio_aq * (exp (-cg_p->z * psi_avg) - 1);

      if (md->use.sur_p->only_counter_ions && ((surf_chrg_eq < 0 && cg_p->z < 0) || (surf_chrg_eq > 0 && cg_p->z > 0)))
				sdl_p->g = -ratio_aq;

      if (sdl_p->g != 0)
				sdl_p->dg = -A_surf * f_sinh * cosh (f_psi) / (cg_p->eq * F_C_MOL);
      else
				sdl_p->dg = -cg_p->z;

			// save g for species
      for (i = 0; i < md->s_x->Count(); i++)
      {
				s = (*md->s_x)[i];

				if (Equal(cg_p->z, s->z, md->G_TOL))
				{
					spdl_p = (*s->diff_layer)[count_charge];

					spdl_p->charge = u->surface_charge;
					spdl_p->count_g = k;
					spdl_p->g_moles = 0.0;
					spdl_p->dg_g_moles = 0.0;
				}
      }
    }

		count_charge++;
  }

	return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::CalcAllDonnanMusic()
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge;
	bool converge;
  //char name[MAX_LENGTH];
  LDBLE new_g, f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;
  LDBLE cz, cm, cp, la;

	if (md->use.sur_p == NULL)
    return true;

  f_sinh = sqrt ((LDBLE)8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * (LDBLE)1000.0) * md->tk_x * md->mu_x);
  cz = cm = 1.0;
  cp = 1.0;

	// calculate g for each surface...
  converge = true;
  count_charge = 0;
	Unknown *u;
  for (j = 0; j < count_unknowns; j++)
  {
		u = (*unknown_list)[j];

    if (u->type != SURFACE_CB)
      continue;

    md->surface_charge_ptr = u->surface_charge;

		// sum eq of each charge number in solution...
		count_g = u->surface_charge->g->Count();

    for (i = 0; i < count_g; i++)
      gd->charge_group[i]->eq = 0.0;

		Species *s_x;
    for (i = 0; i < md->s_x->Count(); i++)
    {
			s_x = (*md->s_x)[i];

      if (s_x->type > HPLUS)
				continue;

      for (k = 0; k < count_g; k++)
      {
				if (Equal(gd->charge_group[k]->z, s_x->z, md->G_TOL))
				{
					gd->charge_group[k]->eq += s_x->z * s_x->moles;
					break;
				}
      }
    }

    // find surface charge from potential...
    A_surf = u->surface_charge->specific_area * u->surface_charge->grams;
		la = (*(*unknown_list)[j + 2]->master)[0]->s->la;
    f_psi = la * md->LOG_10;	// -FPsi/RT
    f_psi = f_psi / 2;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;
    psi_avg = CalcPSIAvg(surf_chrg_eq);

    // fill in g's
    ratio_aq = md->surface_charge_ptr->mass_water / md->mass_water_aq_x;

		SurfaceDiffLayer *sdl_p;
    for (k = 0; k < count_g; k++)
    {
			sdl_p = (*u->surface_charge->g)[k];

      sdl_p->charge = gd->charge_group[k]->z;
      new_g = exp (gd->charge_group[k]->z * psi_avg) - 1;

      if (new_g < -ratio_aq)
				new_g = -ratio_aq + md->G_TOL * (LDBLE)1e-5;

      if (fabs (new_g) >= 1)
      {
				if (fabs ((new_g - sdl_p->g) / new_g) > md->convergence_tolerance)
				  converge = FALSE;
      }
      else
      {
				if (fabs (new_g - sdl_p->g) > md->convergence_tolerance)
					converge = FALSE;
      }

			sdl_p->g = new_g;

      // save g for species
			SpeciesDiffLayer *spdl_p;
      for (i = 0; i < md->s_x->Count(); i++)
      {
				s_x = (*md->s_x)[i];
				spdl_p = (*s_x->diff_layer)[count_charge];

				
				if (Equal(gd->charge_group[k]->z, s_x->z, md->G_TOL))
				{
					spdl_p->charge = u->surface_charge;
					spdl_p->count_g = k;
				}
      }
    }

    count_charge++;
  }

  return (converge);
}
//-----------------------------------------------------------------------------------------------------------
CONVERGE_RESULT ModelEngine::ExecuteSurfaceModel()
{
	//
	// Use extra iterative loop to converge on g_factors
	//

  int i;
	bool result;
  //struct solution *solution_ptr;
  LDBLE prev_aq_x;

	// Allocate space for g factors for diffuse layer in surface complexation
  if (md->dl_type_x != NO_DL)
  {
		Species *s;
		for (i = 0; i < gd->species_list.Count(); i++)
    {
			s = gd->species_list[i];
			s->diff_layer->Clear();
    }
  }

	SurfaceCharge *sc_p;
  if (md->state >= REACTION && md->dl_type_x != NO_DL)
  {
    md->mass_water_bulk_x = md->use.sol_p->mass_water;

		for (i = 0; i < md->use.sur_p->charge->Count(); i++)
    {
			sc_p = (*md->use.sur_p->charge)[i];

			md->mass_water_bulk_x += sc_p->mass_water;

      if (md->use.sur_p->debye_lengths > 0)
				sc_p->mass_water = 0.0;
    }
  }

  if (!Prepare())
		return ERROR_CR;

  if (!KTemp(md->tc_x))
		return ERROR_CR;

  if (md->use.sur_p->dl_type == DONNAN_DL)
  {
    if (!InitialSurfaceWater ())
			return ERROR_CR;

    if (!CalcInitDonnan ())
			return ERROR_CR;
  }
  else
    if (!CalcInitG ())
			return ERROR_CR;

	/*
	if (md->state >= REACTION && md->use.sur_p->new_def == FALSE)
  {
    set (FALSE);
  }
  else
  {
    set (TRUE);
  }
	*/

	Set(true);

  if (ExecuteModel() == ERROR_CR)
    return (ERROR_CR);

  md->g_iterations = 0;

  if (md->use.sur_p->dl_type == DONNAN_DL)
  {
    do
    {
      md->g_iterations++;
      prev_aq_x = md->mass_water_aq_x;

      if(!KTemp(md->tc_x))
				return  ERROR_CR;

			result = Gammas(md->mu_x);
			//d.PrintGammasToFile("ExecuteSurfaceModel-1", md->s_x);
      if(!result)
				return ERROR_CR;

      if(!Molalities(true))
				return ERROR_CR;

      if(!MBSums ())
				return ERROR_CR;

      if(ExecuteModel() == ERROR_CR)
				return (ERROR_CR);

      if (!md->use.sur_p->related_phases) //&& !md->use.sur_p->related_rate)
				if (!InitialSurfaceWater())
					return ERROR_CR;
    } 
		while ((!CalcAllDonnan() || fabs(1 - prev_aq_x / md->mass_water_aq_x) > 1e-6) && md->g_iterations < md->itmax);
  }
  else
  {
    do
    {
      md->g_iterations++;

      if (!KTemp(md->tc_x))
				return ERROR_CR;

			result = Gammas(md->mu_x);
			//d.PrintGammasToFile("ExecuteSurfaceModel-2", md->s_x);
      if (!result)
				return ERROR_CR;

      if (!Molalities(true))
				return ERROR_CR;

      if (!MBSums())
				return ERROR_CR;

      if (ExecuteModel() == ERROR_CR)
				return (ERROR_CR);

      if (!md->use.sur_p->related_phases)
				InitialSurfaceWater();
    }
    while (!CalcAllG() && md->g_iterations < md->itmax);
  }

  if (md->g_iterations >= md->itmax)
		return ERROR_CR;

  return OK_CR;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddSolution(Solution *sol_p, LDBLE extensive, LDBLE intensive)
{
	//
	// Accumulate solution data in master->totals and _x variables.
	//
	// extensive is multiplication factor for solution
	// intensive is fraction of all multiplication factors for all solutions
	//

	int i;
	String value, comment;

	// Add solution to global variables
  md->tc_x += sol_p->tc * intensive;
  md->ph_x += sol_p->ph * intensive;
  md->solution_pe_x += sol_p->solution_pe * intensive;
  md->mu_x += sol_p->mu * intensive;
  md->ah2o_x += sol_p->ah2o * intensive;
  md->density_x += sol_p->density * intensive;

  md->total_h_x += sol_p->total_h * extensive;
  md->total_o_x += sol_p->total_o * extensive;
  md->cb_x += sol_p->cb * extensive;
  md->mass_water_aq_x += sol_p->mass_water * extensive;

#ifdef DEBUG_AddSolution
	d.addsolution_count++;
	fprintf(d.addsolution_f, "%d 1 tc_x: %.20e ph_x: %.20e solution_pe_x: %.20e ", d.addsolution_count, md->tc_x, md->ph_x, md->solution_pe_x);
	fprintf(d.addsolution_f, "mu_x: %.20e ah2o_x: %.20e density_x: %.20e ", md->mu_x, md->ah2o_x, md->density_x);
	fprintf(d.addsolution_f, "total_h_x: %.20e total_o_x: %.20e cb_x: %.20e mass_water_aq_x: %.20e\n", md->total_h_x, md->total_o_x, md->cb_x, md->mass_water_aq_x);
#endif

	Master *m;
	Conc *c;

	// Copy totals data into primary master species
  for (i = 0; i < sol_p->totals->Count(); i++)
  {
		c = (*md->use.sol_p->totals)[i];

		m = MasterPrimarySearch(c->name);

		if (m == NULL)
			return false;

		m->total += c->moles * extensive;

#ifdef DEBUG_AddSolution		
		fprintf(d.addsolution_f, "%d 2 c->name: %s c->moles: %.20e extensive: %.20e m->total: %.20e\n", d.addsolution_count, c->name.CharPtr(), c->moles, extensive, m->total);
#endif
  }

	// Accumulate initial guesses for activities
	MasterActivity *ma;
  for (i = 0; i < md->use.sol_p->ma->Count(); i++)
  {
		ma = (*md->use.sol_p->ma)[i];

		if (!ma->description.IsEmpty())
		{
			m = gd->master_list.Search(&ma->description, true);

			if (m != NULL)
			{
				m->s->la += ma->la * intensive;

#ifdef DEBUG_AddSolution		
				fprintf(d.addsolution_f, "%d 3 ma->description: %s ma->la: %.20e intensive: %.20e m->s->la: %.20e\n", d.addsolution_count, ma->description.CharPtr(), ma->la, intensive, m->s->la);
#endif
			}
		}
  }

	return true;
}
//-----------------------------------------------------------------------------------------------------------
Master *ModelEngine::MasterPrimarySearch(String str)
{
	char *ptr, token[DEFAULT_STR_LENGTH];
	String element;

	str.Copy(token);
	ptr = &token[0];

	// find element name
	GetElement(&ptr, element);

	// return master species
	return gd->master_list.Search(&element, true);
}
//-----------------------------------------------------------------------------------------------------------
/*
bool ModelEngine::AddReaction (Irrev *irrev_ptr, int step_number, LDBLE step_fraction)
{
	//
	// Add irreversible reaction
	//

  int i;
  char c;
  Master *master_ptr;

	// Calculate and save reaction
  if (irrev_ptr->elts->Count() <= 0)
  {
    if (!ReactionCalc(irrev_ptr))
      return false;
  }

	// Step size
  if (!md->incremental_reactions)
  {
    if (irrev_ptr->count_steps > 0)
    {
      if (step_number > irrev_ptr->count_steps)
				md->step_x = irrev_ptr->steps[irrev_ptr->count_steps - 1];
      else
				md->step_x = irrev_ptr->steps[step_number - 1];
    }
    else if (irrev_ptr->count_steps < 0)
    {
      if (step_number > -irrev_ptr->count_steps)
				md->step_x = irrev_ptr->steps[0];
      else
				md->step_x = irrev_ptr->steps[0] * ((LDBLE) step_number) / ((LDBLE) (-irrev_ptr->count_steps));
    }
    else
      md->step_x = 0.0;
  }
  else
  {
    // Incremental reactions
    if (irrev_ptr->count_steps > 0)
    {
      if (step_number > irrev_ptr->count_steps)
				md->step_x = irrev_ptr->steps[irrev_ptr->count_steps - 1];
      else
				md->step_x = irrev_ptr->steps[step_number - 1];
    }
    else if (irrev_ptr->count_steps < 0)
    {
      if (step_number > -irrev_ptr->count_steps)
				md->step_x = 0;
      else
				md->step_x = irrev_ptr->steps[0] / ((LDBLE) (-irrev_ptr->count_steps));
    }
    else
      md->step_x = 0.0;
  }

	// Convert units
	//ToDo: Alterar essa parte
  c = irrev_ptr->units[0];
  if (c == 'm')
  {
    md->step_x *= 1e-3;
  }
  else if (c == 'u')
  {
    md->step_x *= 1e-6;
  }
  else if (c == 'n')
  {
    md->step_x *= 1e-9;
  }

	// Add reaction to totals
	ElementOfSpecies *eos_p;
  for (i = 0; i < irrev_ptr->elts->Count(); i++)
  {
		eos_p = (*irrev_ptr->elts)[i];

    master_ptr = eos_p->e->primary;

    if (master_ptr->s == gd->s_hplus)
      md->total_h_x += eos_p->coef * md->step_x * step_fraction;
    else if (master_ptr->s == gd->s_h2o)
      md->total_o_x += eos_p->coef * md->step_x * step_fraction;
    else
      master_ptr->total += eos_p->coef * md->step_x * step_fraction;
  }

  return true;
}
*/
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddExchange (Exchange *exchange_ptr)
{
	//
	// Accumulate exchange data in master->totals and _x variables.
	//

  int i, j;
  Master *master_ptr;

  if (exchange_ptr == NULL)
    return true;

#ifdef DEBUG_AddExchange
	d.addexchange_count++;
#endif

	// Add element concentrations on exchanger to master species totals
	ExchComp *ec_p;
	ElementOfSpecies *eos_p;
	for (i = 0; i < exchange_ptr->comps->Count(); i++)
  {
		ec_p = (*exchange_ptr->comps)[i];

    for (j = 0; j < ec_p->totals->Count(); j++)
    {
			eos_p = (*ec_p->totals)[j];
      master_ptr = eos_p->e->primary;

      if (master_ptr->s == gd->s_hplus)
			{
				md->total_h_x += eos_p->coef;

#ifdef DEBUG_AddExchange
				fprintf(d.addexchange_f, "%d 1 i: %d j: %d master_ptr->s->name: %s ", d.addexchange_count, i, j, master_ptr->s->name.CharPtr());
				fprintf(d.addexchange_f, "eos_p->coef: %.20e md->total_h_x: %.20e\n", eos_p->coef, md->total_h_x);
#endif
			}
      else if (master_ptr->s == gd->s_h2o)
			{
				md->total_o_x += eos_p->coef;

#ifdef DEBUG_AddExchange
				fprintf(d.addexchange_f, "%d 2 i: %d j: %d master_ptr->s->name: %s ", d.addexchange_count, i, j, master_ptr->s->name.CharPtr());
				fprintf(d.addexchange_f, "eos_p->coef: %.20e md->total_o_x: %.20e\n", eos_p->coef, md->total_o_x);
#endif
			}
      else
			{
				master_ptr->total += eos_p->coef;

#ifdef DEBUG_AddExchange
				fprintf(d.addexchange_f, "%d 3 i: %d j: %d master_ptr->s->name: %s ", d.addexchange_count, i, j, master_ptr->s->name.CharPtr());
				fprintf(d.addexchange_f, "eos_p->coef: %.20e master_ptr->total: %.20e\n", eos_p->coef, master_ptr->total);
#endif
			}

    }
  }

	if (exchange_ptr->new_def)
  {
		for (i = 0; i < gd->master_list.Count(); i++)
		{
			master_ptr = gd->master_list[i];

	    if (master_ptr->type == EX && master_ptr->total > 0)
			{
				master_ptr->s->la = log10 ((LDBLE)0.1 * master_ptr->total);

#ifdef DEBUG_AddExchange
				fprintf(d.addexchange_f, "%d 4 i: %d master_ptr->name: %s master_ptr->total: %.20e master_ptr->s->la: %.20e\n", d.addexchange_count, i, master_ptr->name.CharPtr(), master_ptr->total, master_ptr->s->la);
#endif			
			}
		}
	}
	else
	{
		for (i = 0; i < exchange_ptr->comps->Count(); i++)
    {
			ec_p = (*exchange_ptr->comps)[i];

      ec_p->master->s->la = ec_p->la;
      md->cb_x += ec_p->charge_balance;

#ifdef DEBUG_AddExchange
			fprintf(d.addexchange_f, "%d 4 i: %d ec_p->name: %s ec_p->master->s->la: %.20e ec_p->la: %.20e ec_p->charge_balance: %.20e md->cb_x: %.20e\n", d.addexchange_count, i, ec_p->name.CharPtr(), ec_p->master->s->la, ec_p->la, ec_p->charge_balance, md->cb_x);
#endif	
    }
	}

	return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddSurface(Surface *surface_ptr)
{
	//
	// Accumulate surface data in master->totals and _x variables.
	//

  int i, j;
  Master *master_ptr;

  if (surface_ptr == NULL)
    return true;

	// Add element concentrations on surface to master species totals
  md->dl_type_x = surface_ptr->dl_type;

	SurfaceComp *sc_p;
	ElementOfSpecies *eos_p;
	SurfaceCharge *sch_p;
	for (i = 0; i < surface_ptr->comps->Count(); i++)
  {
		sc_p = (*surface_ptr->comps)[i];

    if (surface_ptr->type == NO_EDL)
      md->cb_x += sc_p->cb;

		if (!surface_ptr->new_def)
      sc_p->master->s->la = sc_p->la;

		// Add surface and specifically sorbed elements
		for (j = 0; j < sc_p->totals->Count(); j++)
    {
			eos_p = (*sc_p->totals)[j];

      master_ptr = eos_p->e->primary;

      if (master_ptr == NULL)
				return false;

			if (master_ptr->s == gd->s_hplus)
				md->total_h_x += eos_p->coef;
      else if (master_ptr->s == gd->s_h2o)
				md->total_o_x += eos_p->coef;
      else
				master_ptr->total += eos_p->coef;
    }
  }

  if (surface_ptr->type != DDL && surface_ptr->type != CD_MUSIC)
    return true;

	for (i = 0; i < surface_ptr->charge->Count(); i++)
  {
		sch_p = (*surface_ptr->charge)[i];

    if (surface_ptr->type == DDL || surface_ptr->type == CD_MUSIC)
      md->cb_x += sch_p->charge_balance;

    if (!surface_ptr->new_def)
    {
      master_ptr = SurfaceGetPSIMaster(sch_p->name, SURF_PSI);
      master_ptr->s->la = sch_p->la_psi;
    }

		// Add diffuse layer elements (including water in Debye layer)
    if (surface_ptr->dl_type != NO_DL && !surface_ptr->new_def)
    {
			for (j = 0; j < sch_p->diffuse_layer_totals->Count(); j++)
      {
				eos_p = (*sch_p->diffuse_layer_totals)[j];

				master_ptr = eos_p->e->primary;

				if (master_ptr->s == gd->s_hplus)
					md->total_h_x += eos_p->coef;
				else if (master_ptr->s == gd->s_h2o)
					md->total_o_x += eos_p->coef;
				else
					master_ptr->total += eos_p->coef;
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddGasPhase(GasPhase *gas_phase_ptr)
{
 
	// Accumulate gas data in master->totals and _x variables.
  int i;

  GasComp *gas_comp_ptr;
  Master *master_ptr;
	ElementOfSpecies *eos_p;

  if (gas_phase_ptr == NULL)
    return true;
 
	// calculate reaction
  eos_list->Clear();
  parent_count = 0;

  for (i = 0; i < gas_phase_ptr->comps->Count(); i++)
  {
		gas_comp_ptr = (*gas_phase_ptr->comps)[i];

		CopyToTempEOSList(gas_comp_ptr->phase->eos_list, gas_comp_ptr->moles);
  }

	// Sort elements in reaction and combine
  CombineElements();

	// Add gas elements to totals
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos_p = (*eos_list)[i];

    master_ptr = eos_p->e->primary;

    if (master_ptr->s == gd->s_hplus)
      md->total_h_x += eos_p->coef;
    else if (master_ptr->s == gd->s_h2o)
      md->total_o_x += eos_p->coef;
    else
      master_ptr->total += eos_p->coef;
  }
  return (OK);
}
//-----------------------------------------------------------------------------------------------------------

/*
bool ModelEngine::ReactionCalc(Irrev *irrev_ptr)
{
	//
	// Go through irreversible reaction initially to
	// determine a list of elements and amounts in 
	// the reaction.
	//

  int i; //, j;
  //LDBLE coef;
	char token[DEFAULT_STR_LENGTH], *ptr;
  Phase *phase_ptr;
	NameCoef *nc_p;
	ElementOfSpecies *eos_p;

	// Go through list and generate list of elements and coefficient of elements in reaction
  eos_list->Clear();
  parent_count = 0;

  for (i = 0; i < irrev_ptr->list->Count(); i++)
  {
		nc_p = (*irrev_ptr->list)[i];
		
		phase_ptr = gd->phase_list.Search(nc_p->name, true);

		// Reactant is a pure phase, copy formula into token
    if (phase_ptr != NULL)
			CopyToTempEOSList(phase_ptr->eos_list, nc_p->coef);
    else
    {
			nc_p->name.Copy(token);
      ptr = &(token[0]);
      GetElementsInSpecies(&ptr, nc_p->coef);
    }
  }

	// Check that all elements are in database
  for (i = 0; i < eos_list->Count(); i++)
  {
		eos_p = (*eos_list)[i];

    if (eos_p->e->master == NULL)
			return false;
  }

  SaveEOSList(irrev_ptr->elts);
  return true;
}
*/
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SumDiffuseLayer(SurfaceCharge *surface_charge_ptr1)
{
  int i, j, count_g;
  LDBLE mass_water_surface;
  LDBLE molality, moles_excess, moles_surface;
	SurfaceCharge *sc_p;
	Species *s;
	SpeciesDiffLayer *sdl_p;
	SurfaceDiffLayer *surdl_p;

  if (md->use.sur_p == NULL)
    return true;

	// Find position of component in list of components
  i = 0;

  for (j = 0; j < md->use.sur_p->charge->Count(); j++)
  {
		sc_p = (*md->use.sur_p->charge)[j];

    if (sc_p == surface_charge_ptr1)
    {
      i = j;
      break;
    }
  }

  if (j >= md->use.sur_p->charge->Count())
		return false;

	// Loop through all surface components, calculate each H2O surface (diffuse layer),
	// H2O aq, and H2O bulk (diffuse layers plus aqueous).
  eos_list->Clear();
  parent_count = 0;

  mass_water_surface = surface_charge_ptr1->mass_water;

  for (j = 0; j < md->s_x->Count(); j++)
  {
		s = (*md->s_x)[j];
		sdl_p = (*s->diff_layer)[i];

    if (s->type > HPLUS)
      continue;

    molality = Under(s->lm);

    count_g = sdl_p->count_g;
		surdl_p = (*surface_charge_ptr1->g)[count_g];

		moles_excess = md->mass_water_aq_x * molality * surdl_p->g;
    moles_surface = mass_water_surface * molality + moles_excess;

		// Accumulate elements in diffuse layer
		CopyToTempEOSList(s->eos_list, moles_surface);
  }

	CopyToTempEOSList(gd->s_h2o->eos_list, mass_water_surface / md->gfw_water);

  CombineElements();

  return true;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::CopyUse(int i)
{
	int index;
	Solution *copy_sol, *ori_sol;
	PPAssemblage *copy_ppa, *ori_ppa;
	Exchange *copy_exc, *ori_exc;
	Surface *copy_sur, *ori_sur;
	GasPhase *copy_gas, *ori_gas;
	SSAssemblage *copy_ssa, *ori_ssa;
	Use *use = &md->use;

	Solution *temp;
	if (use->sol_list->Count() <= i)
	{
		for (index = use->sol_list->Count(); index <= i; index++)
		{
			temp = use->sol_list->AddNew();
			temp->gd = gd;
			temp->InitPE();
		}
	}

	copy_sol = (*use->sol_list)[i];
	ori_sol = (*use->sol_list)[use->sol_number];
	ori_sol->CopyTo(copy_sol);

	// Find pure phase assemblage
	if (use->ppa_in)
  {
		if (use->ppa_list->Count() <= i)
		{
			for (index = use->ppa_list->Count(); index <= i; index++)
				use->ppa_list->AddNew();
		}

		copy_ppa = (*use->ppa_list)[i];
		ori_ppa = (*use->ppa_list)[use->ppa_number];
		ori_ppa->CopyTo(copy_ppa);
  }

	// Find exchange
  if (use->exc_in)
  {
		if (use->exc_list->Count() <= i)
		{
			for (index = use->exc_list->Count(); index <= i; index++)
				use->exc_list->AddNew();
		}

		copy_exc = (*use->exc_list)[i];
		ori_exc = (*use->exc_list)[use->exc_number];
		ori_exc->CopyTo(copy_exc);
  }

	// Find surface
  md->dl_type_x = NO_DL;
  if (use->sur_in)
  {
		if (use->sur_list->Count() <= i)
		{
			for (index = use->sur_list->Count(); index <= i; index++)
				use->sur_list->AddNew();
		}

		copy_sur = (*use->sur_list)[i];
		ori_sur = (*use->sur_list)[use->sur_number];
		ori_sur->CopyTo(copy_sur);
  }

	// Find gas
  if (use->gas_in)
  {
		if (use->gas_list->Count() <= i)
		{
			for (index = use->gas_list->Count(); index <= i; index++)
				use->gas_list->AddNew();
		}

		copy_gas = (*use->gas_list)[i];
		ori_gas = (*use->gas_list)[use->gas_number];
		ori_gas->CopyTo(copy_gas);
  }

	// Find solid solution
  if (use->ssa_in)
  {
		if (use->ssa_list->Count() <= i)
		{
			for (index = use->ssa_list->Count(); index <= i; index++)
				use->ssa_list->AddNew();
		}

		copy_ssa = (*use->ssa_list)[i];
		ori_ssa = (*use->ssa_list)[use->ssa_number];
		ori_ssa->CopyTo(copy_ssa);
  }

}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::SetReaction(int i)
{

	md->use.sol_p = md->use.sol_list->Element(i);
  if (md->use.sol_p == NULL)
		return false; //ToDo: must stop program. Severe error

	// Find pure phase assemblage
  if (md->use.ppa_in)
  {
    md->use.ppa_p = md->use.ppa_list->Element(i);
    if (md->use.ppa_p == NULL)
			return false; //ToDo: must stop program. Severe error
  }

	// Find exchange
  if (md->use.exc_in)
  {
    md->use.exc_p = md->use.exc_list->Element(i);
    if (md->use.exc_p == NULL)
			return false; //ToDo: must stop program. Severe error
  }

	// Find surface
  md->dl_type_x = NO_DL;
  if (md->use.sur_in)
  {
    md->use.sur_p = md->use.sur_list->Element(i);
    if (md->use.sur_p == NULL)
			return false; //ToDo: must stop program. Severe error
  }

	// Find gas
  if (md->use.gas_in)
  {
    md->use.gas_p = md->use.gas_list->Element(i);
    if (md->use.gas_p == NULL)
			return false; //ToDo: must stop program. Severe error
  }

	// Find s_s_assemblage
  if (md->use.ssa_in)
  {
    md->use.ssa_p = md->use.ssa_list->Element(i);
    if (md->use.ssa_p == NULL)
			return false; //ToDo: must stop program. Severe error
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::StepSaveSurf(Surface *dest)
{
	//
	// Save surface for intermediate calculation
	// Amt of surface may have changed due to reaction or surface related
	// to kinetic reactant.
	//
	// input:  n_user is user solution number of target
	//

	int i, j, k;
  Surface *surface_ptr;
	SurfaceComp *sc_p;
	ElementOfSpecies *eos_p;

/*
 *   Malloc space for solution data
 */
  if (md->use.sur_p == NULL)
    return;

  //surface_duplicate (md->use.sur_p->n_user, n_user);
	md->use.sur_p->CopyTo(dest);

  surface_ptr = dest;

	Master *master;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		master = gd->master_list[i];

    if (master->s->type != SURF)
      continue;

		for (j = 0; j < surface_ptr->comps->Count(); j++)
    {
			sc_p = (*surface_ptr->comps)[j];

			for (k = 0; k < sc_p->totals->Count(); k++)
      {
				eos_p = (*sc_p->totals)[k];
	
				if (eos_p->e == master->e)
				{
					if (master->total <= MIN_TOTAL)
						eos_p->coef = (LDBLE)MIN_TOTAL;
					else
						eos_p->coef = master->total;

					break;
				}
      }
    }
  }

  return;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::StepSaveExch(Exchange *dest)
{
	//
	// Save solution composition into structure solution with user number
	// n_user.
	//
	// input:  n_user is user solution number of target
	//
  int i, j, k;
  bool found;
  Exchange *exchange_ptr;
	ExchComp *ec_p;
	ElementOfSpecies *eos_p;

/*
 *   Malloc space for solution data
 */
  if (md->use.exc_p == NULL)
    return;

	md->use.exc_p->CopyTo(dest); //exchange_duplicate (use.exchange_ptr->n_user, n_user);
  exchange_ptr = dest; //exchange_bsearch (n_user, &n);

	Master *master;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		master = gd->master_list[i];

    if (master->s->type != EX)
      continue;

    found = false;

		for (j = 0; j < exchange_ptr->comps->Count(); j++)
    {
			ec_p = (*exchange_ptr->comps)[j];

			for (k = 0; k < ec_p->totals->Count(); k++)
      {
				eos_p = (*ec_p->totals)[k];

				if (eos_p->e == master->e)
				{
					if (!found)
					{
						found = true;

						if (master->total <= MIN_TOTAL)
							eos_p->coef = (LDBLE)MIN_TOTAL;
						else
							eos_p->coef = master->total;

						break;
					}
					else
						eos_p->coef = 0;
				}
      }
    }
  }

  return;
}
//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddPPAssemblage(PPAssemblage *pp_assemblage_ptr)
{
	/*
 *   Add a small amount of each phase if necessary to insure
 *   all elements exist in solution.
 */

  int i, j;
  LDBLE amount_to_add, total;
	char token[DEFAULT_STR_LENGTH];
  char *ptr;
  PurePhase *pure_phase_ptr;
  Master *master_ptr;
	ElementOfSpecies *eos_p;

  if (CheckPPAssemblage(pp_assemblage_ptr))
    return true;

#ifdef DEBUG_AddPPAssemblage
	d.addppassemblage_count++;
#endif
	/*
	 *   Go through list and generate list of elements and
	 *   coefficient of elements in reaction
	 */
  eos_list->Clear();
  parent_count = 0;

	/*
	 *   Check that all elements are in solution for phases with greater than zero mass
	 */
	for (j = 0; j < pp_assemblage_ptr->pure_phases->Count(); j++)
  {
		pure_phase_ptr = (*pp_assemblage_ptr->pure_phases)[j];

#ifdef DEBUG_AddPPAssemblage
		fprintf(d.addppassemblage_f, "%d 1 pure_phase_ptr->name: %s ", d.addppassemblage_count, pure_phase_ptr->name.CharPtr());
#endif

	  eos_list->Clear();
		parent_count = 0;

    amount_to_add = 0.0;
    pure_phase_ptr->delta = 0.0;

		if (!pure_phase_ptr->add_formula.IsEmpty())
    {
			pure_phase_ptr->add_formula.Copy(token);
      ptr = &(token[0]);
			GetElementsInSpecies(&ptr, 1.0);

#ifdef DEBUG_AddPPAssemblage
			fprintf(d.addppassemblage_f, "1\n");
#endif
    }
    else
    {
			CopyToTempEOSList(pure_phase_ptr->phase->eos_list, 1.0);

#ifdef DEBUG_AddPPAssemblage
			fprintf(d.addppassemblage_f, "2\n");
#endif
    }

    if (pure_phase_ptr->moles > 0.0)
    {
      for (i = 0; i < eos_list->Count(); i++)
      {
				eos_p = (*eos_list)[i];

				master_ptr = eos_p->e->primary;

				if (master_ptr->s == gd->s_hplus)
					continue;
				else if (master_ptr->s == gd->s_h2o)
					continue;
				else if (master_ptr->total > MIN_TOTAL)
					continue;
				else
				{
					total = (-master_ptr->total + (LDBLE)1e-10) / eos_p->coef;

					if (amount_to_add < total)
						amount_to_add = total;

#ifdef DEBUG_AddPPAssemblage
					fprintf(d.addppassemblage_f, "%d 2 i: %d eos_p->coef: %.20e master_ptr->total: %.20e amount_to_add: %.20e\n", d.addppassemblage_count, i, eos_p->coef, master_ptr->total, amount_to_add);
#endif
				}
      }

			if (pure_phase_ptr->moles < amount_to_add)
			{
				amount_to_add = pure_phase_ptr->moles;

#ifdef DEBUG_AddPPAssemblage
				fprintf(d.addppassemblage_f, "%d 3 amount_to_add: %.20e\n", d.addppassemblage_count, amount_to_add);
#endif
			}
    }

    if (amount_to_add > 0.0)
    {
      pure_phase_ptr->moles -= amount_to_add;
			pure_phase_ptr->delta = amount_to_add;

#ifdef DEBUG_AddPPAssemblage
			fprintf(d.addppassemblage_f, "%d 4 pure_phase_ptr->moles: %.20e pure_phase_ptr->delta: %.20e\n", d.addppassemblage_count, pure_phase_ptr->moles, pure_phase_ptr->delta);
#endif

			/*
			 *   Add reaction to totals
			 */
			for (i = 0; i < eos_list->Count(); i++)
      {
				eos_p = (*eos_list)[i];

				master_ptr = eos_p->e->primary;

				if (master_ptr->s == gd->s_hplus)
				{
					md->total_h_x += eos_p->coef * amount_to_add;

#ifdef DEBUG_AddPPAssemblage
					fprintf(d.addppassemblage_f, "%d 5 i: %d eos_p->coef: %.20e amount_to_add: %.20e md->total_h_x: %.20e\n", d.addppassemblage_count, i, eos_p->coef, amount_to_add, md->total_h_x);
#endif
				}
				else if (master_ptr->s == gd->s_h2o)
				{
					md->total_o_x += eos_p->coef * amount_to_add;

#ifdef DEBUG_AddPPAssemblage
					fprintf(d.addppassemblage_f, "%d 6 i: %d eos_p->coef: %.20e amount_to_add: %.20e md->total_o_x: %.20e\n", d.addppassemblage_count, i, eos_p->coef, amount_to_add, md->total_o_x);
#endif
				}
				else
				{
					master_ptr->total += eos_p->coef * amount_to_add;

#ifdef DEBUG_AddPPAssemblage
					fprintf(d.addppassemblage_f, "%d 7 i: %d eos_p->coef: %.20e amount_to_add: %.20e master_ptr->total: %.20e\n", d.addppassemblage_count, i, eos_p->coef, amount_to_add, master_ptr->total);
#endif
				}
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------------------------------------
bool ModelEngine::AddSSAssemblage(SSAssemblage *s_s_assemblage_ptr)
{
/*
 *   Accumulate solid_solution data in master->totals and _x variables.
 */
  int i, j, k;
  LDBLE amount_to_add, total;
  SS *s_s_ptr;
	SSComp *ssc_p;
	ElementOfSpecies *eos_p;
  Master *master_ptr;
	char token[DEFAULT_STR_LENGTH];
  char *ptr;

  if (s_s_assemblage_ptr == NULL)
    return true;

  eos_list->Count();
  parent_count = 0;

	/*
	 *   Check that all elements are in solution for phases with greater than zero mass
	 */
	for (i = 0; i < s_s_assemblage_ptr->ss_list->Count(); i++)
  {
		s_s_ptr = (*s_s_assemblage_ptr->ss_list)[i];

		eos_list->Count();
		parent_count = 0;

		for (j = 0; j < s_s_ptr->comps_list->Count(); j++)
    {
			ssc_p = (*s_s_ptr->comps_list)[j];

      amount_to_add = 0.0;
      ssc_p->delta = 0.0;

      if (ssc_p->moles > 0.0)
      {
				ssc_p->phase->formula.Copy(token);
				ptr = &(token[0]);
				GetElementsInSpecies(&ptr, 1.0);

				for (k = 0; k < eos_list->Count(); k++)
				{
					eos_p = (*eos_list)[k];

					master_ptr = eos_p->e->primary;

					if (master_ptr->s == gd->s_hplus)
						continue;
					else if (master_ptr->s == gd->s_h2o)
						continue;
					else if (master_ptr->total > MIN_TOTAL_SS)
						continue;
					else
					{
						total = (-master_ptr->total + (LDBLE)1e-10) / eos_p->coef;

						if (amount_to_add < total)
							amount_to_add = total;
					}
				}
      }

			if (ssc_p->moles < amount_to_add)
				amount_to_add = ssc_p->moles;

			if (amount_to_add > 0.0)
      {
				ssc_p->moles -= amount_to_add;
				ssc_p->delta = amount_to_add;

				/*
				 *   Add reaction to totals
				 */
				for (k = 0; k < eos_list->Count(); k++)
				{
					eos_p = (*eos_list)[k];

					master_ptr = eos_p->e->primary;

					if (master_ptr->s == gd->s_hplus)
						md->total_h_x += eos_p->coef * amount_to_add;
					else if (master_ptr->s == gd->s_h2o)
						md->total_o_x += eos_p->coef * amount_to_add;
					else
						master_ptr->total += eos_p->coef * amount_to_add;
				}
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::ResetData()
{
	count_unknowns = 0;

	//ToDo: Check to see if this is really necessary...
	ArrNewMax(0);
	ResidualNewMax(0);
	DeltaNewMax(0);

	/*
	arr_max = 0;
	residual_max = 0;
	delta_max = 0;
	*/
	/*
		delete [] arr;
		delete [] residual;
		delete [] delta;

		arr = NULL;
		residual = NULL;
		delta = NULL;

		arr_capacity = 0;
		residual_capacity = 0;
		delta_capacity = 0;
	*/
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SaveExchange(Exchange *dest)
{
	/*
	 *   Save exchanger assemblage into structure exchange with user 
	 *   number n_user.
	 */
  int i, j;
  //int count_comps;
  LDBLE charge;

  if (md->use.exc_p == NULL)
    return;

	md->use.exc_p->CopyTo(dest);

	dest->new_def = false;
	dest->solution_equilibria = false;

	dest->comps->Clear();

	Unknown *u;
	ExchComp *ec_p;
	Master *m;
	SpeciesInfo *s_i;
	for (i = 0; i < count_unknowns; i++)
	{
		u = (*unknown_list)[i];

    if (u->type == EXCH)
    {
			ec_p = dest->comps->AddNew();
			m = (*u->master)[0];

			u->exch_comp->CopyTo(ec_p);

      ec_p->master = m;
      ec_p->la = m->s->la;
			u->exch_comp->formula_totals->CopyTo(ec_p->formula_totals);
      ec_p->moles = 0.;

			/*
			 *   Save element concentrations on exchanger
			 */
      eos_list->Clear();
      parent_count = 0;
      charge = 0.0;
			for (j = 0; j < md->species_info_list->Count(); j++)
      {
				s_i = (*md->species_info_list)[j];

				if (s_i->master_s == m->s)
				{
					CopyToTempEOSList(s_i->s->eos_list, s_i->s->moles);
					charge += s_i->s->moles * s_i->s->z;
				}
			}

			/*
			 *   Keep exchanger related to phase even if none currently in solution
			 */
			if (!u->exch_comp->phase_name.IsEmpty() && eos_list->Count() == 0)
				CopyToTempEOSList(m->s->eos_list, (LDBLE)1e-20);

			/*
			 *   Store list
			 */
      ec_p->charge_balance = charge;
			SaveEOSList(ec_p->totals);

			/* update unknown pointer */
      u->exch_comp = ec_p;
    }
  }
	
  md->use.exc_p = NULL;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SaveGasPhase(GasPhase *dest)
{
/*
 *   Save gas composition into structure gas_phase with user 
 *   number n_user.
 */
  int i;
  //char token[MAX_LENGTH];

  if (md->use.gas_p == NULL)
    return;

/*
 *   Store in gas_phase
 */

	md->use.gas_p->CopyTo(dest);

  dest->new_def = false;
  dest->solution_equilibria = false;
  //temp_gas_phase.n_solution = -99;

	/*
	 *   Update amounts
	 */
	dest->comps->Clear();
	GasComp *gc_p, *gc_o_p;
	for (i = 0; i < md->use.gas_p->comps->Count(); i++)
  {
		gc_o_p = (*md->use.gas_p->comps)[i];
		gc_p = dest->comps->AddNew();

    gc_p->moles = gc_o_p->phase->moles_x;
  }

	/* update unknown pointer */
  if (md->gas_unknown != NULL)
    md->gas_unknown->gas_phase = dest;

  md->use.gas_p = NULL;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SaveSurface(Surface *dest)
{
	/*
	 *   Save surface data into structure surface with user 
	 *   number n_user.
	 */
  int i, j, last_charge;
  int count_charge;
  LDBLE charge;

  if (md->use.sur_p == NULL)
    return;

	/*
	 *   Store data for structure surface
	 */
	md->use.sur_p->CopyTo(dest);

  dest->new_def = false;
  dest->dl_type = md->dl_type_x;
  dest->solution_equilibria = false;

	/*
	 *   Allocate space to pointer comps
	 */

	dest->comps->Clear();
	dest->charge->Clear();

	/*
   *  Initial entry of surface sites is random
   *  Charge balance numbering follows the initial entry
   *  Surface sites are then sorted alphabetically
   *  Now when we save, the site order differs from the charge order
   *  last_charge sets up logic to renumber charge balance equations.
   */
  last_charge = -1;

	Unknown *x;
	SurfaceComp *sc_p;
	SurfaceCharge *sch_p;
	Master *m;
	SpeciesInfo *s;

	count_charge = 0;
	for (i = 0; i < count_unknowns; i++)
  {
		x = (*unknown_list)[i];
		m = (*x->master)[0];

    if (x->type == SURFACE)
    {
			sc_p = dest->comps->AddNew();

			x->surface_comp->CopyTo(sc_p);

      sc_p->master = m;
      sc_p->la = m->s->la;
      sc_p->moles = 0.;

      if (x->surface_comp->charge == last_charge)
				sc_p->charge = count_charge - 1;
      else
				sc_p->charge = count_charge;

      last_charge = x->surface_comp->charge;

			/*
			 *   Save element concentrations on surface
			 */
      eos_list->Clear();
      parent_count = 0;
      charge = 0.0;

			for (j = 0; j < md->species_info_list->Count(); j++)
      {
				s = (*md->species_info_list)[j];
	
				if (s->master_s == m->s)
				{
					CopyToTempEOSList(s->s->eos_list, s->s->moles);
					charge += s->s->moles * s->s->z;
				}
      }

      SaveEOSList(sc_p->totals);
      sc_p->totals->CopyTo(sc_p->formula_totals);
      sc_p->cb = charge;

      /* update unknown pointer */
      x->surface_comp = sc_p;
    }
    else if (x->type == SURFACE_CB && md->use.sur_p->type == DDL)
    {
			sch_p = dest->charge->AddNew();

			x->surface_charge->CopyTo(sch_p);
      sch_p->charge_balance = x->f;
      sch_p->mass_water =	x->surface_charge->mass_water;
      sch_p->diffuse_layer_totals->Clear();
      sch_p->g->Clear();

			/*
			 *   Added code to save g
			 */

			if (x->surface_charge->g->Count() > 0)
				x->surface_charge->g->CopyTo(sch_p->g);

      sch_p->la_psi = m->s->la;

			/*
			 *   Store moles from diffuse_layer
			 */
      if (md->dl_type_x != NO_DL)
      {
				SumDiffuseLayer(x->surface_charge);
				SaveEOSList(sch_p->diffuse_layer_totals);
      }

      /* update unknown pointer */
      x->surface_charge = sch_p;
      x->surface_comp = (*unknown_list)[i - 1]->surface_comp;
    }
    else if (x->type == SURFACE_CB && md->use.sur_p->type == CD_MUSIC)
    {
			sch_p = dest->charge->AddNew();

			x->surface_charge->CopyTo(sch_p);

      if (md->dl_type_x != NO_DL)
				sch_p->charge_balance = (x->surface_charge->sigma0 + x->surface_charge->sigma1 + x->surface_charge->sigma2 + x->surface_charge->sigmaddl) * (x->surface_charge->specific_area * x->surface_charge->grams) / F_C_MOL;
      else
				sch_p->charge_balance = (x->surface_charge->sigma0 + x->surface_charge->sigma1 + x->surface_charge->sigma2) * (x->surface_charge->specific_area * x->surface_charge->grams) / F_C_MOL;

			sch_p->mass_water = x->surface_charge->mass_water;
      sch_p->diffuse_layer_totals->Clear();
      sch_p->g->Clear();

			/*
			 *   Added code to save g
			 */

      if (x->surface_charge->g->Count() > 0)
				x->surface_charge->g->CopyTo(sch_p->g);

      sch_p->la_psi = m->s->la;

			/*
			 *   Store moles from diffuse_layer
			 */
      if (md->dl_type_x != NO_DL)
      {
				SumDiffuseLayer(x->surface_charge);
				SaveEOSList(sch_p->diffuse_layer_totals);
      }

      /* update unknown pointer */
      x->surface_charge = sch_p;
      x->surface_comp = (*unknown_list)[i - 1]->surface_comp;
    }
  }

  md->use.sur_p = NULL;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::FreeModelAllocs()
{
  //int i;

	//unknown_list->SetNewMax(0);
	count_unknowns = 0;

	//ToDo: Check to see if this is really necessary...
	ArrNewMax(0);
	ResidualNewMax(0);
	DeltaNewMax(0);

	/*
	arr_max = 0;
	residual_max = 0;
	delta_max = 0;
	*/
	/*
  delete [] arr;
	arr = NULL;
	arr_capacity = 0;

  delete [] delta;
	delta = NULL;
	delta_capacity = 0;

  delete [] residual;
	residual = NULL;
	residual_capacity = 0;
	*/

  md->s_x->Clear();

  md->sum_mb1->Clear();
  md->sum_mb2->Clear(); 
  md->sum_jacob0->Clear(); 
  md->sum_jacob1->Clear();
  md->sum_jacob2->Clear(); 
  md->sum_delta->Clear(); 

	gd->charge_group.Clear();
  
	return;
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SaveSolutionResults()
{
	for (int i = 0; i < count_unknowns; i++)
  {
    if ((*unknown_list)[i] == md->alkalinity_unknown)
    {
			(*unknown_list)[i]->molality_result = (*unknown_list)[i]->f / md->mass_water_aq_x;
			(*unknown_list)[i]->moles_result = (*unknown_list)[i]->f;
			continue;
    }
    
		if ((*unknown_list)[i] == md->ph_unknown)
      continue;
    
		if ((*unknown_list)[i] == md->pe_unknown)
      continue;
    
		if ((*unknown_list)[i] == md->charge_balance_unknown)
    {
			(*unknown_list)[i]->molality_result = (*unknown_list)[i]->sum / md->mass_water_aq_x;
			(*unknown_list)[i]->moles_result = (*unknown_list)[i]->sum;
      continue;
    }
    
		if ((*unknown_list)[i]->type == SOLUTION_PHASE_BOUNDARY)
    {
			(*unknown_list)[i]->molality_result = (*unknown_list)[i]->sum / md->mass_water_aq_x;
			(*unknown_list)[i]->moles_result = (*unknown_list)[i]->sum;
      continue;
    }
    else if ((*unknown_list)[i]->type == MB)
    {
			(*unknown_list)[i]->molality_result = (*unknown_list)[i]->sum / md->mass_water_aq_x;
			(*unknown_list)[i]->moles_result = (*unknown_list)[i]->sum;
      continue;
    }

		(*unknown_list)[i]->molality_result = 0.0;
		(*unknown_list)[i]->moles_result = 0.0;
  }
}
//-----------------------------------------------------------------------------------------------------------
void ModelEngine::SavePPResults(PPAssemblage *ppa)
{
	Unknown *x;
	PurePhase *pp;
	int i, j;

	for(i = 0; i < unknown_list->Count(); i++)
	{
		x = (*unknown_list)[i];

		if (x->type == PP)
		{
			pp = NULL;
			for(j = 0; j < ppa->pure_phases->Count(); j++)
				if (ppa->pure_phases->Element(j)->name == x->pure_phase->name)
				{
					pp = ppa->pure_phases->Element(j);
					break;
				}

			if (pp != NULL)
				pp->moles = x->moles;				
		}
	}
}

//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
