#include "tidy.h"
#include "debug.h"

#include <math.h>
#include <conio.h>

//===========================
// Constructors & Destructors
//===========================
Tidy::Tidy(GlobalData *gd, ModelData *md):
Check(gd),
SSAssemblageModel(gd, md)
{
	rt_list  = NULL;
	rt_list  = new Reaction;

	this->gd = gd;
	this->md = md;
}

Tidy::~Tidy()
{
	if (rt_list != NULL) delete rt_list;
}

//===========================
// Public Functions
//===========================
bool Tidy::TidySpecies()
{
	Species *s;
	int i, j;

	if (!CheckSpecies())
		return false;

	for (i = 0; i < gd->species_list.Count(); i++)
  {
		s = gd->species_list[i];

    s->number = i;

    s->primary   = NULL;
    s->secondary = NULL;

		if (s->check_equation)
    {
      rt_list->AddSpeciesTRXN(s);

      if (!CheckEq(true))
				return false;
    }
  }

	Master *m;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];
		if (m->name[0] != '[' && m->name.FindUpper(1) != -1)
			return false;
    
    m->number = i;
    if (m->primary)
      m->s->primary = m;
    else
      m->s->secondary = m;

    if (m->name == "C")
      gd->s_co3 = m->s;

		if (!m->gfw_formula.IsEmpty() && !ComputeGFW(m->gfw_formula, m->gfw))
				return false;
  }
/*
 *   Write equations for all master species in terms of primary
 *   master species, set coefficient of element in master species
 */
  for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];

		rt_list->token_list->Clear();
    if (m->s->primary != NULL)
    {
			rt_list->AddTRXN(m->s->rxn, 1.0, false);
      rt_list->AddTRXN(m->s->rxn, -1.0, true);
    }
    else
    {
      rt_list->AddTRXN(m->s->rxn, 1.0, false);
      if (!RewriteEqToPrimary())
				return false;
    }

		m->rxn_primary->Reset();
		rt_list->CopyReactions(m->rxn_primary);

    if (!CoefInMaster(m, m->coef))
			return false;
  }

	for (i = 0; i < gd->species_list.Count(); i++)
  {
		if (i == 115)
			printf("");

    s = gd->species_list[i];

		rt_list->token_list->Clear();

    if (s->primary != NULL || s->secondary != NULL)
    {
      rt_list->AddTRXN(s->rxn, 1.0, false);
      rt_list->AddTRXN(s->rxn, -1.0, true);
    }
    else
    {
      rt_list->AddTRXN(s->rxn, 1.0, false);

			if (!RewriteEqToSecondary())
				return false;
    }

		s->rxn_s->Reset();
		rt_list->CopyReactions(s->rxn_s);

    if (!CalcAlk(s->rxn_s, s->alk))
			return false;

		s->co2 = 0.0;
    for (j = 1; j < rt_list->token_list->Count(); j++)
    {
			if (rt_list->token_list->Element(j)->s == gd->s_co3)
      {
				s->co2 = rt_list->token_list->Element(j)->coef;
				break;
      }
    }
  }

	Element *e;
	for (i = 0; i < gd->element_list.Count(); i++)
  {
		e = gd->element_list[i];

		e->master = gd->master_list.Search(&e->name, true); //master_bsearch (elements[i]->name);
    if (e->master == NULL)
			return false;

		e->primary = SearchForMasterPrimary(e->name);
    if (e->primary == NULL)
			return false;
  }

	Master *m_primary;
	for (i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];

    if (!m->primary)
    {
      m_primary = m->s->secondary->e->primary;
      
			if (m_primary == NULL || m_primary->s->secondary == NULL)
				return false;
    }
  }

	ElementOfSpecies *eos;
	for (i = 0; i < gd->species_list.Count(); i++)
  {
		s = gd->species_list[i];

    if (s->e_sec_list->Count() > 0)
    {
      s->h = 0.0;
      s->o = 0.0;

			for (j = 0; j < s->e_sec_list->Count(); j++)
      {
				eos = (*s->e_sec_list)[j];

				if (eos->e->primary == NULL)
					continue;

				if (eos->e->primary->s == gd->s_hplus) 
					s->h += eos->coef;
				else if (eos->e->primary->s == gd->s_h2o) 
					s->o += eos->coef;
				else if (!s->mole_balance.IsEmpty())
				{
					m = eos->e->master;

					if (m->primary && m->s->secondary != NULL)
						m = m->s->secondary;

					if (m->coef != 1)
						eos->coef /= m->coef;
				}
      }
      if (s->type == EX)
      {
				for (j = 0; j < s->e_sec_list->Count(); j++)
		    {
					eos = (*s->e_sec_list)[j];

					if (eos->e->primary->s->type == EX)
					{
						s->equiv = eos->coef;
						break;
					}
				}
      }
    }
  }

	return true;
}

bool Tidy::TidyPhases()
{
  int i;
  bool replaced;

	Check cd(gd);

	eos_list->Clear();

	Phase *p;
	for (i = 0; i < gd->phase_list.Count(); i++)
  {
		p = gd->phase_list[i];

		cd.SelectLOGKExpression(p->logk, p->rxn->logk);

		if (p->rxn->token_list->Count() > 0)
		{
			p->rxn->token_list->Element(0)->name = p->name;    
			p->rxn->token_list->Element(0)->s = NULL;
		}
  }

  for (i = 0; i < gd->phase_list.Count(); i++)
  {
		rt_list->token_list->Clear();
		p = gd->phase_list[i];

    rt_list->AddPhaseTRXN(p->rxn, 1.0, false);

		if (rt_list->token_list->Count() > 0)
			rt_list->token_list->Element(0)->name = p->name;

		if (!ReplaceSolidsGases(replaced))
			return false;

    p->rxn->Reset();
    rt_list->CopyReactions(p->rxn);

		rt_list->TRXNReverseK();
		RewriteEqToSecondary();
		rt_list->TRXNReverseK();

		p->rxn_s->Reset();
    rt_list->CopyReactions(p->rxn_s);

		if (p->check_equation)
    {
      if (!replaced)
				rt_list->AddPhaseTRXN(p, p->rxn);
      else
				rt_list->AddPhaseTRXN(p, p->rxn_s);

			if (!CheckEq(false))
				return false;
    }
  }

	return true;
}

bool Tidy::TidySolution(Solution *sol)
{
	String token;
	sol->totals->Sort(0, false);

	Conc *c;
	for (int i = 0; i < sol->totals->Count(); i++)
	{
		c = (*sol->totals)[i];

		if (c->name == "H(1)" || c->name == "E")
		{
			c->moles = 0.0;
			continue;
		}

		if (gd->master_list.Search(&c->name, true) == NULL)
			return false;

		if (c->n_pe < 0)
			c->n_pe = sol->default_pe;

    token = c->name;
		token.ToLower();

		if (c->units == _unknown_unit_)
      c->units = sol->units;
    else if (token == "alk")
		{
			if (!CheckUnitForAlk(c->units, sol->units))
				return false;
		}
		else
		{
			if (!CheckUnit(c->units, sol->units))
				return false;
    }
	}

	return true;
}

bool Tidy::TidyPurePhase(PPAssemblage *ppa)
{
  int j, l, first;
  Phase *p;
	PurePhase *pp;
	ElementOfSpecies *eos_p;
  LDBLE coef;
  char *ptr, token[DEFAULT_STRING_LENGTH];

	eos_list->Clear();
	parent_count = 0;

	coef = 1.0;

	for (j = 0; j < ppa->pure_phases->Count(); j++)
  {
		pp = (*ppa->pure_phases)[j];

		p = gd->phase_list.Search(&pp->name, true);
    
		if (p == NULL)
			return false;
    else
    {
			pp->phase = p;
			CopyToTempEOSList(p->eos_list, coef);
    }
		
		if (!pp->add_formula.IsEmpty())
    {
			first = eos_list->Count();
			
			p = gd->phase_list.Search(&pp->add_formula, true);

			if (p != NULL)
				pp->add_formula = p->formula;

			pp->add_formula.Copy(token);
			ptr = &token[0];

			GetElementsInSpecies(&ptr, coef);

			for (l = first; l < eos_list->Count(); l++)
			{
				eos_p = (*eos_list)[l];

				if (eos_p->e->master == NULL)
					return false;
			}
    }

		SaveEOSList(ppa->eos_list);
	}

  return true;
}

bool Tidy::CheckObrigatorySpecies()
{
	String h_one("H(1)");

	if (gd->s_hplus == NULL && gd->s_h3oplus != NULL)
  {
    gd->s_hplus = gd->s_h3oplus;
    gd->s_h3oplus = NULL;
  }

  if (gd->s_h2o == NULL || gd->s_h2o->primary == NULL ||	gd->s_h2o->secondary == NULL)
		return false;

	if ((gd->s_hplus == NULL && gd->s_h3oplus == NULL) || (gd->s_hplus != NULL && gd->s_h3oplus != NULL))
		return false;

	if (gd->s_hplus->primary == NULL || gd->s_hplus->secondary == NULL)
		return false;

	if (gd->s_eminus == NULL || gd->s_eminus->primary == NULL)
		return false;

  if (gd->s_h2 == NULL || gd->s_o2 == NULL)
		return false;

	gd->e_h_one = gd->element_list.Store(&h_one, false);

	if (gd->e_h_one == NULL)
		return false;

	return true;
}

//===========================
// Protected Functions
//===========================
bool Tidy::CheckEq(bool association)
{
  int i;
	LDBLE coef;
	String name_of_element;
	char t_ptr[DEFAULT_STRING_LENGTH], *ptr;
  LDBLE sumcharge;

  parent_count = 0;
	eos_list->Clear();

	if (!Equal(rt_list->token_list->Element(0)->coef, -1.0, TOLERANCE))
    return false;

  sumcharge = 0.0;
	ReactionToken *token;
	for (i = 0; i < rt_list->token_list->Count(); i++)
  {
		token = rt_list->token_list->Element(i);

    sumcharge += token->coef * token->z;
		token->name.Copy(t_ptr);
		ptr = &t_ptr[0];
    if (!GetElementsInSpecies(&ptr, token->coef))
      return false;
  }

	eos_list->Sort();
	if (!CombineElements())
    return false;

	if (!Equal(sumcharge, 0.0, TOLERANCE))
		return false;


	for (i = 0; i < eos_list->Count(); i++)
	{
		coef = (*eos_list)[i]->coef;
		name_of_element = (*eos_list)[i]->e->name;

		if (!Equal(coef, 0.0, TOLERANCE) && name_of_element != "e")
			return false;
	}

	return true;
}

bool Tidy::RewriteEqToPrimary()
{
  int j, add_count;
	bool repeat;

  repeat = true;
  add_count = 0;

	while (repeat)
  {
    repeat = false;

    if (++add_count > MAX_ADD_EQUATIONS)
			return false;

		for (j = 1; j < rt_list->token_list->Count(); j++)
    {
      if (rt_list->token_list->Element(j)->s->primary == NULL)
      {
				rt_list->AddTRXN(rt_list->token_list->Element(j)->s->rxn, rt_list->token_list->Element(j)->coef, true);
				repeat = true;
				break;
      }
    }
  }

  rt_list->TRXNCombine();

	return true;
}

bool Tidy::CoefInMaster(Master *m, LDBLE &coef)
{
	String element;
	char elt_name[DEFAULT_STRING_LENGTH], *ptr;;

  coef = 0.0;
	m->name.Copy(elt_name);
	ptr = &elt_name[0];
  if (!GetElement(&ptr, element))
		return false;

	for (int i = 0; i < m->s->eos_list->Count(); i++)
  {
    if (element == (*m->s->eos_list)[i]->e->name)
    {
      coef = (*m->s->eos_list)[i]->coef;
      break;
    }
  }

	return true;
}

bool Tidy::RewriteEqToSecondary()
{
  bool repeat;
	int i, add_count;
  ReactionToken *token;

	add_count = 0;
  repeat = true;

	while (repeat)
  {
    repeat = false;

		if (++add_count > MAX_ADD_EQUATIONS)
			return false;

    for (i = 1; i < rt_list->token_list->Count(); i++)
    {
			token = rt_list->token_list->Element(i);

      if (token->s == NULL)
				return false;

			if (token->s->secondary == NULL && token->s->primary == NULL)
      {
				rt_list->AddTRXN(token->s->rxn, token->coef, true);
				repeat = true;
				break;
      }
    }
  }

	rt_list->TRXNCombine();
  return true;
}

bool Tidy::CalcAlk(Reaction *reaction, LDBLE &alk)
{
  int i;
  Master *m;

  alk = 0.0;

	for (i = 1; i < reaction->token_list->Count(); i++)
  {
    m = reaction->token_list->Element(i)->s->secondary;

    if (m == NULL)
      m = reaction->token_list->Element(i)->s->primary;

    if (m == NULL)
			return false;

    alk += reaction->token_list->Element(i)->coef * m->alk;
  }

  return true;
}

Master * Tidy::SearchForMasterPrimary(String &name)
{
	char token[DEFAULT_STRING_LENGTH], *ptr;
	String element;

	name.Copy(token);
  ptr = &token[0];
  GetElement(&ptr, element);

	return gd->master_list.Search(&element, true);
}

bool Tidy::ReplaceSolidsGases(bool &replace)
{
  LDBLE coef;
  int i, add_count;
	bool repeat;
  ReactionToken *rt;
  Phase *phase_ptr;
  String token;

  add_count = 0;
  repeat = true;
  replace = false;

	while (repeat == true)
  {
    repeat = false;

		if (++add_count > MAX_ADD_EQUATIONS)
			return false;

    for (i = 1; i < rt_list->token_list->Count(); i++)
    {
      rt = rt_list->token_list->Element(i);

      if (rt->s == NULL)
      {
				phase_ptr = gd->phase_list.Search(&rt->name);

				if (phase_ptr == NULL)
				{
					token = rt->name;
					token.Replace ("(g)", "");
					token.Replace ("(s)", "");
					token.Replace ("(G)", "");
					token.Replace ("(S)", "");
					phase_ptr = gd->phase_list.Search(&token);
				}
				
				if (phase_ptr == NULL)
					return false;

				coef = rt->coef;

				rt_list->AddPhaseTRXN(phase_ptr->rxn, coef, false);

				rt->name = phase_ptr->rxn->token_list->Element(0)->name;
				rt->s = phase_ptr->rxn->token_list->Element(0)->s;
				rt->coef = -coef * phase_ptr->rxn->token_list->Element(0)->coef;
	
				repeat = true;		
				replace = true;

				rt_list->TRXNCombine();
				break;
      }
    }
  }

  rt_list->TRXNCombine ();

  return true;
}

//-----------------------------------------------------------------------------------------------------------
bool Tidy::TidyGasPhase(GasPhase *gas)
{
  int j;

	GasComp *gc_p;
  for (j = 0; j < gas->comps->Count(); j++)
  {
		gc_p = (*gas->comps)[j];

		gc_p->phase = gd->phase_list.Search(&gc_p->name, true);
    if (gc_p->phase == NULL)
			return false;
		
    if (gas->type == PRESSURE) // Fixed pressure
    {			
			if (gas->solution_equilibria)
				return false;
				
			// calculate moles
			if (gc_p->p_read != NAN)
			  gc_p->moles = gc_p->p_read * gas->volume / R_LITER_ATM / gas->temperature;
			else
			  return false;
    }
    else // Fixed volume
    {			
			if (!gas->solution_equilibria)
			{
			  if (gc_p->p_read != NAN)
			    gc_p->moles = gc_p->p_read * gas->volume / R_LITER_ATM / gas->temperature;
			  else
					return false;
      }
    }
  }

	if (gas->solution_equilibria)
		gas->new_def = true;
	else
		gas->new_def = false;
	
  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool Tidy::TidySurface(Surface *sur)
{
	//
	// After all of data are read, fill in master species for surface comps Sort surface
	//

  int i, j;
	char *ptr1, token[DEFAULT_STR_LENGTH];
	LDBLE specific_area, grams;
	SurfaceComp *sc_p;
	Master *m_p;
	ElementOfSpecies *eos_p;
	SurfaceCharge *sch_p;

	sur_p = sur;

	for (i = 0; i < sur_p->comps->Count(); i++)
  {
		sc_p = (*sur_p->comps)[i];

		// Find master species for each surface
		for (j = 0; j < sc_p->totals->Count(); j++)
    {
			eos_p = (*sc_p->totals)[j];

			if ((m_p = eos_p->e->master) == NULL)
				return false;

			if (m_p->type != SURF)
				continue;

			// Set flags
			sc_p->master = m_p;

			// Calculate moles of sites
			if (sur_p->sites_units == SITES_DENSITY && sc_p->phase_name.IsEmpty())
			{
				specific_area = (*sur_p->charge)[sc_p->charge]->specific_area;
				grams = (*sur_p->charge)[sc_p->charge]->grams;

				sc_p->moles = sc_p->moles * (LDBLE)1.0e18 * specific_area * grams / AVOGADRO;

				// Calculate totals
				eos_list->Clear();
				parent_count = 0;

				sc_p->formula.Copy(token);
				ptr1 = &token[0];
				GetElementsInSpecies(&ptr1, sc_p->moles);

				SaveEOSList(sc_p->formula_totals);
				SaveEOSList(sc_p->totals);
			}
	
			if (sur_p->type == CD_MUSIC)
			{
				sch_p = (*sur_p->charge)[sc_p->charge];				
				sch_p->charge_balance += sc_p->moles * sc_p->formula_z;
			}
	
			break;
    }
  }

	// Sort components
	sur_p->comps->Sort(0, true);

	return true;
}
//-----------------------------------------------------------------------------------------------------------
bool Tidy::TidyMinExchange(Exchange *exc)
{
	//
	// If exchanger is related to mineral, exchanger amount is set in proportion
	//

	int j, k, jj;
	char *ptr, token[DEFAULT_STR_LENGTH];
  LDBLE conc;
	ExchComp *excc_p;
	ElementOfSpecies *eos_p;
	Master *m_p;
	PurePhase *pp_p;

	exc_p = exc;
	ppa_p = md->use.ppa_p;

	if (exc_p == NULL)
		return false;

	for (j = 0; j < exc_p->comps->Count(); j++)
  {
		excc_p = (*exc_p->comps)[j];

		if (excc_p->phase_name.IsEmpty())
			continue;

    excc_p->master = NULL;

    // First find exchange master species

    for (k = 0; k < excc_p->totals->Count(); k++)
    {
			eos_p = (*excc_p->totals)[k];

			// Find master species
			if ((m_p = eos_p->e->master) == NULL)
				return false;

			if (m_p->type != EX)
				continue;

			excc_p->master = m_p;
			break;
    }
    
		if (excc_p->master == NULL)
			return false;

    // Now find the mineral on which exchanger depends...  
		if (ppa_p == NULL)
			return false;

		for (k = 0; k < ppa_p->pure_phases->Count(); k++)
    {
			pp_p = (*ppa_p->pure_phases)[k];

			if (excc_p->phase_name.Compare(pp_p->name, true) == 0)
			  break;
    }

		if (k == ppa_p->pure_phases->Count())
			return false;

    // use database name for phase
    excc_p->phase_name = pp_p->phase->name;

    // make exchanger concentration proportional to mineral ...
    conc = pp_p->moles * excc_p->phase_proportion;

		eos_list->Clear();
    parent_count = 0;

		excc_p->formula.Copy(token);
		ptr = &token[0];

    if (!GetElementsInSpecies(&ptr, conc))
			return false;

		excc_p->totals->Clear();
    SaveEOSList(excc_p->totals);

		// make sure exchange elements are in phase
		eos_list->Clear();
    parent_count = 0;

		excc_p->formula.Copy(token);
		ptr = &token[0];

		if (!GetElementsInSpecies(&ptr, -excc_p->phase_proportion))
			return false;
      
		pp_p->phase->formula.Copy(token);
		ptr = &token[0];
		if (!GetElementsInSpecies(&ptr, 1.0))
			return false;

		if(!CombineElements())
			return false;

		for (jj = 0; jj < eos_list->Count(); jj++)
    {
			eos_p = (*eos_list)[jj];

			if (eos_p->e->primary->s->type != EX && eos_p->coef < 0)
				return false;
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool Tidy::TidyMinSurface(Surface *sur)
{
	//
	// If surface is related to mineral, surface amount is set in proportion
	//

  int j, k, jj;
  LDBLE conc;
	char *ptr, token[DEFAULT_STR_LENGTH];

	sur_p = sur;
	ppa_p = md->use.ppa_p;

	SurfaceComp *sc_p, *sc_p2;
	Master *m_p;
	ElementOfSpecies *eos_p, *eos_p2;
	PurePhase *pp_p;
	SurfaceCharge *sch_p;
  for (j = 0; j < sur_p->comps->Count(); j++)
  {
		sc_p = (*sur_p->comps)[j];
		sch_p = (*sur_p->charge)[sc_p->charge];

		if (sc_p->phase_name.IsEmpty())
			continue;

    sc_p->master = NULL;

    // First find surface master species

		for (k = 0; k < sc_p->totals->Count(); k++)
    {
			eos_p = (*sc_p->totals)[k];

			// Find master species 				
			m_p = eos_p->e->master;

			if (m_p == NULL)
				return false;

			if (m_p->type != SURF)
				continue;

			sc_p->master = m_p;
			break;
    }
    
		if (sc_p->master == NULL)
			return false;

    // Now find the mineral on which surface depends...		
    if (ppa_p == NULL)
			return false;

		for (k = 0; k < ppa_p->pure_phases->Count(); k++)
    {
			pp_p = (*ppa_p->pure_phases)[k];

			if (sc_p->phase_name.Compare(pp_p->name, true) == 0)
			  break;
    }

		if (k == ppa_p->pure_phases->Count())
			return false;

		if (pp_p->phase == NULL)
			return false;

    // use database name for phase
    sc_p->phase_name = pp_p->phase->name;

    // make surface concentration proportional to mineral ...
    conc = pp_p->moles * sc_p->phase_proportion;

		eos_list->Clear();
    parent_count = 0;

		sc_p->formula.Copy(token);
		ptr = &token[0];

    if (!GetElementsInSpecies(&ptr, conc))
			return false;

		sc_p->totals->Clear();
		SaveEOSList(sc_p->totals);

    // area
		sch_p->grams = pp_p->moles;

		// make sure surface elements are in phase
		// logically necessary for mass balance and to avoid negative concentrations when dissolving phase
		eos_list->Clear();
    parent_count = 0;

		pp_p->phase->formula.Copy(token);
		ptr = &token[0];

		if (!GetElementsInSpecies(&ptr, 1.0))
			return false;
      
		for (jj = 0; jj < sur_p->comps->Count(); jj++)
    {
			sc_p2 = (*sur_p->comps)[jj];

			if (sc_p->charge != sc_p2->charge)
				continue;

			if (sur_p->type == CD_MUSIC)
			{
				sc_p2->formula.Copy(token);
				ptr = &token[0];

				if (!GetElementsInSpecies(&ptr, -sc_p2->phase_proportion))
					return false;
			}
			else
			{
				if (sc_p2->master->s->z != 0.0)
					return false;

				sc_p2->master->s->name.Copy(token);
				ptr = &token[0];

				if (!GetElementsInSpecies(&ptr, -sc_p2->phase_proportion))
					return false;
			}
    }

		if (!CombineElements())
			return false;

    for (jj = 0; jj < eos_list->Count(); jj++)
    {
			eos_p2 = (*eos_list)[jj];

			if (eos_p2->e->primary->s->type != SURF && 
				  eos_p2->coef < 0 && 
					eos_p2->e->primary->s != gd->s_hplus && 
					eos_p2->e->primary->s != gd->s_h2o)
				return false;
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
bool Tidy::TidySSAssemblage(SSAssemblage *ssa)
{
	int j, k;
  LDBLE nb, nc, n_tot, xb, xc, dnb, dnc, a0, a1;
  LDBLE xb2, xb3, xb4, xc2, xc3;
  LDBLE moles;
	
	if (ssa == NULL)
		return true;

	SS *ss_p;
	SSComp *ssc_p, *ssc_p0, *ssc_p1;
  for (j = 0; j < ssa->ss_list->Count(); j++)
  {
		ss_p = (*ssa->ss_list)[j];

    for (k = 0; k < ss_p->comps_list->Count(); k++)
    {
			ssc_p = (*ss_p->comps_list)[k];

			ssc_p->phase = gd->phase_list.Search(&ssc_p->name, true);

			if (ssc_p->phase == NULL)
				return false;
			
			ssc_p->phase->moles_x = 0;
			ssc_p->phase->fraction_x = 0;
		}
				
		if (ssc_p->moles == NAN)
			return false;

		if (!SSCalcA0A1(ss_p))
			return false;

		n_tot = 0;
		for (k = 0; k < ss_p->comps_list->Count(); k++)
		{
			ssc_p = (*ss_p->comps_list)[k];

			moles = ssc_p->moles;

			if (ssc_p->moles <= 0.0)
			{
				moles = (LDBLE)MIN_TOTAL_SS;
				ssc_p->initial_moles = moles;
			}

			n_tot += moles;
		}

		for (k = 0; k < ss_p->comps_list->Count(); k++)
		{
			ssc_p = (*ss_p->comps_list)[k];

			moles = ssc_p->moles;

			if (ssc_p->moles <= 0.0)
			{
				moles = (LDBLE)MIN_TOTAL_SS;
		  }

			ssc_p->fraction_x = moles / n_tot;
			ssc_p->log10_fraction_x = log10 (moles / n_tot);
		}
				
		a0 = ss_p->a0;
		a1 = ss_p->a1;

		if (a0 != 0.0 || a1 != 0)
		{
			ssc_p0 = (*ss_p->comps_list)[0];
			ssc_p1 = (*ss_p->comps_list)[1];

		  ss_p->dn = (LDBLE)1.0 / n_tot;

			nc = ssc_p0->moles;					
			if (nc == 0)
				nc = (LDBLE)MIN_TOTAL_SS;
						
			nb = ssc_p1->moles;					
			if (nb == 0)
			  nb = (LDBLE)MIN_TOTAL_SS;
						
			xc = nc / n_tot;
			xb = nb / n_tot;

			ssc_p0->log10_lambda = xb * xb * (a0 - a1 * (3 - 4 * xb)) / md->LOG_10;
			ssc_p1->log10_lambda = xc * xc * (a0 + a1 * (4 * xb - 1)) / md->LOG_10;

			xc2 = xc * xc;
			xc3 = xc2 * xc;
			xb2 = xb * xb;
			xb3 = xb2 * xb;
			xb4 = xb3 * xb;

			dnb = -2 * a0 * xb * xc2 - 8 * a1 * xb2 * xc2 + 6 * a1 * xb * xc2 - 4 * a1 * xc * xb4 - 8 * a1 * xb3 * xc2 - 4 * a1 * xb2 * xc3 - 2 * a0 * xc * xb2 - 8 * a1 * xc * xb3 + 6 * a1 * xc * xb2 + 1;
			ssc_p0->dnb = dnb / n_tot;
				  
			dnc = 2 * a0 * xb3 + 2 * a0 * xc * xb2 + 8 * a1 * xb4 + 8 * a1 * xc * xb3 - 2 * a1 * xb3 - 6 * a1 * xc * xb2;				  
			ssc_p0->dnc = -xb / nc + dnc / n_tot;
			ssc_p0->dn = (LDBLE)1.0 / n_tot;

			dnb = 2 * a0 * xb * xc2 + 2 * a0 * xc3 + 8 * a1 * xb2 * xc2 + 8 * a1 * xb * xc3 - 2 * a1 * xb * xc2 - 6 * a1 * xc3;
			ssc_p1->dnb = -xc / nb + dnb / n_tot;
			
			dnc = -2 * a0 * xc * xb2 - 8 * a1 * xc * xb3 + 2 * a1 * xc * xb2 - 2 * a0 * xb * xc2 - 8 * a1 * xb2 * xc2 + 6 * a1 * xb * xc2 + 1;
			ssc_p1->dnc = dnc / n_tot;
				  
			if (!SSPrep(ss_p->tk, ss_p))
				return false;
			
			ssc_p1->dn = (LDBLE)1.0 / n_tot;
		}
		else
		{
			ss_p->dn = (LDBLE)1.0 / n_tot;
					
			for (k = 0; k < ss_p->comps_list->Count(); k++)
			{
				ssc_p = (*ss_p->comps_list)[k];

				ssc_p->log10_lambda = 0;
				
				moles = ssc_p->moles;
					  
				if (moles <= 0.0)
					moles = (LDBLE)MIN_TOTAL_SS;
					    
				ssc_p->dnb = (n_tot - moles) / (moles * n_tot);
				ssc_p->dn = (LDBLE)1.0 / n_tot;
			}
		}
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
