#include "datastructures.h"

Solution::Solution(GlobalData *gd)
{
	pe = NULL;
	totals = NULL;
	ma = NULL;
	species_gamma = NULL;

	pe = new List<PEData>;
	totals = new List<Conc>;
	ma = new List<MasterActivity>;
	species_gamma = new List<MasterActivity>;

	this->gd = gd;

	Reset();
}

Solution::Solution()
{
	pe = NULL;
	totals = NULL;
	ma = NULL;
	species_gamma = NULL;

	pe = new List<PEData>;
	totals = new List<Conc>;
	ma = new List<MasterActivity>;
	species_gamma = new List<MasterActivity>;

	this->gd = NULL;;

	Reset();
}

Solution::Solution(Solution *copy)
{
	pe = NULL;
	totals = NULL;
	ma = NULL;
	species_gamma = NULL;

	pe = new List<PEData>;
	totals = new List<Conc>;
	ma = new List<MasterActivity>;
	species_gamma = new List<MasterActivity>;

	copy->CopyTo(this);
}

Solution::~Solution()
{
	delete pe;
	delete totals; 
	delete ma;
	delete species_gamma;
}

void Solution::Reset()
{
	tc = 25.0;
	ph = 7.0;
	solution_pe = 4.0;
	mu = (LDBLE)1e-7;
	ah2o = 1.0;
	density = 1.0;
	total_h = 0.0;
	total_o = 0.0;
	cb = 0.0;
	mass_water = 1.0;
	total_alk = 0.0;

	ph_conc_ptr = NULL;
	pe_conc_ptr = NULL;

	units = _mol_kgw_;

	default_pe = 0;

	number = -1;

	InitPE(gd);
	totals->Clear();
	ma->Clear();
	species_gamma->Clear();
}

void Solution::InitPE(GlobalData *gd)
{
	String e_minus("e-");

	if(gd == NULL)
		return;

	PEData *new_pe, *pe_0;

	pe->Clear();
	new_pe = pe->AddNew();

	if (gd->s_eminus != NULL && gd->s_eminus->rxn->token_list->Count() > 0)
		gd->s_eminus->rxn->CopyTo(new_pe->rxn);
	else
	{
		gd->s_eminus = gd->species_list.Store(&e_minus, false);
		gd->s_eminus->z = -1.0;

		pe_0 = (*pe)[0];

		ReactionToken *rt_ptr;

		rt_ptr = pe_0->rxn->token_list->AddNew();
		rt_ptr->s = gd->s_eminus;
		rt_ptr->coef = -1.0;

		rt_ptr = pe_0->rxn->token_list->AddNew();
		rt_ptr->s = gd->s_eminus;
		rt_ptr->coef = -1.0;
	}
}

Conc * Solution::SetConc(CONC_TYPE type, LDBLE value, String phase, LDBLE sat_index, bool charge)
{
	Conc *new_totals = totals->AddNew();

	if (type == _ph_)
		ph_conc_ptr = new_totals;
	else if (type == _pe_)
		pe_conc_ptr = new_totals;

	new_totals->type = type;
	new_totals->input_conc = value;
	new_totals->charge = charge;
	new_totals->phase_si = sat_index;
	if (!charge)
		new_totals->equation_name = phase;
	else
		new_totals->equation_name = "charge";

	return new_totals;
}

void Solution::SetPE(PXData &ped)
{
	if (ped.charge || !ped.phase.IsEmpty())
	{
		pe_conc_ptr = SetConc(_pe_, ped.value, ped.phase, ped.sat_index, ped.charge);

		solution_pe = pe_conc_ptr->input_conc;
		if (!pe_conc_ptr->equation_name.IsEmpty())
			pe_conc_ptr->name = "E";
	}
	else
		solution_pe = ped.value;
}

void Solution::SetPH(PXData &phd)
{
	if (phd.charge || !phd.phase.IsEmpty())
	{
		ph_conc_ptr = SetConc(_ph_, phd.value, phd.phase, phd.sat_index, phd.charge);

		ph = ph_conc_ptr->input_conc;
		if (!ph_conc_ptr->equation_name.IsEmpty())
			ph_conc_ptr->name = "H(1)";
	}
	else
		ph = phd.value;
}

int Solution::SetRedox(RedoxData &rd)
{
	PEData *new_pe, *pe_t;
	
	new_pe = new PEData;
	new_pe->element_1 = rd.element_1;
	new_pe->element_2 = rd.element_2;
	new_pe->valence_1 = rd.valence_1;
	new_pe->valence_2 = rd.valence_2;

	for (int index = 0; index < pe->Count(); index++)
		if (*pe->Element(index) == *new_pe)
		{
			pe_t = pe->Element(index);
			break;
		}

	if (pe_t == NULL)
		pe->AddNew(default_pe, new_pe); //will save a COPY of new_pe in the list

	delete new_pe;

	return default_pe;
}

int Solution::SetElementConc(ConcData &cd, int &MasterSpeciesID, LDBLE &gfw)
{
	int index;
	Conc *new_totals = totals->AddNew(index, &cd.element);

	new_totals->id = cd.id;
	new_totals->type = _conc_;
	new_totals->name = cd.element;
	new_totals->input_conc = cd.value;
	if (cd.use[0]) new_totals->as = cd.as;
	if (cd.use[1]) new_totals->units = cd.unit;
	if (cd.use[2]) new_totals->gfw = cd.gfw;
	if (cd.use[3]) new_totals->n_pe = SetRedox(cd.redox);
	if (cd.use[4])
	{
		new_totals->charge = cd.charge;

		if (!cd.charge)
		{
			new_totals->equation_name = cd.phase;
			new_totals->phase_si = cd.sat_index;
		}
	}	

	//Use the Master Species List ID for this Concentration ID between the interfaces
	Master *m = gd->master_list.Search(&new_totals->name, MasterSpeciesID, true);
	if (m == NULL) return -1;

	if (!m->primary) 
	{
		if (m->gfw_formula.Compare("0.0") == 0)
			gfw = 0.0;
		else
		{
			m = gd->master_list.Search(&m->gfw_formula, true);
			if (m == NULL) return -1;

			gfw = m->gfw;
		}
	}
	else gfw = m->gfw;

	return index;
}

void Solution::ChangeElementConc(char *element, LDBLE new_conc)
{
	Conc *c = totals->Search(element);

	if (c != NULL) 
		c->input_conc = new_conc;
}

void Solution::ChangeElementConc(int index, LDBLE new_conc)
{
	Conc *c = totals->Element(index);

	if (c != NULL) 
		c->input_conc = new_conc;
}

void Solution::CopyTo(Solution *dest)
{
	dest->tc = tc;
	dest->ph = ph;
	dest->solution_pe = solution_pe;
	dest->mu = mu;
	dest->ah2o = ah2o;
	dest->density = density;
	dest->total_h = total_h;
	dest->total_o = total_o;
	dest->cb = cb;
	dest->mass_water = mass_water;
	dest->total_alk = total_alk;

	dest->units = units;

	dest->default_pe = default_pe;

	dest->number = number;

	pe->CopyTo(dest->pe);
	totals->CopyTo(dest->totals);
	ma->CopyTo(dest->ma);
	species_gamma->CopyTo(dest->species_gamma);

	dest->gd = gd;
}

void Solution::Print(FILE *file)
{
	fprintf(file, "Solution info:\n\n");

	fprintf(file, "TC: %.20e\n", tc);
	fprintf(file, "pH: %.20e\n", ph);
	fprintf(file, "pE: %.20e\n", solution_pe);
	fprintf(file, "Mass of Water (kg): %.20e\n", mass_water);
	fprintf(file, "Density of Water: %.20e\n\n", density);

	fprintf(file, "Totals:\n");

	Conc *c;
	for (int i = 0; i < totals->Count(); i++)
	{
		c = totals->Element(i);

		fprintf(file, "(%d) Name: %s %s\n", i, c->name.CharPtr(), (c->charge?"(charge)":""));
		fprintf(file, "(%d) Type: %d\n", i, c->type);
		fprintf(file, "(%d) Input concentration: %.20e\n", i, c->input_conc);
		fprintf(file, "(%d) gfw: %.20e\n", i, c->gfw);
		fprintf(file, "(%d) As: %s\n", i, c->as);
		fprintf(file, "(%d) Equation: %s\n", i, c->equation_name.CharPtr());
		fprintf(file, "(%d) Phase SI: %.20e\n\n", i, c->phase_si);
	}
}

/*
Solution *Solution::Copy()
{
	return new Solution(this);
}
*/