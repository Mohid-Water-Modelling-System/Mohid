#include "modeldata.h"
#include <float.h>
#include <math.h>

ModelData::ModelData(GlobalData *gd)
{
	this->gd = gd;

	LOG_10 = log((LDBLE)10.0);

	censor = 0.0;

	pe_x = NULL;

	min_value = (LDBLE)1e-10;
	aqueous_only = 0;

	alkalinity_unknown = NULL;
	carbon_unknown = NULL;
	ph_unknown = NULL;
	pe_unknown = NULL;
	charge_balance_unknown = NULL;
	solution_phase_boundary_unknown = NULL;
	mu_unknown = NULL;
	mb_unknown = NULL;
	ah2o_unknown = NULL;
	mass_hydrogen_unknown = NULL;
	mass_oxygen_unknown = NULL;
	pure_phase_unknown = NULL;
	gas_unknown = NULL;
	s_s_unknown = NULL;
	surface_unknown = NULL;

	s_x = NULL;

	mb_unknowns = NULL;

	sum_mb1 = NULL; 
	sum_mb2 = NULL;
	sum_jacob0 = NULL;
	sum_jacob1 = NULL;
	sum_jacob2 = NULL;
	sum_delta = NULL;
	
	species_info_list = NULL;

	negative_concentrations = false;

	pe_x = new List<PEData>;
	s_x = new ListOfPointers<Species>;
 	mb_unknowns = new List<UnknownInfo>;
	sum_mb1 = new List<STCoef>; 
	sum_mb2 = new List<STCoef>;
	sum_jacob0 = new List<STCoef>;
	sum_jacob1 = new List<STCoef>;
	sum_jacob2 = new List<STCoef>;
	sum_delta = new List<STCoef>;
	species_info_list = new SpeciesList;

	/*
	normal.SetPointerToArray(&normal_p);
	ineq_array.SetPointerToArray(&ineq_array_p);
	zero.SetPointerToArray(&zero_p);
	res.SetPointerToArray(&res_p);
	delta1.SetPointerToArray(&delta1_p);
	x_arg_e.SetPointerToArray(&x_arg);
	res_arg_e.SetPointerToArray(&res_arg);
	scratch_e.SetPointerToArray(&scratch);
	back_eq.SetPointerToArray(&back_eq_p);
	cu.SetPointerToArray(&cu_p);
	iu.SetPointerToArray(&iu_p);
	is.SetPointerToArray(&is_p);
	*/
	diagonal_scale = false;
	mass_water_switch = false;
	delay_mass_water = false;
	calculating_deriv = false;

	pe_step_size = 10.0;
	step_size = 100.0;
	pp_scale = 1.0;
	pp_column_scale = 1.0;

	ineq_tol = pow ((LDBLE) 10, (LDBLE) -DBL_DIG);
	convergence_tolerance = (LDBLE)1e-8; //If user wants more precision, swap to 1e-12 (Must put this as optional)

	same_model = false;
	same_temperature = false;

	itmax = 100;
}

ModelData::~ModelData()
{
	delete pe_x;
	delete s_x;
	delete mb_unknowns;
	delete sum_mb1;
	delete sum_mb2;
	delete sum_jacob0;
	delete sum_jacob1;
	delete sum_jacob2;
	delete sum_delta;
	delete species_info_list;
}

void ModelData::Reset()
{
	use.Reset();

	pe_x->Clear();
	s_x->Clear();
	mb_unknowns->Clear();
	sum_mb1->Clear(); 
	sum_mb2->Clear();
	sum_jacob0->Clear();
	sum_jacob1->Clear();
	sum_jacob2->Clear();
	sum_delta->Clear();

	species_info_list->Clear();

	alkalinity_unknown = NULL;
	carbon_unknown = NULL;
	ph_unknown = NULL;
	pe_unknown = NULL;
	charge_balance_unknown = NULL;
	solution_phase_boundary_unknown = NULL;
	mu_unknown = NULL;
	mb_unknown = NULL;
	ah2o_unknown = NULL;
	mass_hydrogen_unknown = NULL;
	mass_oxygen_unknown = NULL;
	pure_phase_unknown = NULL;
	gas_unknown = NULL;
	s_s_unknown = NULL;
	surface_unknown = NULL;
	exchange_unknown = NULL;

	//gd->x.Clear();
}

void ModelData::ResetXVariables()
{
  //new_x = FALSE;

  tc_x = 0.0;
  ph_x = 0.0;
  solution_pe_x = 0.0;
  mu_x = 0.0;
  ah2o_x = 0.0;
  density_x = 0.0;
  total_h_x = 0.0;
  total_o_x = 0.0;
  cb_x = 0.0;
  mass_water_aq_x = 0.0;

	max_unknowns = 0;

  //units_x = moles_per_kilogram_string;

	Master *m;
	for (int i = 0; i < gd->master_list.Count(); i++)
  {
		m = gd->master_list[i];

    m->total = 0.0;
    m->total_primary = 0.0;
    m->s->la = 0.0;
  }
}