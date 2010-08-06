#ifndef modeldataH
#define modeldataH

#include "use.h"
#include "specieslist.h"
#include "array.h"

class ModelData
{
public:
	ModelData(GlobalData *gd);
	~ModelData();

public:
	void Reset();
	void ResetXVariables();

public:
	Use use;

	List<class PEData> *pe_x;
	ListOfPointers<class Species> *s_x;
	int default_pe_x;

	LDBLE LOG_10,
				gfw_water,
				convergence_tolerance,
				min_value,
				pp_scale,
				pp_column_scale,
				ineq_tol,
				censor,
				xd,
				alpha,
				z;

	int max_unknowns;

	Unknown *alkalinity_unknown,
					*carbon_unknown,
					*ph_unknown,
					*pe_unknown,
					*charge_balance_unknown,
					*solution_phase_boundary_unknown,
					*mu_unknown,
					*mb_unknown,
					*ah2o_unknown,
					*mass_hydrogen_unknown,
					*mass_oxygen_unknown,
					*pure_phase_unknown,
					*gas_unknown,
					*s_s_unknown,
					*surface_unknown,
					*exchange_unknown;

	LDBLE mass_water_aq_x,
				mass_water_surfaces_x,
				mass_water_bulk_x,				
				mu_x,
				tk_x,
				tc_x,
				ph_x,
				solution_pe_x,
				ah2o_x,
				density_x,
				cb_x,
				total_ions_x,
				total_o_x,
				total_h_x,
				step_x;

	LDBLE	total_alkalinity,
				total_carbon,
				total_co2;
	
	int state,
			iterations,
			itmax,
			aqueous_only,
			g_iterations,
			count_total_steps,
			reaction_step;

	DIFFUSE_LAYER_TYPE dl_type_x; //Set on InitialSolution

	bool diagonal_scale,
			 initial_solution_isotopes,
			 incremental_reactions,
			 same_model,
			 same_temperature,
			 mass_water_switch,
			 delay_mass_water,
			 stop_program,
			 remove_unstable_phases,
			 numerical_deriv,
			 full_pitzer,
			 calculating_deriv,
			 negative_concentrations,
			 gas_in;

	SpeciesList *species_info_list;

	List<UnknownInfo> *mb_unknowns;
	List<STCoef> *sum_mb1, 
							 *sum_mb2,
							 *sum_jacob0,
							 *sum_jacob1,
							 *sum_jacob2,
							 *sum_delta;

	SurfaceCharge *surface_charge_ptr;

	LDBLE pe_step_size,
				pe_step_size_now,
				step_size,
				step_size_now;

	/*
	Array<LDBLE> normal,
							 ineq_array,
							 zero, 
							 res, 
							 delta1,
							 cu,
							 x_arg_e,
							 res_arg_e,
							 scratch_e;

	LDBLE *normal_p,
				*ineq_array_p,
			  *zero_p,
				*res_p,
				*delta1_p,
				*cu_p,
				*x_arg,
				*res_arg,
				*scratch;

	Array<int> back_eq,
						 iu,
						 is;

	int *back_eq_p,
			*iu_p,
			*is_p;
			*/

	LDBLE G_TOL;

public: //duvidosas
	LDBLE kin_time_x;
	int run_reactions_iterations;

private:
	GlobalData *gd;

};

#endif