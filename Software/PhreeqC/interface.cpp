#include "model.h"

#define IntToBool(a) (a == 1?true:false)

#define INT_CONCENTRATION  1
#define INT_PHASE          2
#define INT_GAS            3
#define INT_SURFACE        4
#define INT_SPECIES        5
#define INT_EXCHANGE       6
#define INT_PUREPHASESOLID 7
#define INT_PUREPHASEGAS   8

#define INT_PHASE_SOLID 1
#define INT_PHASE_GAS   2

ConcData cd;
PPData pd;
PXData px;
ExchCompData ed;

List<Model> instances;
Model *model;
Solution *sol;

void pm_start(int *instance_id, int *status)
{
	*status = 0;

	try
	{
		model = instances.AddNew();

		if (model == NULL)
			return;

		*status      = 1;
		*instance_id = instances.Count() - 1;

	}
	catch(...)
	{
		*instance_id = -1;
	}
}

void pm_model_option_printalways(int *instance_id, int *print_always, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->SetOptions(0, (*print_always == 1?true:false));

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_data(int *instance_id, LDBLE *temperature, int *units, LDBLE *density, LDBLE *mass_water, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		sol->tc         = *temperature;
		sol->units      = *units;
		sol->mass_water = *mass_water;
		sol->density    = *density;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_redox(int *instance_id, char *element_1, LDBLE *valence_1, char *element_2, LDBLE *valence_2, int *status, int c1, int c2)
{
	*status = 1;

	//ToDo: Must create the code to do this on the solution object (See read_solution on the original PhreeqC read.c file)
}

void pm_solution_ph(int *instance_id, LDBLE *value, int *charge, int *use_phase, char *phase_name, LDBLE *si, int *status, int c)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		px.value = *value;

		if (*charge == 1)
			px.charge = true;
		else 
			px.charge = false;

		if (*use_phase == 1)
		{
			px.phase.SetString(phase_name);
			px.sat_index = *si;
		}
		else
		{
			px.phase.SetString("");
			px.sat_index = 0.0;
		}

		sol->SetPH(px);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_pe(int *instance_id, LDBLE *value, int *charge, int *use_phase, char *phase_name, LDBLE *si, int *status, int c)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		px.value = *value;

		if (*charge == 1)
			px.charge = true;
		else 
			px.charge = false;

		if (*use_phase == 1)
		{
			px.phase.SetString(phase_name);
			px.sat_index = *si;
		}
		else
		{
			px.phase.SetString("");
			px.sat_index = 0.0;
		}

		sol->SetPE(px);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_add(int *instance_id, char *name, LDBLE *value, int *status, int c)
{
	*status = 0;

	try
	{
		cd.Reset();

		cd.element.SetString(name);
		cd.value = *value;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_use(int *instance_id, int *use_charge, int *use_phase, int *use_as, int *use_gfw, int *use_unit, int *use_redox, int *status)
{
	*status = 0;

	try
	{
		*status = 1;

		//as
		if (*use_as == 1)
			cd.use[0] = true;
		else 
			cd.use[0] = false;

		//unit
		if (*use_unit == 1)
			cd.use[1] = true;
		else 
			cd.use[1] = false;

		//gfw
		if (*use_gfw == 1)
			cd.use[2] = true;
		else 
			cd.use[2] = false;

		//redox
		if (*use_redox == 1)
			cd.use[3] = true;
		else 
			cd.use[3] = false;

		//charge
		if (*use_charge == 1)
			cd.charge = true;
		else
			cd.charge = false;

		//charge or phase
		if (*use_charge == 1 || *use_phase)
			cd.use[4] = true;
		else 
			cd.use[4] = false;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_as(int *instance_id, char *as, int *status, int c)
{
	*status = 0;

	try
	{
		cd.as.SetString(as);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_gfw(int *instance_id, LDBLE *gfw, int *status)
{
	*status = 0;

	try
	{
		cd.gfw = *gfw;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_unit(int *instance_id, int *unit, int *status)
{
	*status = 0;

	try
	{
		cd.unit = *unit;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_phase(int *instance_id, char *phase_name, LDBLE *si, int *status)
{
	*status = 0;

	try
	{
		cd.phase.SetString(phase_name);
		cd.sat_index = *si;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_redox(int *instance_id, char *element_1, LDBLE *valence_1, char *element_2, LDBLE *valence_2, int *status, int c1, int c2)
{
	*status = 0;

	try
	{
		cd.redox.element_1.SetString(element_1);
		cd.redox.valence_1 = *valence_1;

		cd.redox.element_1.SetString(element_1);
		cd.redox.valence_1 = *valence_1;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_conc_save(int *instance_id, int *conc_id, int *ma_index, LDBLE *gfw, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		int index;
		LDBLE gfw_;
		*conc_id = sol->SetElementConc(cd, index, gfw_);

		if (*conc_id == -1)
		{
			*status = 0;
			return;
		}

		*ma_index = index;
		*gfw = gfw_;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_use_ppa(int *instance_id, int *use, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		model->UsePPA((*use == 1 ? true : false));

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_use_exa(int *instance_id, int *use, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		model->UseExchange((*use == 1 ? true : false));

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_use_sa(int *instance_id, int *use, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		model->UseSurface((*use == 1 ? true : false));

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_use(int *instance_id, int *use_ppa, int *use_gas, int *use_ssa, int *use_sa, int *use_exa, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		model->UsePPA((*use_ppa == 1 ? true : false));
		model->UseGas((*use_gas == 1 ? true : false));
		model->UseSSA((*use_ssa == 1 ? true : false));
		model->UseSurface((*use_sa == 1 ? true : false));
		model->UseExchange((*use_exa == 1 ? true : false));

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}
void pm_ppa_pp(int *instance_id, char *name, char *add_formula, LDBLE *si, LDBLE *moles, int *force_equality, int *dissolve_only, int *pp_id, int *status,	int c1, int c2)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		PPData pd;

		pd.name.SetString(name);
		pd.add_formula.SetString(add_formula);
		pd.si = *si;
		pd.moles = *moles;
		pd.force_equality = (*force_equality == 1 ? true : false);
		pd.dissolve_only = (*dissolve_only == 1 ? true : false);

		int index;
		*pp_id = model->io_ppassemblage->AddPurePhase(pd, index);
//		*pp_index = index;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_exa_exchanger(int *instance_id, char *formula, int *ex_type, char *has, LDBLE *amount, int *index, int *status, int c1, int c2)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;		

		ed.formula.SetString(formula);
		ed.name.SetString(has);
		ed.amount = *amount;
		ed.type = *ex_type;

		*index = model->io_exchange->AddExchanger(&ed);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_sa_options(int *instance_id, int *no_edl, int *diffuse_layer, LDBLE *diffuse_layer_thickness, int *only_counter_ions, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if (*no_edl == 1)
			model->io_surface->SetNoElectrostatic();
		else if (*diffuse_layer == 1)
		{
			model->io_surface->SetDiffuseLayer(*diffuse_layer_thickness);
			model->io_surface->SetOnlyCounterIons((*only_counter_ions == 1 ? true : false));
		}

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_sa_surface(int *instance_id, int *s_type, char *name, char *that_has, LDBLE *specific_area, LDBLE *sites, LDBLE *mass, int *use_area, int *use_mass, int *status, int c1, int c2)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		SurCompData sd;

		if (*s_type == 0)
		{
			sd.formula.SetString(name);
			sd.phase_name.SetString("");
			sd.moles = *sites;
			sd.phase_proportion = 0.0;
		}
		else
		{
			sd.phase_name.SetString(name);
			sd.formula.SetString("");
			sd.phase_proportion = *sites;
			sd.moles = 0.0;
		}
		
		sd.rate_name.SetString(""); //ToDo: Not implemente s_type = 2 yet

		sd.area = *specific_area;
		sd.grams = *mass;

		sd.use_moles = true;
		sd.use_area = (*use_area == 1 ? true : false);
		sd.use_grams = (*use_mass == 1 ? true : false);

		model->io_surface->SetSurComp(sd);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_set_temperature(int *instance_id, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->tc = *new_value;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_set_water_mass(int *instance_id, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->mass_water = *new_value;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_set_ph(int *instance_id, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->ph = *new_value;

		if (model->io_solution->ph_conc_ptr != NULL)
			model->io_solution->ph_conc_ptr->input_conc = *new_value; //ToDo: Not sure if a "ph total" is created if not inputed by user. Check this

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_set_pe(int *instance_id, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->solution_pe = *new_value;

		if (model->io_solution->pe_conc_ptr != NULL)
			model->io_solution->pe_conc_ptr->input_conc = *new_value; //ToDo: Not sure if a "pe total" is created if not inputed by user. Check this

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_solution_set_conc_value(int *instance_id, int *index, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->ChangeElementConc(*index, *new_value);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_exchange_value (int *instance_id, int *index, LDBLE *value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_exchange->comps->Element(*index)->amount = *value;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_phase_value (int *instance_id, int *index, LDBLE *value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_ppassemblage->pure_phases->Element(*index)->moles = *value;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_masterspecies_value(int *instance_id, int *index, LDBLE *value, int *MassOrConcentration, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if (*MassOrConcentration == 0)
			*value = model->GetMasterSpeciesMolality(*index);
		else if (*MassOrConcentration == 1)
			*value = model->GetMasterSpeciesMoles(*index);
		else
			return;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_species_value(int *instance_id, int *index, LDBLE *value, int *MassOrConcentration, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if (*MassOrConcentration == 0)
			*value = model->GetSpeciesMolality(*index);
		else if (*MassOrConcentration == 1)
			*value = model->GetSpeciesActivity(*index);
		else
			return;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_exchange_value(int *instance_id, int *index, LDBLE *value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		*value = model->GetExchangeMoles(*index);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_phase_value(int *instance_id, int *index, LDBLE *value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;
		
		*value = model->GetPhaseMoles(*index);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_species_index(int *instance_id, char *name, int *index, int *status, int c)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		*index = model->GetSpeciesIndex(name);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_get_masterspecies_index(int *instance_id, char *name, int *index, int *status, int c)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		*index = model->GetMasterSpeciesIndex(name);

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_run_model(int *instance_id, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if(!model->RunModel(false, false))
		{
			//model->PrintUserInputToFile();
			return;
		}

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_run_model(int *instance_id, int *print_input, int *print_result, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if(!model->RunModel(IntToBool(*print_input), IntToBool(*print_result)))
		{
			model->PrintUserInputToFile();
			return;
		}

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_run_model(int *instance_id, int *print_input, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if(!model->RunModel(IntToBool(*print_input), false))
		{
			printf("\nPhreeqCLIB - The model has not converged.\n");
			return;
		}

		*status = 1;
	}
	catch(...)
	{
		printf("\nPhreeqCLIB - An undefined error has happened.\n");
		return;
	}
}

void pm_setup_model(int *instance_id, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if (!model->SetupModel())
			return;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}


void pm_read_database(int *instance_id, char *database_path, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		if (!model->LoadDatabase(database_path))
			return;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_kill(int *instance_id, int *status)
{
	*status = 0;

	try
	{
		if (!instances.Delete(*instance_id))
			return;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_input_value(int *instance_id, int *input_id, int *group, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		switch(*group)
		{
		case INT_CONCENTRATION:
			model->io_solution->ChangeElementConc(*input_id, *new_value);
			break;
		case INT_PHASE:
			model->io_ppassemblage->pure_phases->Element(*input_id)->moles = *new_value;
			break;
		case INT_EXCHANGE:
			model->io_exchange->comps->Element(*input_id)->amount = *new_value;
			break;
		default:
			break;
		}
			
		*status = 1;
	}
	catch(...)
	{		
		return;
	}
}

void pm_get_result_value(int *instance_id, int *result_id, int *group, LDBLE *new_value, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		switch(*group)
		{
		case INT_CONCENTRATION:
			*new_value = model->GetMasterSpeciesMoles(*result_id);
			break;
		case INT_PHASE:
			*new_value = model->GetPhaseMoles(*result_id);
			break;
		case INT_EXCHANGE:
			*new_value = model->GetExchangeMoles(*result_id);
			break;
		case INT_SPECIES:
			*new_value = model->GetSpeciesMoles(*result_id);
			break;
		default:
			break;
		}
			
		*status = 1;
	}
	catch(...)
	{		
		return;
	}
}

void pm_get_data(int *instance_id, LDBLE *mass_of_water, LDBLE *ph, LDBLE *pe, int *status)
{
		*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		*mass_of_water = model->md->mass_water_aq_x;
		*ph = -model->gd->s_hplus->la;
		*pe = -model->gd->s_eminus->la;
			
		*status = 1;
	}
	catch(...)
	{		
		return;
	}
}

void pm_set_required_data(int *instance_id, LDBLE *mass_of_water, LDBLE *temp, LDBLE *density, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		model->io_solution->tc = *temp;
		model->io_solution->mass_water = *mass_of_water;

		model->io_solution->density = *density;
		model->io_solution->units = _mol_kgw_;

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_ph(int *instance_id, int *charge, LDBLE *ph, int *status)
{
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		px.value = *ph;

		if (*charge == 1)
			px.charge = true;
		else 
			px.charge = false;

/*
		if (*use_phase == 1)
		{
			//Not working for now...
			//px.phase.SetString(phase_name);
			//px.sat_index = *si;

			px.phase.SetString("");
			px.sat_index = 0.0;
		}
		else
		{*/
			px.phase.SetString("");
			px.sat_index = 0.0;
		//}

		sol->SetPH(px);

		model->io_solution->ph = *ph;
		if (model->io_solution->ph_conc_ptr != NULL)
			model->io_solution->ph_conc_ptr->input_conc = *ph; //ToDo: Not sure if a "ph total" is created if not inputed by user. Check this

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}

void pm_set_pe(int *instance_id, int *charge, LDBLE *pe_value, int *status)
{						 	
	*status = 0;

	try
	{
		model = instances.Element(*instance_id);
		if (model == NULL)
			return;

		sol = model->io_solution;

		px.value = (LDBLE)*pe_value;

		if (*charge == 1)
			px.charge = true;
		else 
			px.charge = false;

/*		if (*use_phase == 1)
		{
			//Not working for now...
			//px.phase.SetString(phase_name);
			//px.sat_index = *si;

			px.phase.SetString("");
			px.sat_index = 0.0;
		}
		else
		{*/
			px.phase.SetString("");
			px.sat_index = 0.0;
		//}

		sol->SetPE(px);

		model->io_solution->solution_pe = *pe_value;
		if (model->io_solution->pe_conc_ptr != NULL)
			model->io_solution->pe_conc_ptr->input_conc = *pe_value; //ToDo: Not sure if a "pe total" is created if not inputed by user. Check this

		*status = 1;
	}
	catch(...)
	{
		return;
	}
}