#include "model.h"
#include <math.h>
//#include "checkdatabase.h"

//===================================================================================
// Constructor and Destructor
//===================================================================================
Model::Model():
ParseEngine(NULL), IOEngine(), Messages("phreeqc_input.err")
{
	t = NULL;

	gd = NULL;
	db = NULL;
	md = NULL;

	me = NULL;

	error = NULL;

	gd = new GlobalData;
	db = new Database(gd);
	md = new ModelData(gd);

	t = new Tidy(gd, md);

	me = new ModelEngine(gd, md);
	
	ParseEngine::SetGlobalData(gd);

	IOEngine::SetGlobalData(gd);
	IOEngine::SetModelData(md);
	IOEngine::SetModelEngine(me);

	solution_r = NULL;

	output_results = NULL;
	options = new ModelOptions;

#ifdef DEBUG_MOHID
	String filename("speciesinfo.txt");
	me->d.OpenSpeciesInfoFile(filename);
	me->d.OpenGammasFile();
	me->d.OpenDeltaFile();
	me->d.OpenArrayFile("arr.txt");
	me->d.OpenLAFile();
	me->d.OpenDebugFile();
	me->d.OpenSetFile();
	me->d.OpenReviseGuessesFile();
	me->d.OpenBuildModelFile();
	me->d.OpenMolalitiesFile();
	me->d.OpenMBSumsFile();
	me->d.OpenBuildSSAssemblageFile();
	me->d.OpenBuildJacobianSumsFile();
	me->d.OpenIneqFile();
	me->d.OpenResetFile();
	me->d.OpenSetupSolutionFile();
#endif

#ifdef DEBUG_AddSolution
	me->d.OpenAddSolutionFile();
#endif
#ifdef DEBUG_AddExchange
	me->d.OpenAddExchangeFile();
#endif
#ifdef DEBUG_AddPPAssemblage
	me->d.OpenAddPPAssemblageFile();
#endif
}


Model::Model(Model *copy):
ParseEngine(NULL), IOEngine(), Messages("phreeqc_input.err")
{
	t = NULL;

	gd = NULL;
	db = NULL;
	md = NULL;

	me = NULL;

	error = NULL;

	gd = new GlobalData;
	db = new Database(gd);
	md = new ModelData(gd);

	t = new Tidy(gd, md);

	me = new ModelEngine(gd, md);

	ParseEngine::SetGlobalData(gd);
	IOEngine::SetGlobalData(gd);
	IOEngine::SetModelData(md);
	IOEngine::SetModelEngine(me);

	solution_r = NULL;

	output_results = NULL;
	options = new ModelOptions;

#ifdef DEBUG_MOHID
	String filename("speciesinfo.txt");
	me->d.OpenSpeciesInfoFile(filename);
	me->d.OpenGammasFile();
	me->d.OpenDeltaFile();
	me->d.OpenArrayFile("arr.txt");
	me->d.OpenLAFile();
	me->d.OpenDebugFile();
	me->d.OpenSetFile();
	me->d.OpenReviseGuessesFile();
	me->d.OpenBuildModelFile();
	me->d.OpenMolalitiesFile();
	me->d.OpenMBSumsFile();
	me->d.OpenBuildSSAssemblageFile();
	me->d.OpenBuildJacobianSumsFile();
	me->d.OpenIneqFile();
	me->d.OpenResetFile();
	me->d.OpenSetupSolutionFile();
#endif

#ifdef DEBUG_AddSolution
	me->d.OpenAddSolutionFile();
#endif
#ifdef DEBUG_AddExchange
	me->d.OpenAddExchangeFile();
#endif
#ifdef DEBUG_AddPPAssemblage
	me->d.OpenAddPPAssemblageFile();
#endif
}

Model::~Model()
{
#ifdef DEBUG_AddSolution
	me->d.CloseAddSolutionFile();
#endif
#ifdef DEBUG_AddExchange
	me->d.CloseAddExchangeFile();
#endif
#ifdef DEBUG_AddPPAssemblage
	me->d.CloseAddPPAssemblageFile();
#endif

	delete gd;
	delete t;

#ifdef DEBUG_MOHID
	me->d.CloseSpeciesInfoFile();
	me->d.CloseGammasFile();
	me->d.CloseDeltaFile();
	me->d.CloseArrayFile();
	me->d.CloseLAFile();
	me->d.CloseDebugFile();
	me->d.CloseSetFile();
	me->d.CloseReviseGuessesFile();
	me->d.CloseBuildModelFile();
	me->d.CloseMolalitiesFile();
	me->d.CloseMBSumsFile();
	me->d.CloseBuildSSAssemblageFile();
	me->d.CloseBuildJacobianSumsFile();
	me->d.CloseIneqFile();
	me->d.CloseResetFile();
	me->d.CloseSetupSolutionFile();
#endif

	delete me;

	delete md;
	delete db;
	
	if (output_results != NULL)
		fclose(output_results);

	delete options;
	delete error;
}

//===================================================================================
// Public functions
//===================================================================================
bool Model::LoadDatabase(String file_name)
{
	if (!db->Load(file_name))
		return false;

	return true;
}

bool Model::SetupModel()
{
	if (!db->TidyDatabase())
		return false;

	IOEngine::InitSolutionPE();

	if (options->print_always == 1)
	{
		//output_results = fopen("phreeqc_output.txt", "wt");
	}

	return true;
}

bool Model::RunModel(bool print_input, bool print_result)
{
	try
	{
		solution = io_solution;

		if (md->use.ppa_in)
		{
			ppassemblage = io_ppassemblage;
			//ppassemblage->pure_phases->Sort();
		}

		if (md->use.gas_in)
			gasphase = io_gasphase;

		if (md->use.sur_in)
		{
			surface = io_surface;
			SetupSurface();
		}

		if (md->use.ssa_in)
			ssassemblage = io_ssassemblage;

		if (md->use.exc_in)
		{
			exchange = io_exchange;
			SetupExchange();
		}

		if (print_input)
		{
			output_results = fopen("phreeqc_output.txt", "wt");
			PrintUserInputToFile(output_results);
			fclose(output_results);
			output_results = NULL;
		}

		if (!ResetModel())
			return false;

		
		// Calculate distribution of species for initial solutions

		if (!RunInitialSolutionCalculations())
			return false;

		// Calculate distribution for exchangers
		if (exc_exists && !RunInitialExchangersCalculations())
				return false;

		// Calculate distribution for surfaces
		if (sur_exists && !RunInitialSurfacesCalculations())
				return false;

		// Calculate initial gas composition
		if (md->use.gas_in && !RunInitialGasPhaseCalculations())
				return false;

		// Calculate reactions
	#ifdef DEBUG_MOHID
		me->d.gammas_count = 0;
	#endif


	#ifdef DEBUG_MOHID
			me->d.debug_status = true;
	#endif
		if (!RunReactions())
			return false;
	#ifdef DEBUG_MOHID
			me->d.debug_status = false;
	#endif

		//me->SaveSolutionResults();	

		//if (print_result)
		//	PrintResultToFile(output_results);

		me->SavePPResults(io_ppassemblage); //Save the PP results (moles and deltas)
	}
	catch(...)
	{
		return false;
	}

	return true;
}

//===================================================================================
// Functions called directly by RunModel
//===================================================================================
bool Model::RunInitialSolutionCalculations()
{
	//
	// make initial solution calculations
	//

	CONVERGE_RESULT converge; 
	Use *use = &md->use;

  md->state = INITIAL_SOLUTION;

	if (!SetUse())
		return false;

	md->dl_type_x = NO_DL;
	use->sol_p = (*use->sol_list)[0];

  if (!me->Prepare())
		return false;
 
	if (!me->KTemp(md->use.sol_p->tc))
		return false;

	if (!me->Set(true))
		return false;

	converge = me->ExecuteModel();

	if (converge == ERROR_CR) 
		return false;
		
	if (converge != OK_CR && !md->diagonal_scale)
  {
		md->diagonal_scale = true;
		me->Set(true);
		converge = me->ExecuteModel();
		md->diagonal_scale = false;
  }
     
	if (converge != OK_CR || !me->CheckResiduals())
		return false;

  if (!me->SumSpecies ()) 
		return false;    

	me->SaveSolution(use->sol_p); //ToDo: Check if this is ok.
	solution_r = use->sol_p;

  return true;
}
//-----------------------------------------------------------------------------------
bool Model::RunInitialExchangersCalculations()
{
	//
	// Number of solution to use must be passed in "solution number"
	//

  //int i;
	CONVERGE_RESULT converge; //, converge_residuals;
	Use *use;
  //int n;
  //char token[2 * MAX_LENGTH];
	Exchange *exc;

	use = &md->use;
  md->state = INITIAL_EXCHANGE;
	
  SetUse();

  md->dl_type_x = NO_DL;

	if (!md->use.exc_in)
		return true;

	exc = (*use->exc_list)[0]; 

	if (exc->solution_equilibria)
  {
		use->exc_p = exc; 
		use->sol_p = (*use->sol_list)[0];		

		if (!me->Prepare())
			return false;

		if (!me->KTemp(use->sol_p->tc))
			return false;

		if (!me->Set(true))
			return false;
	
		converge = me->ExecuteModel();

		if (converge == ERROR_CR)
			return false;

		if (!me->CheckResiduals())
			return false;
	
		if (!me->SumSpecies())
			return false;

		md->species_info_list->Sort(gd->s_hplus);
    
		me->SaveExchange(exc); //xexchange_save (n_user); //ToDo: Verificar se é necessária
  }
	
  return true;
}
//-----------------------------------------------------------------------------------
bool Model::RunInitialGasPhaseCalculations()
{
	//
	// Go through list of gas_phases, make initial calculations
	// for any marked "new" that are defined to be in equilibrium with a 
	// solution.
	//

  int i, j;
	CONVERGE_RESULT converge; //, converge1;
  //int n, last, n_user, print1;
  //char token[2 * MAX_LENGTH];
	GasPhase *gas_p;
	GasComp *gas_comp_ptr;
	Phase *phase_ptr;
	ReactionToken *rxn_ptr;
  LDBLE lp;
	Use *use;

	use = &md->use;
  md->state = INITIAL_GAS_PHASE;
  
	SetUse();

  md->dl_type_x = NO_DL;

	if (!md->use.gas_in)
		return true;

	gas_p = (*use->gas_list)[0];

	if (gas_p->solution_equilibria)
  {
		use->sol_p = (*use->sol_list)[0];

		if (!me->Prepare())
			return false;

		if (!me->KTemp(md->use.sol_p->tc))
			return false;

		if (!me->Set(true))
			return false;
	
		converge = me->ExecuteModel();

		if (converge == ERROR_CR)
			return false;

		if (!me->CheckResiduals())
			return false;
	
		use->gas_p = gas_p;

		gas_p->total_p = 0;
		gas_p->total_moles = 0;

		for (i = 0; i < gas_p->comps->Count(); i++)
    {
			gas_comp_ptr = (*gas_p->comps)[i];
			phase_ptr = gas_comp_ptr->phase;
	
			if (phase_ptr->in)
			{
				lp = -phase_ptr->lk;
	  
				for (j = 1; j < phase_ptr->rxn_x->token_list->Count(); j++)
				{
					rxn_ptr = phase_ptr->rxn_x->token_list->Element(j);
					lp += rxn_ptr->s->la * rxn_ptr->coef;
				}

				gas_comp_ptr->phase->p_soln_x = exp (lp * md->LOG_10);
	  		gas_p->total_p += gas_comp_ptr->phase->p_soln_x;
				gas_comp_ptr->phase->moles_x = gas_comp_ptr->phase->p_soln_x * gas_p->volume / ((LDBLE)R_LITER_ATM * md->tk_x);
				gas_comp_ptr->moles = gas_comp_ptr->phase->moles_x;
				gas_p->total_moles += gas_comp_ptr->moles;
			}
			else
			{
				gas_comp_ptr->phase->moles_x = 0;
			}
		}
  
		me->SaveGasPhase(gas_p); //xgas_save (n_user); //ToDo: Verificar se é necessária		
  }

	return true;
}
//-----------------------------------------------------------------------------------
bool Model::RunReactions()
{
	//
	// Make all reaction calculation which could include:
  //    equilibrium with a pure-phase assemblage,
  //    equilibrium with an exchanger,
  //    equilibrium with an surface,
  //    equilibrium with a gas phase,
  //    equilibrium with a solid solution assemblage,
  //    or irreversible reaction.
	//
	int count_steps;

  md->state = REACTION;

	if (!SetUse())
		return true;

	// Find maximum number of steps
	count_steps = 1;

	me->CopyUse(2);

	/*
	if (md->use.irr_in && md->use.irr_p != NULL)
		if (abs(md->use.irr_p->count_steps) > count_steps)
			count_steps = abs(md->use.irr_p->count_steps);
	*/

	md->count_total_steps = count_steps;

	// save data for saving solutions
  //memcpy (&save_data, &save, sizeof (struct save));

	for(md->reaction_step = 1; md->reaction_step <= count_steps; md->reaction_step++)
	{
    if (md->reaction_step > 1)
    {
      me->CopyUse(2); //CopyUse must be done
    }

		me->SetInitialMoles(2);

		//Run reaction step

		me->RunReactions(2, 1.0);

		// REACTION is a keyword that will be "deactivated" on this first version
    //if (md->reaction_step < count_steps)
    //{
      //me->Saver();
    //}

	}
	
	//Saver(); //Check if it's necessary, or if do some calculation that is necessary

	solution_r = md->use.sol_list->Element(2);
	return true;
}
//-----------------------------------------------------------------------------------
bool Model::PrintResults(String filename)
{
	me->PrintResults(filename);
	return true;
}
//===================================================================================
// Private functions
//===================================================================================
bool Model::SetUse()
{
	Use *use = &md->use;

	use->ResetPointers();

  if (md->state < REACTION)
    return true;
	
	if (!use->CheckIn())
		return false;

	use->sol_p = (*use->sol_list)[use->sol_number];

	if (use->ppa_in)
		use->ppa_p = (*use->ppa_list)[use->ppa_number];

	if (use->exc_in)
		use->exc_p = (*use->exc_list)[use->exc_number];

	md->dl_type_x = NO_DL;

	if (use->sur_in)
		use->sur_p = (*use->sur_list)[use->sur_number];

	if (use->gas_in)
		use->gas_p = (*use->gas_list)[use->gas_number];

	if (use->ssa_in)
		use->ssa_p = (*use->ssa_list)[use->ssa_number];

	return true;
}

//-----------------------------------------------------------------------------------

bool Model::RunInitialSurfacesCalculations()
{
	//
	// Go through list of surfaces, make initial calculations
	// for any marked "new" that are defined to be in equilibrium with a 
	// solution.
	//

  //int i;
  //int n, last, n_user, print1;
  //char token[2 * MAX_LENGTH];
	Use *use;
	Surface *sur;

	use = &md->use;

  md->state = INITIAL_SURFACE;

  SetUse ();

	if (!md->use.sur_in)
		return true;

	sur = (*use->sur_list)[0];

  if (sur->solution_equilibria)
  {
    use->sur_p = sur;
    md->dl_type_x = md->use.sur_p->dl_type;
		use->sol_p = (*use->sol_list)[0];

    me->SetAndRunWrapper(1, 0.0, 1); //(-1, FALSE, FALSE, -1, 0.0); //ToDo: The 1 argument enter in the place of -1.
		md->species_info_list->Sort(gd->s_hplus);
    
		me->SaveSurface(use->sur_p); //xsurface_save (n_user); //ToDo: Check this
  }

  return true;
}

bool Model::ResetModel()
{
	//
	// This functions "reset" all the data to start the model, and tidy solution and others
	//

	// reset global data 
	solution_r = NULL;
	gd->ResetToExecuteModel();

	// reset model data and use data
	md->Reset();

	// reset model engine data
	//me->ResetData();

	// copy original data input to number 0 on the use lists
	SaveDataStructures();

	if (md->use.sur_in && !t->TidySurface((*md->use.sur_list)[0]))
		return false;

	if (md->use.gas_in && !t->TidyGasPhase((*md->use.gas_list)[0]))
		return false;

	if (md->use.ppa_in && !t->TidyPurePhase((*md->use.ppa_list)[0]))
		return false;

	if (md->use.ssa_in && !t->TidySSAssemblage((*md->use.ssa_list)[0]))
		return false;

	if (md->use.exc_in && !t->TidyMinExchange((*md->use.exc_list)[0]))
		return false;

	if (md->use.sur_in && !t->TidyMinSurface((*md->use.sur_list)[0]))
		return false;

	if (!t->TidySolution((*md->use.sol_list)[0]))
		return false;

	if (!t->CheckObrigatorySpecies())
		return false;

	return true;
}

void Model::SaveDataStructures()
{
	Use *use = &md->use;

	use->ClearLists();

	Solution *new_solution = use->sol_list->AddNew();
	solution->CopyTo(new_solution);
	new_solution->number = 0;

	if (use->exc_in)
	{
		Exchange *new_exchange = use->exc_list->AddNew();
		exchange->CopyTo(new_exchange);
		new_exchange->number = 0;
	}

	if (use->sur_in)
	{
		Surface *new_surface = use->sur_list->AddNew();
		surface->CopyTo(new_surface);
		new_surface->number = 0;
	}

	if (use->gas_in)
	{
		GasPhase *new_gas = use->gas_list->AddNew();
		gasphase->CopyTo(new_gas);
		new_gas->number = 0;
	}

	if (use->ssa_in)
	{
		SSAssemblage *new_ssassemblage = use->ssa_list->AddNew();
		ssassemblage->CopyTo(new_ssassemblage);
		new_ssassemblage->number = 0;
	}

	if (use->ppa_in)
	{
		PPAssemblage *new_ppassemblage = use->ppa_list->AddNew();
		ppassemblage->CopyTo(new_ppassemblage);
		new_ppassemblage->number = 0;
	}
}

void Model::SetupSurface()
{
	LDBLE conc;
	char token[DEFAULT_STR_LENGTH];
	char *ptr;
	String token1, name;
	int count_charge = 0, j;

	SurfaceCharge *sur_charge;
	SurfaceComp *sur_comp;
	for (int i = 0; i < surface->comps->Count(); i++)
	{
		sur_comp = (*surface->comps)[i];

		if (!sur_comp->use_moles)
			conc = 1.0;
		else
			conc = sur_comp->moles;

		eos_list->Clear();
		parent_count = 0;

		// Accumulate elements in elt_list
		sur_comp->formula.Copy(token);
		ptr = token;
		GetElementsInSpecies(&ptr, conc);

		// save formula for adjusting number of exchange sites
		ptr = token;
		GetSpeciesNameAndZ(&ptr, token1, sur_comp->formula_z);

		// Search for charge structure
		ptr = token;
		GetElement(&ptr, name);

		SaveEOSList(sur_comp->formula_totals);
		SaveEOSList(sur_comp->totals);

		name = name(0, "_");

		for (j = 0; j < count_charge; j++)
		{
			sur_charge = (*surface->charge)[j];
			if (name == sur_charge->name)
				break;
		}

		if (j >= count_charge)
		{
			sur_charge = surface->charge->AddNew();

			sur_charge->name = name;

			if (sur_comp->phase_name.IsEmpty() && sur_comp->rate_name.IsEmpty())
			{
				sur_charge->specific_area = 600.0;
				sur_charge->grams = 0.0;
			}
			else
			{
				sur_charge->specific_area = 0.0;
				sur_charge->grams = 1.0;
			}

			sur_comp->charge = j;

			if (sur_comp->use_area)
				sur_charge->specific_area = sur_comp->area;

			if (sur_comp->use_grams)
				sur_charge->grams = sur_comp->grams;
		}
	}	
}

void Model::SetupExchange()
{
	LDBLE conc;
	char token[DEFAULT_STR_LENGTH];
	char *ptr;
	String token1, name;
	ExchComp *exc;

	for (int i = 0; i < exchange->comps->Count(); i++)
	{
		exc = (*exchange->comps)[i];

		if (exc->type > 0)
			conc = 1.0;
		else
			conc = exc->amount;

		eos_list->Clear();
		parent_count = 0;

		// Accumulate elements in elt_list
		exc->formula.Copy(token);
		ptr = token;
		GetElementsInSpecies(&ptr, conc);

		// save formula for adjusting number of exchange sites
		ptr = token;
		GetSpeciesNameAndZ(&ptr, token1, exc->formula_z);

		SaveEOSList(exc->formula_totals);
		SaveEOSList(exc->totals);

		exc->moles = conc;
		exc->charge_balance = 0.0;
	}
}

void Model::PrintUserInputToFile()
{
	Open(); //Messages::Open -> open file to write user input

	solution->Print(message_file);
	Write("\n");

	if (md->use.exc_in)
		exchange->Print(message_file);

	if (md->use.ppa_in)
		ppassemblage->Print(message_file);

	Close(); //Messages::Close
}

void Model::SetOptions(int option, bool state)
{
	switch(option)
	{
	case 0: 
		options->print_always = state;
		break;
	}
}

void Model::PrintUserInputToFile(FILE *file)
{
	solution->Print(file);
	Write("\n");

	if (md->use.exc_in)
		exchange->Print(file);

	if (md->use.ppa_in)
		ppassemblage->Print(file);

	Write("\n------------------------------------------------");
}