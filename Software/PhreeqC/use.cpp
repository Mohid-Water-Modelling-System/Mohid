#include "modeldata.h"

//===================================================================================
// Contructs
//===================================================================================
Use::Use()
{
	sol_list = NULL;
	ppa_list = NULL;
	exc_list = NULL;
	sur_list = NULL;
	gas_list = NULL;
	ssa_list = NULL;
	//irr_list = NULL;

	sol_list = new List<Solution>;
	ppa_list = new List<PPAssemblage>;
	exc_list = new List<Exchange>;
	sur_list = new List<Surface>;
	gas_list = new List<GasPhase>;
	ssa_list = new List<SSAssemblage>;
	//irr_list = new List<Irrev>;

	sol_number = 0;
	ppa_number = 0;
	exc_number = 0;
	sur_number = 0;
	gas_number = 0;
	ssa_number = 0;
	irr_number = 0;

	sol_in = true;
	exc_in = false;
	ppa_in = false;
	gas_in = false;
	sur_in = false;
	ssa_in = false;
	irr_in = false;

	Reset(); 
}
//-----------------------------------------------------------------------------------	
Use::~Use()
{
	// Will free memory associated with the lists AND free memory of every item on the lists
	delete sol_list;
	delete ppa_list;
	delete exc_list;
	delete sur_list;
	delete gas_list;
	delete ssa_list;
	//delete irr_list;
}

//===================================================================================
// Interface
//===================================================================================
void Use::Reset()
{
	sol_p = NULL;

	ResetPointers();
	ClearLists();
}
//-----------------------------------------------------------------------------------	
void Use::ResetPointers()
{
	ppa_p = NULL;
	exc_p = NULL;
	sur_p = NULL;
	gas_p = NULL;
	ssa_p = NULL;
	//irr_p = NULL;
}
//-----------------------------------------------------------------------------------	
bool Use::CheckIn()
{
	if (ppa_in || exc_in || sur_in || gas_in || ssa_in || irr_in)
		return true;

	return false;
}
//-----------------------------------------------------------------------------------	
void Use::ClearLists()
{
	sol_list->Clear();
	ppa_list->Clear();
	exc_list->Clear();
	sur_list->Clear();
	gas_list->Clear();
	ssa_list->Clear();
	//irr_list->Clear();
}