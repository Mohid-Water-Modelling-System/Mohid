#include "IOEngine.h"

IOEngine::IOEngine()
{
	md = NULL;
	me = NULL;
	gd = NULL;

	io_solution     = NULL;
	io_ppassemblage = NULL;
	io_gasphase     = NULL;
	io_surface      = NULL;
	io_ssassemblage = NULL;
	io_exchange     = NULL;

	io_solution     = new Solution(gd);
	io_ppassemblage = new PPAssemblage(gd);
	io_gasphase     = new GasPhase;
	io_surface      = new Surface;
	io_ssassemblage = new SSAssemblage;
	io_exchange     = new Exchange;
}

IOEngine::~IOEngine()
{
	delete io_solution;
	delete io_ppassemblage;
	delete io_gasphase;
	delete io_surface;
	delete io_ssassemblage;
	delete io_exchange;
}

void IOEngine::InitSolutionPE() 
{ 
	io_solution->gd = gd;
	io_solution->InitPE();
}

LDBLE IOEngine::GetMasterSpeciesMoles(int index)
{
	if (index < 0 || index > gd->master_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetMasterSpeciesMoles): Master Species not found.");

	return gd->master_list.Element(index)->total;
}

LDBLE IOEngine::GetMasterSpeciesMilligrams(int index)
{
	if (index < 0 || index > gd->master_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetMasterSpeciesMoles): Master Species not found.");

	return gd->master_list.Element(index)->total * gd->master_list.Element(index)->gfw;
}

LDBLE IOEngine::GetMasterSpeciesMolality(int index)
{
	if (index < 0 || index > gd->master_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetMasterSpeciesMolality): Master Species not found.");

	return gd->master_list.Element(index)->total / md->mass_water_aq_x;
}

LDBLE IOEngine::GetSpeciesMolality(int index)
{
	if (index < 0 || index > gd->species_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetSpeciesMolality): Species not found.");

	return gd->species_list.Element(index)->moles / md->mass_water_aq_x;
}

LDBLE IOEngine::GetSpeciesMoles(int index)
{
	if (index < 0 || index > gd->species_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetSpeciesMolality): Species not found.");

	return gd->species_list.Element(index)->moles;
}

LDBLE IOEngine::GetSpeciesActivity(int index)
{
	if (index < 0 || index > gd->species_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetSpeciesActivity): Species not found.");

	Species *s = gd->species_list.Element(index);
	return me->Under(s->lm + s->lg);
}

LDBLE IOEngine::GetExchangeMoles(int index)
{
	if (index < 0 || index > gd->species_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetExchangeMoles): Exchange not found.");	

	return gd->species_list.Element(index)->moles;
}

LDBLE IOEngine::GetExchangeMilligrams(int index)
{
	if (index < 0 || index > gd->species_list.Count())
		throw ExceptionHandler("PhreeqC Object (GetExchangeMoles): Exchange not found.");	

	//ToDo: Check if the gfw value is ok.
	return gd->species_list.Element(index)->moles * gd->species_list.Element(index)->gfw;
}

LDBLE IOEngine::GetPhaseMoles(int index)
{
	if (index < 0 || index > io_ppassemblage->pure_phases->Count())
		throw ExceptionHandler("PhreeqC Object (GetPhaseMoles): Phase not found.");	

	return io_ppassemblage->pure_phases->Element(index)->moles;
}

int IOEngine::GetSpeciesIndex(char *name)
{
	int index;
	gd->species_list.Search(name, index, true);

	if (index < 0)
		throw ExceptionHandler("PhreeqC Object (GetSpeciesIndex): Species not found.");

	return index;
}

int IOEngine::GetMasterSpeciesIndex(char *name)
{
	int index;
	gd->master_list.Search(name, index, true);

	if (index < 0)
		throw ExceptionHandler("PhreeqC Object (GetMasterIndex): Master Species not found.");

	return index;
}