#include "units.h"

bool Units::Check_kgs_(UNITS units)
{
	switch (units)
	{
	case _mol_kgs_:
	case _mmol_kgs_:
	case _umol_kgs_:
	case _g_kgs_:
	case _mg_kgs_:
	case _ug_kgs_:
	case _eq_kgs_:
	case _meq_kgs_:
	case _ueq_kgs_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_l_(UNITS units)
{
	switch (units)
	{
	case _mol_l_:
	case _mmol_l_:
	case _umol_l_:
	case _g_l_:
	case _mg_l_:
	case _ug_l_:
	case _eq_l_:
	case _meq_l_:
	case _ueq_l_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_g_l_(UNITS units)
{
	switch (units)
	{
	case _g_l_:
	case _mg_l_:
	case _ug_l_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_g_kgs_(UNITS units)
{
	switch (units)
	{
	case _g_kgs_:
	case _mg_kgs_:
	case _ug_kgs_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_g_(UNITS units)
{
	switch (units)
	{
	case _g_kgs_:
	case _mg_kgs_:
	case _ug_kgs_:
	case _g_kgw_:
	case _mg_kgw_:
	case _ug_kgw_:
	case _g_l_:
	case _mg_l_:
	case _ug_l_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_mol_l_(UNITS units)
{
	switch (units)
	{
	case _mol_l_:
	case _mmol_l_:
	case _umol_l_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_mol_kgs_(UNITS units)
{
	switch (units)
	{
	case _mol_kgs_:
	case _mmol_kgs_:
	case _umol_kgs_:
		return true;
	default:
		return false;
	}
}

bool Units::Check_eq_l_(UNITS units)
{
	switch (units)
	{
	case _eq_l_:
	case _meq_l_:
	case _ueq_l_:
		return true;
	default:
		return false;
	}
}

LDBLE Units::ConversionFactorForMoles(UNITS units)
{
	switch (units)
	{
	case _mmol_l_:
	case _mg_l_:
	case _meq_l_:
	case _mmol_kgs_:
	case _mg_kgs_:
	case _meq_kgs_:
	case _mmol_kgw_:
	case _mg_kgw_:
	case _meq_kgw_:
		return (LDBLE)1e-3;
	case _umol_l_:
	case _ug_l_:
	case _ueq_l_:
	case _umol_kgs_:
	case _ug_kgs_:
	case _ueq_kgs_:
	case _umol_kgw_:
	case _ug_kgw_:
	case _ueq_kgw_:
		return (LDBLE)1e-6;
	default:
		return 1.0;
	}
}

bool Units::CheckUnit(UNITS unit, UNITS solution_units)
{
	switch(unit)
	{
	case _eq_l_:
	case _meq_l_:
	case _ueq_l_:
	case _eq_kgw_:
	case _meq_kgw_:
	case _ueq_kgw_:
	case _eq_kgs_:
	case _meq_kgs_:
	case _ueq_kgs_:
		return false;
	default:
		return CheckUnitCompatibility(unit, solution_units);
	}
}

bool Units::CheckUnitForAlk(UNITS unit, UNITS solution_units)
{
	switch(unit)
	{
	case _mol_l_:
	case _mmol_l_:
	case _umol_l_:
	case _mol_kgw_:
	case _mmol_kgw_:
	case _umol_kgw_:
	case _mol_kgs_:
	case _mmol_kgs_:
	case _umol_kgs_:
		return false;
	default:
		return CheckUnitCompatibility(unit, solution_units);
	}
}

bool Units::CheckUnitCompatibility(UNITS units, UNITS solution_units)
{
	int units_range, default_units_range;

	switch (units)
	{
	case _mol_l_:
	case _mmol_l_:
	case _umol_l_:
	case _g_l_:
	case _mg_l_:
	case _ug_l_:
	case _eq_l_:
	case _meq_l_:
	case _ueq_l_:
		units_range = 0;
		break;
	case _mol_kgs_:
	case _mmol_kgs_:
	case _umol_kgs_: 
	case _g_kgs_:
	case _mg_kgs_:
	case _ug_kgs_: 
	case _eq_kgs_:
	case _meq_kgs_:
	case _ueq_kgs_:
		units_range = 1;
		break;
	default:
		units_range = 2;
		break;
	}

	switch (solution_units)
	{
	case _mol_l_:
	case _mmol_l_:
	case _umol_l_:
	case _g_l_:
	case _mg_l_:
	case _ug_l_:
	case _eq_l_:
	case _meq_l_:
	case _ueq_l_:
		default_units_range = 0;
		break;
	case _mol_kgs_:
	case _mmol_kgs_:
	case _umol_kgs_: 
	case _g_kgs_:
	case _mg_kgs_:
	case _ug_kgs_: 
	case _eq_kgs_:
	case _meq_kgs_:
	case _ueq_kgs_:
		default_units_range = 1;
		break;
	default:
		default_units_range = 2;
		break;
	}

	if (default_units_range != units_range)
		return false;

	return true;
}
