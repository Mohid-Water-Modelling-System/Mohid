#ifndef unitsH
#define unitsH

#include "constants.h"

class Units
{
public:
	bool Check_l_(UNITS units);
	bool Check_g_l_(UNITS units);
	bool Check_g_kgs_(UNITS units);
	bool Check_g_(UNITS units);
	bool Check_kgs_(UNITS units);
	bool Check_mol_l_(UNITS units);
	bool Check_mol_kgs_(UNITS units);
	bool Check_eq_l_(UNITS units);
	LDBLE ConversionFactorForMoles(UNITS units);
	bool CheckUnit(UNITS unit, UNITS solution_units);
	bool CheckUnitForAlk(UNITS unit, UNITS solution_units);
	bool CheckUnitCompatibility(UNITS units, UNITS solution_units);

};

#endif