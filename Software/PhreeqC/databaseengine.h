#ifndef databaseengineH
#define databaseengineH

#include "definitions.h"
#include "datastructures.h"
#include "line.h"
#include "tidy.h"

class DBEngine:
	public Tidy
{
public:
	DBEngine(GlobalData *gd);
	~DBEngine();

public:
	bool GetValueFromNextToken(LDBLE &value);
	bool ReadDeltaH(LDBLE &delta_h, DELTA_H_UNIT &units);
	bool ReadAnalyticalExpression(LDBLE *logk);
	bool ReadAddLogK(List<NameCoef> *list);
	bool ReadAddConstant(List<NameCoef> *list);
	bool ReadMassBalance(Species *s);

protected:
	bool LoadSecondarySpecies(String str, Species *s);
	bool AddPSIMasterSpecies(String &str);

protected:
	Line line;	
	GlobalData *gd;
	LDBLE gfw_water;
};

#endif