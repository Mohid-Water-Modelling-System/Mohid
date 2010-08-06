#ifndef tidyH
#define tidyH

#include "modeldata.h"
#include "check.h"
#include "datastructures.h"
#include "parseengine.h"
#include "debug.h"
#include "units.h"
#include "ssassemblagemodel.h"

class Tidy:
	public Check,
	public SSAssemblageModel,
	public Units
{
public:
	Tidy(GlobalData *gd, ModelData *md);
	~Tidy();

public:
	bool TidySpecies();
	bool TidyPhases();

	bool TidySolution(Solution *sol);
	bool TidyPurePhase(PPAssemblage *ppa);
	bool TidyGasPhase(GasPhase *gas);
	bool TidySSAssemblage(SSAssemblage *ssa);
	bool TidySurface(Surface *sur);
	bool TidyMinExchange(Exchange *exc);
	bool TidyMinSurface(Surface *sur);

	bool CheckObrigatorySpecies();

protected:
	bool CheckEq(bool association);
	bool RewriteEqToPrimary();
	bool CoefInMaster(Master *m, LDBLE &coef);
	bool RewriteEqToSecondary();
	bool CalcAlk(Reaction *reaction, LDBLE &alk);
	Master * SearchForMasterPrimary(String &name);
	bool ReplaceSolidsGases(bool &replace);

protected:
	Reaction *rt_list;
	Debug d;

private:
	GlobalData *gd;
	ModelData *md;

	Surface *sur_p;
	Exchange *exc_p;
	PPAssemblage *ppa_p;
};

#endif