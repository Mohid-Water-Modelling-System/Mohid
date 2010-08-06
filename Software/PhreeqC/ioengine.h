#ifndef ioengineH
#define ioengineH

#include "modeldata.h"
#include "modelengine.h"

class IOEngine
{
public:
	IOEngine();
	~IOEngine();

public:
	Solution     *io_solution;
	PPAssemblage *io_ppassemblage;
	GasPhase     *io_gasphase;
	Surface      *io_surface;
	SSAssemblage *io_ssassemblage;
	Exchange     *io_exchange;

public:
	void UsePPA(bool use = true) { md->use.ppa_in = use; }
	void UseGas(bool use = true) { md->use.gas_in = use; }
	void UseSurface(bool use = true) { md->use.sur_in = use; }
	void UseSSA(bool use = true) { md->use.ssa_in = use; }
	void UseExchange(bool use = true) { md->use.exc_in = use; }

	bool ThereIsUse() { return (md->use.ppa_in || md->use.gas_in || md->use.sur_in || md->use.ssa_in || md->use.exc_in); }

public:
	LDBLE GetMasterSpeciesMoles(int index);
	LDBLE GetMasterSpeciesMilligrams(int index);
	LDBLE GetMasterSpeciesMolality(int index);
	LDBLE GetSpeciesMoles(int index);
	LDBLE GetSpeciesMolality(int index);
	LDBLE GetSpeciesActivity(int index);
	LDBLE GetExchangeMoles(int index);
	LDBLE GetExchangeMilligrams(int index);
	LDBLE GetPhaseMoles(int index);
	int GetSpeciesIndex(char *name);
	int GetMasterSpeciesIndex(char *name);

protected:
	void InitSolutionPE();
	void SetGlobalData(GlobalData *gd) { this->gd = gd; }
	void SetModelData(ModelData *md) { this->md = md; }
	void SetModelEngine(ModelEngine *me) { this->me = me; }

private:
	ModelData *md;
	ModelEngine *me;
	GlobalData *gd;

	String str;
};

#endif