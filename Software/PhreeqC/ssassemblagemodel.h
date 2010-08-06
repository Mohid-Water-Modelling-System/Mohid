#ifndef ssassemblagemodelH
#define ssassemblagemodelH

#include "modeldata.h"
#include "ssbasemodel.h"

class SSAssemblageModel:
	public SSBaseModel
{
public:
	SSAssemblageModel(GlobalData *gd, ModelData *md);
	~SSAssemblageModel();

public:
	bool TidySSAssemblage();

protected:
	bool SSCalcA0A1(SS *ss);

public:
	SSAssemblage *ssa;

private:
	GlobalData *gd;
	ModelData *md;
};

#endif