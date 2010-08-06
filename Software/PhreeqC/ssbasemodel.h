#ifndef ssbasemodelH
#define ssbasemodelH

#include "datastructures.h"
#include "modeldata.h"
#include "parseengine.h"

class SSBaseModel:
	public ParseEngine
{
public:
	SSBaseModel(GlobalData *gd, ModelData *md);
	~SSBaseModel();

	LDBLE KCalc(LDBLE * logk, LDBLE tempk);


public:
	bool SSPrep(LDBLE t, SS *s_s_ptr);

protected:
	LDBLE a0, a1, kc, kb;

private:
	bool Scan(LDBLE * xx0, LDBLE * xx1);
	LDBLE FSpinodal(LDBLE x);
	LDBLE Halve(LDBLE x0, LDBLE x1, LDBLE tol);
	bool SLNQ(int n, LDBLE * a, LDBLE * delta, int ncols);
	bool SolveMisc(LDBLE * xxc1, LDBLE * xxc2, LDBLE tol);

private:
	GlobalData *gd;
	ModelData *md;

};

#endif