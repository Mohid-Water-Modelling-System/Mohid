#ifndef checkdatabaseH
#define checkdatabaseH

#include "datastructures.h"

class Check
{
public:
	Check(GlobalData *gd);

public:
  bool CheckSpecies();
	bool CheckMaster();
	bool CheckElement();
	bool CheckPhase();

public:
	void SelectLOGKExpression(LDBLE *source_k, LDBLE *target_k);
	bool AddOtherLogk(LDBLE *source_k, List<NameCoef> *add_logk);

private:
	GlobalData *gd;
};

#endif