#include "check.h"

Check::Check(GlobalData *gd)
{
	this->gd = gd;
}

bool Check::CheckSpecies()
{
	Species *s;
	for (int i = 0; i < gd->species_list.Count(); i++)
  {
		s = gd->species_list[i];

		if (s->eos_list->IsEmpty())
			return false;

		if (s->rxn->token_list->Count() <= 0)
			return false;
    else
    {
			SelectLOGKExpression(s->logk, s->rxn->logk);
      //AddOtherLogk(s->rxn->logk, s->add_logk);
    }
  }

  return true;
}

bool Check::CheckMaster()
{
	return true;
}

bool Check::CheckElement()
{
	return true;
}

bool Check::CheckPhase()
{
	return true;
}

void Check::SelectLOGKExpression(LDBLE *source_k, LDBLE *target_k)
{
  int j;
	bool analytic;

  analytic = false;
  for (j = 2; j < 7; j++)
  {
    if (source_k[j] != 0.0)
    {
      analytic = true;
      break;
    }
  }
  if (analytic == true)
  {
    target_k[0] = 0.0;
    target_k[1] = 0.0;
    for (j = 2; j < 7; j++)
    {
      target_k[j] = source_k[j];
    }
  }
  else
  {
    target_k[0] = source_k[0];
    target_k[1] = source_k[1];
    for (j = 2; j < 7; j++)
    {
      target_k[j] = 0.0;
    }
  }
}

bool Check::AddOtherLogk(LDBLE *source_k, List<NameCoef> *add_logk)
{
  if (add_logk->Count() <= 0)
    return true;

	return true;
}