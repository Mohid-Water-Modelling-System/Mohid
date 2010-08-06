#include "databaseengine.h"
#include <math.h>
//#include "checkdatabase.h"
#include <conio.h>

//===========================
// Constructors & Destructors
//===========================
DBEngine::DBEngine(GlobalData *gd):
Tidy(gd, NULL)
{
	this->gd = gd;
}

DBEngine::~DBEngine()
{
}

//===========================
// Public Functions - Line read group
//===========================
bool DBEngine::GetValueFromNextToken(LDBLE &value)
{
	if (!line.NextToken())
		return false;

	value = (LDBLE)line.Token().AsDouble();
	return true;
}

bool DBEngine::ReadDeltaH(LDBLE &delta_h, DELTA_H_UNIT &units)
{
	bool kilo, joul;
  String token;

	//Read delta H
 
	delta_h = 0.0;
	if (!line.NextToken())
		return false;

	if (line.TokenType() == _empty_ || line.TokenType() == _error_)
    return false;
  
	delta_h = (LDBLE)line.Token().AsDouble();

	// Read delta H units
 
	if (!line.NextToken())
		return true;
	
	if (line.TokenType() == _empty_ || line.TokenType() == _error_)
    return true;

	units = _kjoules_;
  kilo = true;
  joul = true;

	if (line.TokenType() == _upper_ || line.TokenType() == _lower_)
  {
		token = line.Token();
		token.ToLower();

		if (token[0] != 'k')
    {
      // convert to kilo
      kilo = false;
      delta_h /= 1000.0;
    }
		if (token.Find("c") != -1)
    {
      // convert to joules 
      delta_h *= (LDBLE)JOULES_PER_CALORIE;
      joul = false;
    }
  }

  if (kilo && joul)
    units = _kjoules_;
  else if (!kilo && joul)
    units = _joules_;
  else if (kilo && !joul)
    units = _kcal_;
  else
    units = _cal_;
  
  return true;
}

bool DBEngine::ReadAnalyticalExpression(LDBLE *logk)
{
	int j, i;
	for (j = 0, i = 0; j < 5; j++)
  {
		logk[j] = 0.0;
		if (line.NextToken())
		{
			i++;
			logk[j] = (LDBLE)line.Token().AsDouble();
		}
	}

	if (i == 0)
		return false;

	return true;
}

bool DBEngine::ReadAddLogK(List<NameCoef> *list)
{
	NameCoef *item = list->AddNew();
	
	if (!line.NextToken() || line.TokenType() == _empty_)
		return false;

	// read name
	item->name = line.Token();

	// read coef
	if (!line.NextToken() || line.TokenType() == _empty_)
		item->coef = 1;
	else
		item->coef = (LDBLE)line.Token().AsDouble();

	return true;
}

bool DBEngine::ReadAddConstant(List<NameCoef> *list)
{
	NameCoef *item = list->AddNew();
	
	if (!line.NextToken() || line.TokenType() == _empty_)
		return false;

	item->name = "XconstantX";
	item->coef = (LDBLE)line.Token().AsDouble();

	return true;
}

bool DBEngine::ReadMassBalance(Species *s)
{
	if (!line.NextToken())
		return false;

	s->mole_balance = line.Token();

	if(!LoadSecondarySpecies(line.Token(), s))
		return false;

	return true;
}



//===========================
// Public Functions - Support group
//===========================


//===========================
// Protected Functions
//===========================
bool DBEngine::LoadSecondarySpecies(String str, Species *s)
{
	eos_list->Clear();

	char str_ptr[DEFAULT_STR_LENGTH], *ptr;
	str.Copy(str_ptr);
	ptr = &str_ptr[0];
	LDBLE coef = 1.0;
	parent_count = 0;
	if(!GetSecondaryElementsInSpecies(&ptr, coef))
		return false;

	eos_list->CopyTo(s->e_sec_list);
	return true;
}


bool DBEngine::AddPSIMasterSpecies(String &str)
{
  Master *m;
	ReactionToken *t;
  String token;
  int plane;
	char str_[DEFAULT_STR_LENGTH], *ptr;
	LDBLE coef = 1.0;

  for (plane = SURF_PSI; plane <= SURF_PSI2; plane++)
  {
    token = str;
    switch (plane)
    {
    case SURF_PSI:
      break;
    case SURF_PSI1:
      token += "b";
      break;
    case SURF_PSI2:
      token += "d";
      break;
    }

    m = gd->master_list.Search(&token);
    if (m == NULL)
    {
      m = gd->master_list.Store(&token);
      m->type = plane;
      m->e = gd->element_list.Store(&token, false);

      m->s = gd->species_list.Search(&token);
      if (m->s == NULL)
      {
				m->s = gd->species_list.Store(&token, false);
				m->s->z = 0.0;
      }

			eos_list->Clear();			
			parent_count = 0;
			token.Copy(str_);
			ptr = &str[0];
			if (!GetElementsInSpecies(&ptr, coef))
				return false;
			eos_list->CopyTo(m->s->eos_list);
      m->s->type = plane;
      m->primary = true;
      
			// Define reaction for psi
			t = m->s->rxn->token_list->AddNew();
			t->s = m->s;
      t->coef = -1.0;

			t = m->s->rxn->token_list->AddNew();
			t->s = m->s;
      t->coef = 1.0;
    }
  }

  return true;
}

