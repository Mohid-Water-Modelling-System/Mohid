#include "datastructures.h"

//===========================
// Constructors & Destructors
//===========================
Reaction::Reaction() 
{ 
	token_list = NULL;
	token_list = new List<ReactionToken>; 

	Reset();
}

Reaction::Reaction(Reaction *copy)
{
	token_list = NULL;
	token_list = new List<ReactionToken>; 

	copy->CopyTo(this);
}

Reaction::~Reaction() 
{ 
	delete token_list; 
}

//===========================
// Public functions - List management
//===========================
void Reaction::Reset() 
{ 
	name = "";

	for (int i = 0; i < 8; i++)
		logk[i] = 0.0;
	
	token_list->Clear(); 
}

Reaction *Reaction::Copy()
{
	return new Reaction(this);
}

void Reaction::CopyTo(Reaction *copy)
{
	int i;

	copy->name = name;

	for (i = 0; i < 8; i++)
		copy->logk[i] = logk[i];

	token_list->CopyTo(copy->token_list);
}

//===========================
// Public Functions
//===========================
void Reaction::AddTRXN(Reaction *r, LDBLE coef, bool combine)
{
  int i;

	if (token_list->Count() == 0)
  {
		for (i = 0; i < 7; i++)
			logk[i] = r->logk[i];
	}
  else
  {
    for (i = 0; i < 7; i++)
      logk[i] += (coef * r->logk[i]);
  }

	ReactionToken *rt_ptr, *new_rt;
		
	for (i = 0; i < r->token_list->Count(); i++)
	{
		rt_ptr = r->token_list->Element(i);

		if (rt_ptr->s == NULL)
			break;

		new_rt = token_list->AddNew();

		new_rt->name = rt_ptr->s->name;
		new_rt->s = rt_ptr->s;
		new_rt->coef = coef * rt_ptr->coef;
	}

  if (combine)
    TRXNCombine();
}

void Reaction::AddPhaseTRXN(Reaction *r, LDBLE coef, bool combine)
{
  int i;

	if (token_list->Count() == 0)
  {
		for (i = 0; i < 8; i++)
			logk[i] = r->logk[i];
	}
  else
  {
    for (i = 0; i < 8; i++)
      logk[i] += (coef * r->logk[i]);
  }

	ReactionToken *rt_ptr, *new_rt;

	for (i = 0; i < r->token_list->Count(); i++)
	{
		rt_ptr = r->token_list->Element(i);
		new_rt = token_list->AddNew();

		if (rt_ptr->s != NULL)
		{
			new_rt->name = rt_ptr->s->name;
			new_rt->s = rt_ptr->s;
		}
		else
		{
			new_rt->name = rt_ptr->name;
			new_rt->s = NULL;
		}

		new_rt->coef = coef * rt_ptr->coef;
	}

  if (combine)
    TRXNCombine();
}

void Reaction::AddPhaseTRXN(Phase *p, Reaction *r) 
{
	char *ptr, token[DEFAULT_STRING_LENGTH];
  String name;
  LDBLE z;

	p->formula.Copy(token);
	ptr = token;
  GetSpeciesNameAndZ(&ptr, name, z);

	token_list->Clear();
	ReactionToken * new_rt = token_list->AddNew();
	new_rt->name = p->formula;
	new_rt->z = z;
  new_rt->s = NULL;
  new_rt->u = NULL;
	
  new_rt->coef = p->rxn->token_list->Element(0)->coef;

	ReactionToken *rt_ptr;
	for (int i = 1; i < r->token_list->Count(); i++)
  {
		rt_ptr = r->token_list->Element(i);
		new_rt = token_list->AddNew();

    new_rt->name = rt_ptr->s->name;
    new_rt->z = rt_ptr->s->z;
    new_rt->s = NULL;
    new_rt->u = NULL;
    new_rt->coef = rt_ptr->coef;
  }
}

void Reaction::AddSpeciesTRXN(Species *s)
{
	token_list->Clear();

	ReactionToken *rt_ptr;
	ReactionToken *new_rt;
  for (int i = 0; i < s->rxn->token_list->Count(); i++)
  {
		rt_ptr = s->rxn->token_list->Element(i);
		new_rt = token_list->AddNew();

    new_rt->name = rt_ptr->s->name;
    new_rt->z = rt_ptr->s->z;
    new_rt->s = rt_ptr->s;
    new_rt->u = NULL;
    new_rt->coef = rt_ptr->coef;
  }
}

void Reaction::TRXNCombine()
{
	//
	// Combines coefficients of tokens that are equal in temporary
	// reaction structure, trxn.
	//

  int j, k;

	// Sort trxn species
	token_list->Sort(1, false);

	// Combine trxn tokens
	j = 1;
  for (k = 2; k < token_list->Count(); k++)
  {
		if ((*token_list)[k]->s != NULL)
    {
      if (j > 0 && (*token_list)[k]->s == (*token_list)[j]->s)
      {
				(*token_list)[j]->coef += (*token_list)[k]->coef;
				
				if (Equal((*token_list)[j]->coef, (LDBLE)0.0, (LDBLE)1e-5))
					j--;
      }
      else
      {
				j++;
				
				if (k != j)
				{
					(*token_list)[j]->name = (*token_list)[k]->name;
					(*token_list)[j]->s = (*token_list)[k]->s;
					(*token_list)[j]->coef = (*token_list)[k]->coef;
				}
      }
    }
    else
    {
      if (j > 0 && (*token_list)[k]->s == (*token_list)[j]->s && (*token_list)[k]->name == (*token_list)[j]->name)
      {
				(*token_list)[j]->coef += (*token_list)[k]->coef;

				if (Equal((*token_list)[j]->coef, (LDBLE)0.0, (LDBLE)1e-5))
					j--;
      }
      else
      {
				j++;

				if (k != j)
				{
					(*token_list)[j]->name = (*token_list)[k]->name;
					(*token_list)[j]->s = (*token_list)[k]->s;
					(*token_list)[j]->coef = (*token_list)[k]->coef;
				}
      }
    }
  }

	//token_list->SetNewCapacity(j + 1);		
	token_list->CompressByIndex(j + 1);
}

void Reaction::TRXNReverseK()
{
	//
	// Changes K from dissociation to association and back
	//

	// Accumulate log k for reaction
  for (int i = 0; i < 8; i++)
    logk[i] = -logk[i];
}

bool Reaction::TRXNSwap(String token, int position)
{
	//
	// Moves specified token to initial position in reaction.
	// Input: token, token name to move to position (position = 0 by default)
	// Return: ERROR, if token not found.
	//         OK, if token moved to initial position.
	//

  int index;

	if (position < 0 || position >= token_list->Count())
		return false;

	// Locate token
	if (token_list->Search(&token, index) == NULL)
		return false;

	// Swap token to first position
	token_list->Swap(position, index);

	// Make coefficient of token -1.0
	LDBLE coef = (LDBLE)-1.0 / token_list->Element(0)->coef;

  return TRXNMultiply(coef);
}

bool Reaction::TRXNMultiply(LDBLE coef)
{
	//
	// Multiplies temporary reaction, trxn,  by a constant
	// Arguments:
	//    input: coef, multiplier.
	//

  int i;

	// Multiply log k for reaction
	for (i = 0; i < 7; i++)
    logk[i] *= coef;

	// Multiply coefficients of reaction
	for (i = 0; i < token_list->Count(); i++)
    (*token_list)[i]->coef *= coef;

	return true;
}

//===========================
// Public Functions - Copy's
//===========================
void Reaction::CopyReactions(Reaction *reaction)
{
  int i;

	for (i = 0; i < 8; i++)
    reaction->logk[i] = logk[i];

	ReactionToken *rt_new, *rt_ptr;
	for (i = 0; i < token_list->Count(); i++)
  {
		rt_ptr = token_list->Element(i);
		rt_new = reaction->token_list->AddNew();

    rt_new->s = rt_ptr->s;
    rt_new->name = rt_ptr->name;
    rt_new->coef = rt_ptr->coef;
  }
}

void Reaction::CopyReactionsToSpecies(Reaction *reaction)
{
	int i;

	reaction->token_list->Clear();
	reaction->token_list->SetNewCapacity(token_list->Count());

	ReactionToken *rt_new, *rt_ptr;	
	for (i = 0; i < token_list->Count(); i++)
	{
		rt_ptr = token_list->Element(i);
		rt_new = reaction->token_list->AddNew();

		rt_new->name = rt_ptr->name;
		rt_new->s = rt_ptr->s;
		rt_new->coef = rt_ptr->coef;
	}

	for (i = 0; i < 8; i++)
		reaction->logk[i] = logk[i];
}

void Reaction::CopyReactionsToExchange(class Reaction *reaction)
{
	reaction->token_list->Clear();
	reaction->token_list->SetNewCapacity(token_list->Count());

	int i;

	ReactionToken *rt_new, *rt_ptr;
	for (i = 0; i < token_list->Count(); i++)
	{
		rt_ptr = token_list->Element(i);
		rt_new = reaction->token_list->AddNew();

		rt_new->s = rt_ptr->s;
		rt_new->coef = rt_ptr->coef;
	}
}

void Reaction::CopyReactionsToPhase(Reaction *reaction)
{
	reaction->token_list->SetNewCapacity(token_list->Count());

	ReactionToken *rt_new;

	rt_new = reaction->token_list->AddNew();
	rt_new->coef = token_list->Element(0)->coef;
  rt_new->s = token_list->Element(1)->s;
	rt_new->name = token_list->Element(1)->name;

	ReactionToken *rt_ptr;
	for (int i = 1; i < token_list->Count(); i++)
  {
		rt_ptr = token_list->Element(i);
		rt_new = reaction->token_list->AddNew();
		
		if (rt_ptr->s == NULL)
			rt_new->name = rt_ptr->name;
		else
			rt_new->name = "";

		rt_new->s = rt_ptr->s;
		rt_new->coef = rt_ptr->coef;
  }
}

//===========================
// Public Functions - Others
//===========================

void Reaction::PrintToFile(FILE *file, int spaces)
{
	/*
	String spaces_str(spaces, " ");
	int i;

	fprintf(file, "%sname: %s\n", spaces_str.CharPtr(), name.CharPtr());
	fprintf(file, "%s  logk: ", spaces_str.CharPtr());
	for (i = 0; i < 8; i++)
		fprintf(file, "%d:(%f)  ",  i, logk[i]);  
	fprintf(file, "\n%s  Tokens: %d\n", spaces_str.CharPtr(), token_list->Count());
	for (i = 0; i < token_list->Count(); i++)
	{
		fprintf(file, "%s   %d:\n", spaces_str.CharPtr(), i);
		(*token_list)[i]->PrintToFile(file, spaces + 5); 
	}
	*/
}