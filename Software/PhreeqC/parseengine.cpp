#include "parseengine.h"


ParseEngine::ParseEngine(GlobalData *gd)
{
	eos_list = NULL;
	t_eos_list = NULL;

	eos_list = new List<ElementOfSpecies>;
	t_eos_list = new List<ElementOfSpecies>;

	this->gd = gd;
}

ParseEngine::~ParseEngine()
{
	delete eos_list;
	delete t_eos_list;
}
//===========================
// Public Functions
//===========================
void ParseEngine::SaveEOSList(List<class ElementOfSpecies> *list)
{
	if (eos_list->Count() > 0)
	{
		eos_list->Sort(0, false);

		CombineElements();

		eos_list->CopyTo(list);
	}
	else
		list->Clear();	
}

bool ParseEngine::GetSpecies(char **str)
{
	ReactionToken *rt_new = rt.token_list->AddNew();

	if (GetCoef(str, rt_new->coef))
		if (GetSpeciesNameAndZ(str, rt_new->name, rt_new->z))
			return true;
	
	return false;
}

bool ParseEngine::ParseEq(String &eq, bool association)
{
	rt.Reset();
	eq.RemoveSpaces();

	char str_[DEFAULT_STR_LENGTH], *ptr;

	long p = eq.Find("=");
	if (p == -1)
		return false;

	String ls(eq(0, p));
	String rs(eq(p + 1, -1));

	ls.Copy(str_);
	ptr = &str_[0];
	while (*ptr != '\0')
	{
		if (!GetSpecies(&ptr))
			return false;

		if (!association)
			rt.token_list->Last()->coef *= -1;
	}

	rs.Copy(str_);
	ptr = &str_[0];
	if (association)
	{
		if (!GetSpecies(&ptr))
			return false;

		rt.token_list->Last()->coef *= -1;
		rt.token_list->Swap(0, rt.token_list->Count() - 1);
	}

	while (*ptr != '\0')
	{
		if (!GetSpecies(&ptr))
			return false;

		if (association)
			rt.token_list->Last()->coef *= -1;
	}	

	//printf("\n%s\n", eq.CharPtr());
	//getch();
	rt.token_list->Sort(1, false);

	ReactionToken *rxnt_token = rt.token_list->Element(0);
	String name = rxnt_token->name;

	name.Replace("(s)", "");
	name.Replace("(S)", "");
	name.Replace("(g)", "");
	name.Replace("(G)", "");

	parent_count = 0;
	name.Copy(str_);
	ptr = &str_[0];
	if(!GetElementsInSpecies(&ptr, rxnt_token->coef))
		return false;

  if (!CombineElements())
		return false;

	for (int i = 0; i < eos_list->Count(); i++)
		(*eos_list)[i]->coef = -(*eos_list)[i]->coef;

	return true;
}

bool ParseEngine::GetElementsInSpecies(char **str, LDBLE coef)
{
	//
	// Makes a list of elements with their coefficients, stores elements
  // in eos_list.  
	// Also uses global variable parent_count.
	//
	// Arguments:
	// **t_ptr    input, point in token string to start looking
	//            output, is next position to start looking
	//   coef     input, coefficient to multiply subscripts by
	//

  int i, count;
  char c, c1;
  LDBLE d;
  String element;

  while (((c = **str) != '+') && (c != '-') && (c != '\0'))
  {
    // close parenthesis
    if (c == ')')
    {
      if (--parent_count < 0)	
				return false;

			(*str)++;
      return true;
    }

    c1 = *((*str) + 1);

    // beginning of element name
    if (isupper ((int) c) || (c == 'e' && c1 == '-') || (c == '['))
    {

			// Get new element and subscript
      if (!GetElement(str, element))
      	return false;
      
			ElementOfSpecies * elt_ptr = AddElementToList(element);

      if (!GetNumber (str, d))
				return false;

			elt_ptr->coef = d * coef;
			
			continue;
    }

		// Open parentheses
    if (c == '(')
    {
      count = eos_list->Count();
      parent_count++;
      (*str)++;

      if (!GetElementsInSpecies(str, coef))
				return false;

      if (!GetNumber(str, d))
				return false;

      for (i = count; i < eos_list->Count(); i++)
				(*eos_list)[i]->coef *= d;			

			continue;
    }

		// Colon
    if (c == ':')
    {
      count = eos_list->Count();
      (*str)++;
      
      if (!GetNumber(str, d))
				return false;

      if (!GetElementsInSpecies(str, coef))
				return false;

      for (i = count; i < eos_list->Count(); i++)
				(*eos_list)[i]->coef *= d;			

      continue;
    }

		return false;
  }

  if (parent_count != 0)
		return false;
  
	return true;
}

bool ParseEngine::GetSecondaryElementsInSpecies(char **str, LDBLE &coef)
{
  int i, count;
  char c, c1;
  LDBLE d;
  String element("");

  while (((c = **str) != '+') && (c != '-') && (c != '\0'))
  {
    // close parenthesis
    if (c == ')')
    {
      parent_count--;
      if (parent_count < 0)
				return false;

      (*str)++;
      return true;
    }

    c1 = *((*str) + 1);
    // beginning of element name
    if (isupper ((int) c) || c == '[' || (c == 'e' && c1 == '-'))
    {
			// Get new element and subscript
      if (!GetSecondary(str, element))
				return false;
      
      ElementOfSpecies * elt_ptr = AddElementToList(element);

      if (!GetNumber(str, d))
      	return false;
      
      elt_ptr->coef = d * coef;

			continue;
    }
 
		// Open parentheses
    if (c == '(')
    {
      count = eos_list->Count();
			/*
      if (c1 == ')')
      {
				sprintf (error_string, "Empty parentheses.");
				warning_msg (error_string);
      }
			*/
      parent_count++;
      (*str)++;
      if (!GetSecondaryElementsInSpecies(str, coef))
				return false;

			if (!GetNumber(str, d))
				return false;
      
      for (i = count; i < eos_list->Count(); i++)
				(*eos_list)[i]->coef *= d;d;

			continue;
    }

		// Colon
    if (c == ':')
    {
      count = eos_list->Count();
      (*str)++;
			if (!GetNumber(str, d))
				return false;

      if (!GetSecondaryElementsInSpecies(str, coef))
				return false;

			for (i = count; i < eos_list->Count(); i++)
				(*eos_list)[i]->coef *= d;			

			continue;
    }

		// Not beginning of element and not opening parentesis
    return false;
  }

  if (parent_count != 0)
    return false;
  
  return true;
}


//===========================
// Protected Functions
//===========================
bool ParseEngine::GetCoef(char **str, LDBLE &coef)
{
  int i;
  char c, c1;
  char *ptr, *ptr1;
  char token[DEFAULT_STR_LENGTH];;

  ptr = *str;		
  c = *ptr;			
  coef = 0.0;

	// No leading sign or number
  if (isalpha ((int) c) || (c == '(') || (c == ')') || (c == '[') || (c == ']'))
  {
    coef = 1.0;
    return true;
  }

	// Leading +, no digits
  c1 = *(ptr + 1);
  if (c == '+' && (isalpha ((int) c1) || (c1 == '(') || (c1 == ')') || (c1 == '[') || (c1 == ']')))
  {
    *str = ++ptr;
    coef = 1.0;
    return true;
  }

	// Leading -, no digits
  if (c == '-' && (isalpha ((int) c1) || (c1 == '(') || (c1 == ')') || (c1 == '[') || (c1 == ']')))
  {
    *str = ++ptr;
    coef = -1.0;
    return true;
  }
  
	i = 0;
	// Has number coefficient
  if (isdigit ((int) c) || c == '+' || c == '-' || c == '.')
  {
    while (isdigit ((int) c) || c == '+' || c == '-' || c == '.')
    {
      token[i++] = c;
      c = *(++ptr);
    }
    token[i] = '\0';
    *str = ptr;
    errno = 0;
    coef = (LDBLE)strtod(token, &ptr1);
    
		if ((errno == ERANGE) || (*ptr1 != '\0'))
			return false;
		
		return true;
  }

	return false;
}

bool ParseEngine::CombineElements(void)
{
  int i, j;

  if (eos_list->Count() < 1)
		return false;

  if (eos_list->Count() == 1)
		return true;

	eos_list->Sort(0, false);

	j = 0;
	ElementOfSpecies *a, *b;

	t_eos_list->Clear();
  for (i = 1; i < eos_list->Count(); i++)
  {
		b = (*eos_list)[i];
		a = (*eos_list)[j];

    if (b->e == a->e)
			a->coef += b->coef;
    else
    {
      j++;

      if (i != j)
      {
				a = (*eos_list)[j];
				a->name = b->e->name;
				a->e = b->e;
				a->coef = b->coef;
      }
    }
  }

	//eos_list->SetNewCapacity(j + 1);
	eos_list->CompressByIndex(j + 1);
  
  return true;
}

ElementOfSpecies *ParseEngine::AddElementToList(String &element)
{
	ElementOfSpecies * eos = eos_list->AddNew();

	eos->name = element;
	eos->e = gd->element_list.Store(&element, false);

	return eos;
}

bool ParseEngine::ComputeGFW(String &str, LDBLE &gfw)
{
	//
	// Input:  string contains a chemical formula
	// Output:  gfw contains the calculated gfw

  int i;
	char token[DEFAULT_STRING_LENGTH], *ptr;

	eos_list->Clear();
  parent_count = 0;

	str.Copy(token);
	ptr = &token[0];

	LDBLE coef = 1.0;
	if (!GetElementsInSpecies(&ptr, coef))
    return false;
  
  gfw = 0.0;

	ElementOfSpecies *eos_p;

  for (i = 0; i < eos_list->Count(); i++)
  {
		eos_p = (*eos_list)[i];

    if (eos_p->e->gfw <= 0.0)
      return false;

		gfw += eos_p->coef * eos_p->e->gfw;
  }

  return true;
}

void ParseEngine::CopyToTempEOSList(List<class ElementOfSpecies> *list, LDBLE coef)
{
	ElementOfSpecies *eos;

  if (list == NULL)
    return;

  for (int i = 0; i < list->Count(); i++)
  {
		eos = eos_list->AddNew();

		eos->e = (*list)[i]->e;
		eos->coef = (*list)[i]->coef * coef;
  }

  return;
}
