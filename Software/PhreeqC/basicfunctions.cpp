#include "basicfunctions.h"
#include <math.h>

//===========================
// Public Functions
//===========================
bool BasicFunctions::GetNumber(char **t_ptr, LDBLE &num)
{
	bool decimal = false;
  char c;
  String token;

  num = 1.0;
  c = **t_ptr;

  if (isdigit ((int) c) || (c == '.'))
  {
    while (isdigit ((int) c) || (c == '.'))
    {
      if (c == '.' && !decimal)
				decimal = true;
			else if (c == '.' && decimal)
				break;

      token += c;
      c = *(++(*t_ptr));
    }

		try
		{
			num = (LDBLE)token.AsDouble();
		}
		catch(...)
		{
			return false;
		}
  }

  return true;
}

bool BasicFunctions::GetElement(char **str, String &element)
{
  char c;

  c = *(*str)++;
  if (c == '\0')
		return false;

	// Load name into char array element
  element = c;
  if (c == '[')
  {
    while ((c = (**str)) != ']')
    {
      element += c;
      (*str)++;
      if ((c = (**str)) == ']')
      {
				element += c;
				(*str)++;
				break;
      }
      else if (**str == '\0')
				return false;
    }
    while (islower ((int) (c = (**str))) || c == '_')
    {
      element += c;
      (*str)++;
    }
  }
  else
  {
    while (islower ((int) (c = (**str))) || c == '_')
    {
      element += c;
      (*str)++;
    }
  }
  
  return true;
}

void BasicFunctions::CorrectPlusOne(String &to_correct)
{
	if (to_correct.Length() >= 3)
		if (to_correct[to_correct.Length() - 1] == '1' && to_correct.IsAmong(to_correct.Length() - 2, "+-"))
			to_correct = to_correct(0, to_correct.Length() - 1);
}

bool BasicFunctions::Equal(const LDBLE a, const LDBLE b, const LDBLE tol)
{
  if (fabs (a - b) <= tol)
    return true;

  return false;
}

bool BasicFunctions::GetSecondary(char **str, String &element)
{
	//
	//      Function reads an element name out of the equation string.
	//      An element name is composed of a capital letter followed by any number
	//      of lower case characters.
	//
	//      Arguments:
	//         str   input, points to position in the equation to begin
	//               output, points to next character of equation after
	//                      element name.
	//         element  input pointer to place to return element character string
	//

  int j;
  char c;
  char *ptr;

  c = *(*str)++;
  if (c == '\0')
    return false;

	// Load name into char array element
  element = c;
  if (c == '[')
  {
    while ((c = (**str)) != ']')
    {
      element += c;
      (*str)++;
      if ((c = (**str)) == ']')
      {
				element += c;
				(*str)++;
				c = (**str);
				break;
      }
      else if ((c = (**str)) == '\0')
				return false;
    }
    while (islower ((int) (c = (**str))) || c == '_')
    {
      element += c;
      (*str)++;
    }
  }
  else
  {
    while (islower ((int) (c = (**str))) || c == '_')
    {
      element += c;
      (*str)++;
    }
  }

	// Check if secondary master species element
  j = element.Length();
  ptr = *str;
  if (c == '(')
  {
    // copy parenthesis
    element += c;
    (*str)++;

    // copy number
    for (;;)
    {
      c = **str;
      if (isdigit ((int) c) || c == '-' || c == '.')
      {
				element += c;
				(*str)++;
      }
      else if (c == '+')
				(*str)++;
      else
				break;
    }

    // go back to before parenthesis
    if (c != ')')
    {
      element.Cut(j);
      *str = ptr;
      // put in closing parenthesis
    }
    else
    {
      element += c;
      (*str)++;
    }
  }

	return true;
}

bool BasicFunctions::GetSpeciesNameAndZ(char **str, String &name, LDBLE &z)
{
  int i, j;
  char c;
  char *ptr, *ptr1;
  String charge;

  ptr = *str;
  i = 0;
	name   = "";
	charge = "";

  // Find end of token or begining of charge
  while (((c = *ptr) != '+') && (c != '-') && (c != '=') && (c != '\0'))
  {
    name += c;
    if (c == '[')
    {
      ptr++;
      while ((c = *ptr) != ']')
      {
				if (c == '\0')
					return false;
				
				name += c;
				ptr++;
      }
      name += c;
    }
    ptr++;
  }

	// Check for an empty string
	if (name.Length() == 0)
    return false;

	// End of token is = or \0, charge is zero
  if (c == '=' || c == '\0')
  {
    *str = ptr;
    z = 0.0;
  }
  else
  {
		// Copy characters into charge until next species or end is detected
    ptr1 = ptr;
    while ((isalpha ((int) (c = *ptr1)) == 0) && (c != '(') && (c != ')') && (c != ']') && (c != '[') && (c != '=') && (c != '\0'))
    {
      charge += c;
      ptr1++;
    }
		
		// Go back to last + or - if not end of side,
		// everything before the last + or - in charge is part of the charge
    if ((c != '=') && (c != '\0'))
    {
	    j = 0;
      while (((c = *ptr1) != '+') && (c != '-'))
      {
				j++;
				ptr1--;
      }

			charge = charge(0, charge.Length() - j);
    }
    *str = ptr1;
		
		// Charge has been written, now need to check if charge has legal format
		if (!GetCharge(charge, z))
			return false;

    name += charge;
  }

	return true;
}

bool BasicFunctions::GetZ(String &token, LDBLE &z)
{
	long minus, plus;

	token.Trim(_all_);

	plus  = token.Find("+");
	minus = token.Find("-");

	if (plus != -1)
	{
		z = (LDBLE)token(plus, -1).AsDouble();
		if (z == 0.0)
			z = 1.0;
	}
	else if (minus != -1)
	{
		z = (LDBLE)token(plus, -1).AsDouble();
		if (z == 0.0)
			z = -1.0;
	}
	else
		z = 0.0;

	return true;
}




//===========================
// Protected Functions
//===========================
bool BasicFunctions::GetCharge(String &charge, LDBLE &z)
{
  int i;
  char c, c1;

	// Charge is zero
	if (charge.IsEmpty() || (c = charge[0]) == '\0')
  {
    z = 0.0;
    return true;
  }

	// Error check for + or - at start of string
  if (c != '+' && c != '-')
    return false;
  
	// Count string of +'s or -'s
  i = 0;
  while (c == (c1 = charge[i++]));
  i--;

  if (c1 == '\0')
  {
    if (c == '-')
      i = -i;
  }
  else
  {
		// + or - followed by a number
		if (charge.Find(".") == -1)
			i = charge.AsInt();
		else
		{
			charge.ClearRightZeros();
			z = (LDBLE)charge.AsDouble();
			return true;
		}
  }

	// Charge is zero, must have had +0 or -0 in eqn
  if (i == 0)
    charge = "";
 
	// Charge is +1 or -1, single + or -
  if (abs(i) == 1)
    charge = c;

	// Abs(z)>1, set charge to + or - plus integer
  if (abs(i) > 1)
		if (i > 0)
		{
			charge = "+";
			charge += i;
		}
		else
			charge = i;

	z = (LDBLE)i;
	return true;
}

