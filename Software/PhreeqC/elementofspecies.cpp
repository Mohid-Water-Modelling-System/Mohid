#include "datastructures.h"

ElementOfSpecies::ElementOfSpecies()
{
	Reset();
}

ElementOfSpecies::ElementOfSpecies(ElementOfSpecies *copy)
{
	copy->CopyTo(this);
}

ElementOfSpecies::~ElementOfSpecies()
{
}

void ElementOfSpecies::Reset(void)
{
	name = "";

	e = NULL; //never delete

	coef = 0.0;
}

/*
ElementOfSpecies *ElementOfSpecies::Copy()
{
	return new ElementOfSpecies(this);
}
*/

void ElementOfSpecies::CopyTo(ElementOfSpecies *copy)
{
	copy->name = name;
	
	copy->e = e;

	copy->coef = coef;
}

void ElementOfSpecies::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());
	fprintf(file, "%s  coef: %f\n", spaces_str.CharPtr(), coef);
	if (e != NULL)
		fprintf(file, "%s  e: %s (%p)\n", spaces_str.CharPtr(), e->name, e);
}