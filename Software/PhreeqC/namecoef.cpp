#include "datastructures.h"

NameCoef::NameCoef()
{
	Reset();
}

NameCoef::NameCoef(NameCoef *copy)
{
	copy->CopyTo(this);
}

NameCoef::~NameCoef()
{
}

void NameCoef::Reset()
{
	name = "";
	coef = 0.0;
}

/*
NameCoef *NameCoef::Copy()
{
	return new NameCoef(*this);
}
*/

void NameCoef::CopyTo(NameCoef *copy)
{
	copy->name = name;
	copy->coef = coef;
}

void NameCoef::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());
	fprintf(file, "%s  coef: %f\n", spaces_str.CharPtr(), coef);
}