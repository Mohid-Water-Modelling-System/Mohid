#include "datastructures.h"

STCoef::STCoef()
{
	Reset();
}

STCoef::STCoef(STCoef *copy)
{
	copy->CopyTo(this);
}

STCoef::~STCoef()
{
}

void STCoef::Reset()
{
	source = NULL;
	target = NULL;
	coef = 1.0;
}

void STCoef::CopyTo(STCoef *copy)
{
	copy->source = source;
	copy->target = target;
	copy->coef = coef;
}

void STCoef::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s source:%s\n", spaces_str.CharPtr(), *source);
	fprintf(file, "%s target:%s\n", spaces_str.CharPtr(), *target);
	fprintf(file, "%s coef:%s\n", spaces_str.CharPtr(), coef);
}