#include "datastructures.h"

UnknownInfo::UnknownInfo()
{
	Reset();
}

UnknownInfo::UnknownInfo(UnknownInfo *copy)
{
	copy->CopyTo(this);
}

UnknownInfo::~UnknownInfo()
{
}

void UnknownInfo::CopyTo(UnknownInfo *copy)
{
	copy->u = u;

	copy->source = source;
	copy->gamma_source = gamma_source;
	copy->coef = coef;
}

void UnknownInfo::Reset()
{
	u = NULL;

	source = NULL;
	gamma_source = NULL;
	coef = 0.0;
}

void UnknownInfo::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%s  coef: %f\n", spaces_str.CharPtr(), coef);
	fprintf(file, "%s  source: %f (%p)\n", spaces_str.CharPtr(), *source, source);
	fprintf(file, "%s  gamma_source: %f (%p)\n", spaces_str.CharPtr(), *gamma_source, gamma_source);
	
	if (u != NULL)
	{
		fprintf(file, "%s  u: %s (%p)\n", spaces_str.CharPtr(), u->name.CharPtr(), u);
		u->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  u: NULL\n");
}