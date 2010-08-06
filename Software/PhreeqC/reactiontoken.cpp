#include "datastructures.h"

ReactionToken::ReactionToken()
{
	Reset();
}

ReactionToken::ReactionToken(ReactionToken *copy)
{
	name = copy->name;

	coef = copy->coef;
	z    = copy->z;

	s = copy->s;
	u = copy->u;
}

ReactionToken::ReactionToken(Species *s, const LDBLE coef)
{
	name = "";

	this->coef = coef;
	this->z    = 0;

	this->s = s;
	this->u = NULL;
}

ReactionToken::~ReactionToken()
{
}

void ReactionToken::Reset()
{
	name = "";

	coef = 0.0;
	z    = 0.0;

	s = NULL;
	u = NULL;
}

/*
ReactionToken *ReactionToken::Copy()
{
	return new ReactionToken(this);
}
*/
void ReactionToken::CopyTo(ReactionToken *copy)
{
	copy->name = name;

	copy->coef = coef;
	copy->z    = z;

	copy->s = s;
	copy->u = u;
}

void ReactionToken::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%sname: %s\n", spaces_str.CharPtr(), name.CharPtr());
	fprintf(file, "%s  coef: %f\n", spaces_str.CharPtr(), coef);  
	fprintf(file, "%s  z: %f\n", spaces_str.CharPtr(), z);  
	if (u != NULL)
		fprintf(file, "%s  u: %s (%p)\n", spaces_str.CharPtr(), u->name.CharPtr(), u);
	else
		fprintf(file, "%s  u: NULL\n", spaces_str.CharPtr());  
}