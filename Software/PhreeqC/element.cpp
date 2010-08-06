#include "datastructures.h"

Element::Element()
{
	Reset();
}

Element::Element(Element *copy)
{
	copy->CopyTo(this);
}

Element::~Element()
{
}

void Element::Clear(void)
{
}

/*
Element *Element::Copy()
{
	return new Element(this);
}
*/
void Element::CopyTo(Element *copy)
{
	copy->name = name;

	copy->master  = master;
	copy->primary = primary;

	copy->gfw = gfw;
}

void Element::PrintToFile(FILE *file, int spaces)
{
	/*
	String spaces_str(spaces, " ");

	fprintf(file, "%s%s\n", spaces_str.CharPtr(), name.CharPtr());
	fprintf(file, "%s  gfw: %f", spaces_str.CharPtr(), gfw);
	if (master != NULL)
		fprintf(file, "%s  master: %s (%p)\n", spaces_str.CharPtr(), master->name.CharPtr(), master);
	else
		fprintf(file, "%s  master: NULL\n", spaces_str.CharPtr());  
	if (primary != NULL)
		fprintf(file, "%s  primary: %s (%p)\n", spaces_str.CharPtr(), primary->name.CharPtr(), primary);
	else
		fprintf(file, "%s  primary: NULL\n", spaces_str.CharPtr());  
		*/
}

void Element::Reset(void)
{
	name = "";

	master  = NULL;
	primary = NULL;

	gfw = 0.0;
}