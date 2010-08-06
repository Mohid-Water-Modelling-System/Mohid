#include "keyword.h"

Keyword::Keyword()
{
	Reset();
}

Keyword::Keyword(Keyword &copy)
{
	copy.CopyTo(*this);
}

Keyword::Keyword(String &name, int id)
{
	this->name = name;
	this->id = id;
}

Keyword::~Keyword()
{
}

void Keyword::Reset()
{
	name = "";
	id = -1;
}

Keyword *Keyword::Copy()
{
	return new Keyword(*this);
}

void Keyword::CopyTo(Keyword &copy)
{
	copy.name = name;
	copy.id = id;
}

