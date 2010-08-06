#ifndef keywordH
#define keywordH

#include "string.h"

class Keyword
{
public:
	Keyword();
	Keyword(Keyword &copy);
	Keyword(String &keyword, int id);
	~Keyword();

public:
	void Reset();
	Keyword *Copy();
	void CopyTo(Keyword &copy);

public:
	String name;
	int id;

};

#endif