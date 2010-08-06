#ifndef tokenlistH
#define tokenlistH

#include "string.h"
#include "list.h"


class TokenListItem
{
public:
	String name;
};

class TokenList
{
public:
	TokenList(const String str = "", const String separator = " ");
	~TokenList();

public:
	void SetSeparator(const String separator = " ");
	void SetString(String &str);
	String *NextToken();
	int NumberOfTokens();
	int TokenIndex();
	void MoveTo(int token_index);
	void Reset();
	bool EOTL();

private:
	void Split();

private:
	String str;
	String separator;
	List<TokenListItem> tokens;
	int token_index;
};

#endif