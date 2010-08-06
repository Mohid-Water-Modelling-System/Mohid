#include "tokenlist.h"

TokenList::TokenList(const String str, const String separator)
{
	this->str = str;
	this->separator = separator;

	Split();
}

TokenList::~TokenList()
{
}

void TokenList::SetSeparator(const String separator)
{
	this->separator = separator;
	Split();
}

void TokenList::SetString(String &str)
{
	this->str = str;
	Split();
}

String *TokenList::NextToken()
{
	String *str;

	if (EOTL())
		return NULL;
	else
	{
		str = &(tokens[token_index]->name);
		token_index++;
	}

	return str;
}

int TokenList::NumberOfTokens()
{
	return tokens.Count();
}

int TokenList::TokenIndex()
{
	return token_index;
}

void TokenList::MoveTo(int token_index)
{
	this->token_index = token_index;
}

void TokenList::Reset()
{
	token_index = 0;
}

bool TokenList::EOTL()
{
	if (token_index >= tokens.Count())
		return true;

	return false;
}

void TokenList::Split()
{
	int p = 0;
	int next_sep;
	int sep_length = separator.Length();
	String t_token;

	tokens.Clear();

	while (p < str.Length())
	{
		next_sep = str.Find(separator);
		if (next_sep == -1)
			next_sep = str.Length();

		t_token = str(p, next_sep - p);
		tokens.AddNew(&t_token);

		p = next_sep + sep_length;
	}

	token_index = 0;
}


