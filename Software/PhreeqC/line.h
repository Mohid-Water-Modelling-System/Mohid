#ifndef lineH
#define lineH

#include "constants.h"
#include "list.h"
#include "keyword.h"

class Line
{
public:
	Line();
	~Line();

public:
	void Reset(void);
	void SetLine(String &text_to_set);
	bool NextToken(String &token_delimiter = String(" "));

public:
	RETURN_TYPE LineType(void) { return line_type; }
	RETURN_TYPE TokenType(void) { return token_type; }
	int KeywordId(void) { return keyword_id; }
	String &Token(void) { return next_token; }

public:
	String text;
	String next_token;

private:
	bool HasOnlySpaces(void);

private:
	int keyword_id;
	RETURN_TYPE line_type;
	RETURN_TYPE token_type;
	long next_token_position;
	List<Keyword> *keyword_list;
};

#endif