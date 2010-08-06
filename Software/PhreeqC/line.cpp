#include "line.h"

Line::Line()
{
	next_token_position = 0;

	keyword_list = NULL;
	keyword_list = new List<Keyword>;

	keyword_list->AddNew("end");
	keyword_list->AddNew("solution_species");
	keyword_list->AddNew("solution_master_species");
	keyword_list->AddNew("phases");	                     
	keyword_list->AddNew("exchange_master_species");
	keyword_list->AddNew("exchange_species");
	keyword_list->AddNew("surface_master_species");
	keyword_list->AddNew("surface_species");

	for (int i = 0; i < keyword_list->Count(); i++)
		(*keyword_list)[i]->id = i + 1;
}

Line::~Line()
{
	if (keyword_list != NULL)
		delete keyword_list;
}

void Line::Reset(void)
{
}

bool Line::NextToken(String &token_delimiter)
{
	token_type = _empty_;
	next_token = "";

	if (text.IsEmpty())
		return false;

	if (next_token_position >= text.Length())
		return false;

	// Check what we have
	if (isupper ((int) text[next_token_position]) || text[next_token_position] == '[')
		token_type = _upper_;
	else if (islower ((int) text[next_token_position]))
		token_type = _lower_;
	else if (isdigit ((int) text[next_token_position]) || text[next_token_position] == '.' || text[next_token_position] == '-')
		token_type = _digit_;
	else if (text[next_token_position] == '\0')
		token_type = _empty_;
	else
		token_type = _unknown_;

	long l;
	for (l = next_token_position; l < text.Length() && !text.IsAmong(l, token_delimiter); l++);

	next_token = text(next_token_position, l - next_token_position);

	next_token_position = l;
	for (; next_token_position < text.Length() && text.IsAmong(next_token_position, token_delimiter); next_token_position++);

	return true;
}

void Line::SetLine(String &text_to_set)
{
	char tab[] = {char(9), '\0'};

	text = text_to_set;
	text.Replace(tab, " ");
	text.Trim(_all_);

	Keyword *k_ptr;

	keyword_id = -1;
	if (text.IsEmpty() || HasOnlySpaces())
		line_type = _empty_;
	else if ((k_ptr = keyword_list->Search(&text, keyword_id, true)) != NULL)
		line_type = _keyword_;
	else if (text.Length() >= 2 && text[0] == '-' && isalpha((int) text[1]))
		line_type = _option_;
	else
		line_type = _none_;

	next_token_position = 0;
}

bool Line::HasOnlySpaces(void)
{
	if (text.IsEmpty())
		return false;

	for (long i = 0; i < text.Length(); i++)
		if (text[i] != ' ')
			return false;

	return true;
}