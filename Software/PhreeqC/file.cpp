#include "file.h"

File::File(Line *line)
{
	if (line == NULL)
		throw 1;

	this->line = line;

	this->eof     = false;
	this->is_open = false;
}

File::File(string &name, Line *line)
{
	if (line == NULL)
		throw 1;

	this->line = line;

	this->eof     = false;
	this->is_open = false;

	OpenFile(name);
}

File::~File()
{
	CloseFile();
}

void File::OpenFile(string &file_name, bool throw_exception)
{
	if (is_open)
		CloseFile();

	file = NULL;
	file = fopen(file_name.c_str(), "rt");
	//file.open(file_name.c_str());

	if (file == NULL || feof(file))
	{
		this->eof     = true;
		this->is_open = false;

		if (throw_exception)
			throw 1; //change this
		else
			return;
	}

	this->eof     = false;
	this->is_open = true;
}

void File::CloseFile(void)
{
	if (is_open)
		fclose(file);

	this->eof     = false;
	this->is_open = false;
}

bool File::GetLine(void)
{
	string line_t, line_d;
	bool ignore_eol = false;
  char c;

	if (feof(file))
		return false;

	line_d = "";

	do
	{
		c = fgetc(file);
		if (c == EOF)
			break;

		if (c == '#')
		{
			/* ignore all chars after # until newline */
			while (!feof(file) && c != '\n')
				c = fgetc(file);	

			if (!ignore_eol && line->text.Length() > 0)
				break;
			else if(ignore_eol)
				ignore_eol = false;
		}
		else if (c == ';')
		{
			ignore_eol = false;
			break;
		}
		else if (c == '\n')
		{
			if (!ignore_eol)
				break;

			ignore_eol = false;
		}
		else if (c == '\\')
			ignore_eol = true;
		else
			line_d += c;
	}
	while (!feof(file));

	if (feof(file))
		this->eof = true;
	else
		this->eof = false;

	line->SetLine(String(line_d.c_str()));
	return true;
}