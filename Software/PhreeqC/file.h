#ifndef fileH
#define fileH

//#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
//#include <fstream>

#include "line.h"

using namespace std;

class File
{
public:
	File(Line *line);
	File(string &name, Line *line);
	~File();

public:
	void Reset(void);
	void OpenFile(string &file_name, bool throw_exception = false);
	void CloseFile(void);
	bool GetLine(void);
	bool IsOpen(void) { return is_open; }
	bool Eof(void) { return eof; }

private:
	bool eof;
	bool is_open;
	FILE *file;
	Line *line;
};

#endif