#ifndef messagesH
#define messagesH

#include "string.h"
#include <stdio.h>
#include <stdlib.h>

class Messages
{
public:
	Messages();
	Messages(const char *message_path);
	~Messages();

public:
	void Reset();
	void SetPath(const char *message_path);
	int Open();
	void Close();
	void Write(const char *message);
	void Write(String &message);

protected:
	String message_path;
	FILE *message_file;
	
};

#endif