#include "messages.h"

Messages::Messages()
{
	message_file = NULL;
	message_path = "";
}

Messages::Messages(const char *message_path)
{
	message_file = NULL;
	SetPath(message_path);
}

Messages::~Messages()
{
	if (message_file != NULL)
		fclose(message_file);
}

void Messages::SetPath(const char *message_path)
{
	this->message_path.SetString(message_path);
}

/*
	Possible return values by Open:
		 0: File message_path not specified
		-1: Can't open message_file
		 1: File opened to write
*/
int Messages::Open()
{
	if (message_file != NULL)
		Close();

	if (message_path.IsEmpty())
		return 0;

	message_file = fopen(message_path.CharPtr(), "wt");
	if (message_file == NULL)
		return -1;

	return 1;
}

void Messages::Close()
{
	if (message_file != NULL)
	{
		fclose(message_file);
		message_file = NULL;
	}
}

void Messages::Write(const char *message)
{
	if (message != NULL && message_file != NULL)
		fprintf(message_file, "%s", message);
}

void Messages::Write(String &message)
{
	Write(message.CharPtr());
}