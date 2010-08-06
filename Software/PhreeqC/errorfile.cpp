#include "errors.h"

ErrorFile::ErrorFile(String error_file)
{
	error_f = NULL;
	error_f = fopen(error_file.CharPtr(), "wt");

	if (error_f == NULL)
		throw exception("\nCan't open error file.\n");

}

ErrorFile::~ErrorFile()
{
	if (error_f != NULL)
		fclose(error_f);
}

void ErrorFile::SaveMessage(String message)
{
	fprintf(error_f, "%s", message.CharPtr());
}

void ErrorFile::SaveMessage(char *message)
{
	fprintf(error_f, "%s", message);
}
