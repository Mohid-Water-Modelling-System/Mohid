#ifndef errorsH
#define errorsH

#include <string.h>
#include <error.h>

#include "string.h"

class ErrorFile
{
public:
	ErrorFile(String error_file);
	~ErrorFile();

public:
	void SaveMessage(String message);
	void SaveMessage(char *message);

private:
	FILE *error_f;

};

class ErrorHandler 
{
public:
	ErrorHandler(const char *message, int error_code = 0)
	{
		SetMessage(message);
		this->error_code = error_code;
	}

	ErrorHandler()
	{
		message = NULL;
		error_code = 0;
	}

	~ErrorHandler()
	{
		if (message != NULL) delete [] message;
	}

	virtual char * Message() = 0;
	virtual int Code() = 0;

protected:
	void SetMessage(const char *message)
	{
		if (message != NULL)
		{
			size_t length = strlen(message);
			this->message = new char [length + 1];
			strcpy(this->message, message);
		}
		else
			this->message = NULL;
}

protected:
	char *message;
	int error_code;	
	int function;

};

typedef enum 
{
	NEGATIVE_CONCENTRATION_MEE = 1
} MODEL_ENGINE_ERRORS;

typedef enum
{
	UNKNOWN_MEF,
	RUN_REACTIONS_MEF, 
	SET_AND_RUN_WRAPPER_MEF
} MODEL_ENGINE_FUNCTION;

class EModelEngine:
	public ErrorHandler
{
public:
	EModelEngine(int code, MODEL_ENGINE_FUNCTION function = UNKNOWN_MEF): ErrorHandler()
	{
		switch(error_code = code)
		{
		case 1:
			SetMessage("negative concentration in system.");
			break;
		case 2:
			SetMessage("numerical method failed on all combinations of convergence parameters");
			break;
		case 3:
			break;
		default:
			SetMessage("no message");
		};

		this->function = function;
	}

	virtual char * Message()
	{
		switch(function)
		{
		case RUN_REACTIONS_MEF:
			return "RunReactions function";
		case SET_AND_RUN_WRAPPER_MEF:
			return "SetAndRunWrapper";
		default: 
			return "unknown function";
		}
	}

	virtual int Code()
	{
		return function;
	}
};

#endif