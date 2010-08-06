#include "exceptionhandler.h"

ExceptionHandler::ExceptionHandler(const char *message)
{
	this->message = message;
}

const char * ExceptionHandler::Message()
{
	return message;
}