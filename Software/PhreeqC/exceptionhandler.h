#pragma once

#ifndef exceptionhandlerH
#define exceptionhandlerH

class ExceptionHandler
{
public:
	ExceptionHandler(const char *message);
	
public:
	const char *Message();

private:
	const char *message;
};

#endif
