#include "stdafx.h"
#include "Generic.h"

void Generic::CatchException(exception e)
{
	cout << e.what() << endl;
	CuWrapper::CuCleanUp();
	exit(1);
}

/// Print is defined in Generic.h