#ifndef basicfunctionsH
#define basicfunctionsH

#include "constants.h"
#include "string.h"
#include "list.h"

class BasicFunctions
{
public:
	bool GetNumber(char **t_ptr, LDBLE &num);
	void CorrectPlusOne(String &to_correct);
	bool Equal(const LDBLE a, const LDBLE b, const LDBLE tol);
	bool GetElement(char **str, String &element);
	bool GetSecondary(char **str, String &element);
	bool GetSpeciesNameAndZ(char **str, String &name, LDBLE &z);
	bool GetZ(String &token, LDBLE &z);


protected:
	bool GetCharge(String &charge, LDBLE &z);


};

#endif