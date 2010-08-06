#ifndef dataH
#define dataH

class Data
{
public:
	void FreeAndNull(void *pointer) 
	{
		delete pointer;
		pointer = NULL;
	}

};

#endif