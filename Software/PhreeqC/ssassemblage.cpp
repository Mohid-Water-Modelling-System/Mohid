#include "datastructures.h"

//===========================================================================================================
// Constructs
//===========================================================================================================
SSAssemblage::SSAssemblage()
{
	ss_list = NULL;
	ss_list = new List<SS>;

	Reset();
}
SSAssemblage::SSAssemblage(SSAssemblage *copy)
{
	ss_list = NULL;
	ss_list = new List<SS>;

	copy->CopyTo(this);
}
//-----------------------------------------------------------------------------------------------------------
SSAssemblage::~SSAssemblage()
{
	delete ss_list;
}
//-----------------------------------------------------------------------------------------------------------
void SSAssemblage::Reset()
{
	number = -1;

	ss_list->Clear();
}

void SSAssemblage::CopyTo(SSAssemblage *copy)
{
	copy->number = number;
	ss_list->CopyTo(copy->ss_list);
}

/*
SSAssemblage *SSAssemblage::Copy()
{
	return new SSAssemblage(this);
}
*/
