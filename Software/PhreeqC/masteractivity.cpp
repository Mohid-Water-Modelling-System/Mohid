#include "datastructures.h"

MasterActivity::MasterActivity()
{
	Reset();
}

MasterActivity::MasterActivity(MasterActivity *copy)
{
	copy->CopyTo(this);
}

MasterActivity::~MasterActivity()
{
}

void MasterActivity::Reset()
{
	description = "";
	la = 0.0;
}

/*
MasterActivity *MasterActivity::Copy()
{
	return new MasterActivity(this);
}
*/

void MasterActivity::CopyTo(MasterActivity * copy)
{
	copy->la = la;
	copy->description = description;
}

