#include "datastructures.h"

ChargeGroup::ChargeGroup()
{
	z = 0.0;
	eq = 0.0;
}

ChargeGroup::ChargeGroup(ChargeGroup *copy)
{
	copy->CopyTo(this);
}

ChargeGroup::~ChargeGroup()
{
}

void ChargeGroup::Reset()
{
	z = 0.0;
	eq = 0.0;
}

void ChargeGroup::CopyTo(ChargeGroup *copy)
{
	copy->z = z;
	copy->eq = eq;
}

/*
ChargeGroup *ChargeGroup::Copy()
{
	return new ChargeGroup(this);
}
*/

