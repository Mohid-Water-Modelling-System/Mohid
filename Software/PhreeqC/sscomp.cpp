#include "datastructures.h"

//===========================================================================================================
// Constructs
//===========================================================================================================
SSComp::SSComp()
{
	Reset();
}
//-----------------------------------------------------------------------------------------------------------
SSComp::SSComp(SSComp *copy)
{
	copy->CopyTo(this);
}
//-----------------------------------------------------------------------------------------------------------
SSComp::~SSComp()
{
}
//-----------------------------------------------------------------------------------------------------------
void SSComp::Reset()
{
	initial_moles = 0;
	delta = 0;

	init_moles = 0;
	moles = (LDBLE)NAN;
	fraction_x = 0;
	log10_lambda = 0;
	log10_fraction_x = 0;
	dn = 0; 
	dnc = 0;
	dnb = 0;

	phase = NULL;

	name = "";
}
//-----------------------------------------------------------------------------------------------------------
void SSComp::CopyTo(SSComp *copy)
{
	copy->initial_moles = initial_moles;
	copy->init_moles = init_moles;
	copy->delta = delta;
	copy->moles = moles;
	copy->fraction_x = fraction_x;
	copy->log10_lambda = log10_lambda;
	copy->log10_fraction_x = log10_fraction_x;
	copy->dn = dn; 
	copy->dnc = dnc;
	copy->dnb = dnb;

	copy->phase = phase;

	copy->name = name;
}
//-----------------------------------------------------------------------------------------------------------
/*
SSComp *SSComp::Copy()
{
	return new SSComp(this);
}
*/
//-----------------------------------------------------------------------------------------------------------
