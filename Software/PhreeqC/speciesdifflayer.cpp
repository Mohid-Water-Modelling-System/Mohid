#include "datastructures.h"	

SpeciesDiffLayer::SpeciesDiffLayer()
{
	charge = NULL;
	charge = new SurfaceCharge;

	Reset();
}

SpeciesDiffLayer::SpeciesDiffLayer(SpeciesDiffLayer *copy)
{
	charge = NULL;
	charge = new SurfaceCharge;

	copy->CopyTo(this);
}

SpeciesDiffLayer::~SpeciesDiffLayer()
{
	delete charge;
}

void SpeciesDiffLayer::Reset()
{
	g_moles = 0.0;
	dg_g_moles = 0.0;
	dx_moles = 0.0;
	drelated_moles = 0.0;
	dh2o_moles = 0.0;

	charge->Reset();

	count_g = 0;
}

/*
SpeciesDiffLayer *SpeciesDiffLayer::Copy()
{
	return new SpeciesDiffLayer(this);
}
*/

void SpeciesDiffLayer::CopyTo(SpeciesDiffLayer *copy)
{
	copy->g_moles = g_moles;
	copy->dg_g_moles = dg_g_moles;
	copy->dx_moles = dx_moles;
	copy->drelated_moles = drelated_moles;
	copy->dh2o_moles = dh2o_moles;

	charge->CopyTo(copy->charge);

	copy->count_g = count_g;
}

void SpeciesDiffLayer::PrintToFile(FILE *file, int spaces)
{
	String spaces_str(spaces, " ");

	fprintf(file, "%sg_moles: %f\n", spaces_str.CharPtr(), g_moles);
	fprintf(file, "%sdg_g_moles: %f\n", spaces_str.CharPtr(), dg_g_moles);
	fprintf(file, "%sdx_moles: %f\n", spaces_str.CharPtr(), dx_moles);
	fprintf(file, "%sdrelated_moles: %f\n", spaces_str.CharPtr(), drelated_moles);
	fprintf(file, "%sdh2o_moles: %f\n", spaces_str.CharPtr(), dh2o_moles);
	
	fprintf(file, "%scount_g: %f\n", spaces_str.CharPtr(), count_g);

	if (charge != NULL)
	{
		fprintf(file, "%s  CHARGE: (%p)\n", spaces_str.CharPtr(), charge);
		charge->PrintToFile(file, spaces + 5);
	}
	else
		fprintf(file, "%s  CHARGE: NULL\n", spaces_str.CharPtr(), charge);
}