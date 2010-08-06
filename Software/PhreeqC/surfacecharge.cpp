#include "datastructures.h"

SurfaceCharge::SurfaceCharge()
{
	diffuse_layer_totals = NULL;
	g = NULL;

	diffuse_layer_totals = new List<ElementOfSpecies>;
	g = new List<SurfaceDiffLayer>;

	Reset();
}

SurfaceCharge::SurfaceCharge(SurfaceCharge *copy)
{
	diffuse_layer_totals = NULL;
	g = NULL;

	diffuse_layer_totals = new List<ElementOfSpecies>;
	g = new List<SurfaceDiffLayer>;

	copy->CopyTo(this);
}

SurfaceCharge::~SurfaceCharge()
{
	delete diffuse_layer_totals;
	delete g;
}

void SurfaceCharge::Reset()
{
  name = "";

  specific_area = 0.0;
  grams = 0.0;
  charge_balance = 0.0;
  mass_water = 0.0;
  
	diffuse_layer_totals->Clear();
  g->Clear();
  
	la_psi = 0.0;
	la_psi1 = 0.0;
	la_psi2 = 0.0;
  psi = 0.0;
	psi1 = 0.0;
	psi2 = 0.0;

	capacitance[0] = 1.0;
	capacitance[1] = 5.0;

  sigma0 = 0.0;
	sigma1 = 0.0;
	sigma2 = 0.0;
	sigmaddl = 0.0;
}

/*
SurfaceCharge *SurfaceCharge::Copy()
{
	return new SurfaceCharge(this);
}
*/

void SurfaceCharge::CopyTo(SurfaceCharge *copy)
{
  copy->name = name;

  copy->specific_area = specific_area;
  copy->grams = grams;
  copy->charge_balance = charge_balance;
  copy->mass_water = mass_water;
  
	diffuse_layer_totals->CopyTo(copy->diffuse_layer_totals);
  g->CopyTo(copy->g);
  
	copy->la_psi = la_psi;
	copy->la_psi1 = la_psi1;
	copy->la_psi2 = la_psi2;
  copy->psi = psi;
	copy->psi1 = psi1;
	copy->psi2 = psi2;

	for (int i = 0; i < 2; i++)
		copy->capacitance[i] = capacitance[i];

  copy->sigma0 = sigma0;
	copy->sigma1 = sigma1;
	copy->sigma2 = sigma2;
	copy->sigmaddl = sigmaddl;
}

void SurfaceCharge::PrintToFile(FILE *file, int spaces)
{
}