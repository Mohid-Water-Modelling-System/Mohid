#include "datastructures.h"

SurfaceDiffLayer::SurfaceDiffLayer()
{
	Reset();
}

SurfaceDiffLayer::SurfaceDiffLayer(SurfaceDiffLayer *copy)
{
	copy->CopyTo(this);
}

SurfaceDiffLayer::~SurfaceDiffLayer()
{
}

void SurfaceDiffLayer::Reset()
{
  charge = 0.0;
  g = 0.0;
  dg = 0.0;
  psi_to_z = 0.0;
}

/*
SurfaceDiffLayer *SurfaceDiffLayer::Copy()
{
	return new SurfaceDiffLayer(this);
}
*/

void SurfaceDiffLayer::CopyTo(SurfaceDiffLayer *copy)
{
  copy->charge = charge;
  copy->g = g;
  copy->dg = dg;
  copy->psi_to_z = psi_to_z;
}

