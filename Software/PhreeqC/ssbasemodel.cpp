#include "ssbasemodel.h"
#include <math.h>

//===========================================================================================================
// Constructs
//===========================================================================================================
SSBaseModel::SSBaseModel(GlobalData *gd, ModelData *md):
ParseEngine(gd)
{
	this->gd = gd;
	this->md = md;
}
//-----------------------------------------------------------------------------------------------------------
SSBaseModel::~SSBaseModel()
{
}
//-----------------------------------------------------------------------------------------------------------

//===========================================================================================================
// Interface
//===========================================================================================================
bool SSBaseModel::SSPrep(LDBLE t, SS *ss_p)
{
  int i, j, k, divisions;
	bool converged;
  LDBLE r, rt, ag0, ag1, crit_pt;
  LDBLE xc, tc;
  LDBLE x, x0, x1, xsm1, xsm2, xb1, xb2;
  LDBLE xc1, xc2;
  LDBLE facb1, faca1, spim1, xblm1, acrae, acrael, xliapt, xliapm;
  LDBLE xaly, xaly1, xaly2;
  LDBLE faca, facb, spialy, facal, facbl;
  LDBLE tol;

  tol = (LDBLE)1e-6;
  r = R_KJ_DEG_MOL;
  rt = r * t;
  
	a0 = ss_p->ag0 / rt;
  a1 = ss_p->ag1 / rt;
  
	ss_p->a0 = a0;
  ss_p->a1 = a1;

  ag0 = a0 * rt;
  ag1 = a1 * rt;

  kc = exp (KCalc ((*ss_p->comps_list)[0]->phase->rxn->logk, t) * md->LOG_10);
  kb = exp (KCalc ((*ss_p->comps_list)[1]->phase->rxn->logk, t) * md->LOG_10);
  
	crit_pt = fabs (a0) + fabs (a1);

  ss_p->miscibility = false;
  ss_p->spinodal = false;

  xsm1 = 0.5;
  xsm2 = 0.5;

  xb1 = 0.5;
  xb2 = 0.5;

  xc1 = 0;
  xc2 = 0;

  if (crit_pt >= tol)
  {
    if (fabs (a1) < tol)
    {
      xc = 0.5;
      tc = ag0 / (2 * r);
    }
    else
    {
      xc = (LDBLE)0.5 + (pow ((LDBLE)(ag0 * ag0 + 27 * ag1 * ag1), (LDBLE)0.5) - ag0) / (18 * ag1);
      tc = (12 * ag1 * xc - 6 * ag1 + 2 * ag0) * (xc - xc * xc) / r;
    }

    if (tc >= t)
    {
      x0 = 0;
      x1 = 1;

      if (Scan(&x0, &x1))
      {
				xsm1 = Halve(x0, x1, tol);
				ss_p->spinodal = true;

				x0 = x1;
				x1 = 1;

				if (Scan(&x0, &x1))
					xsm2 = Halve(x0, x1, tol);
				else
					return false;
      }
    }
  }

  if (ss_p->spinodal)
  {
    converged = false;

    for (i = 1; i < 3; i++)
    {
			divisions = (int) pow (10., i);

			for (j = 0; j < divisions; j++)
			{
				for (k = divisions; k > 0; k--)
				{
					xc1 = (LDBLE) j / divisions + (LDBLE)0.001;
					xc2 = (LDBLE) k / divisions;
					converged = SolveMisc(&xc1, &xc2, tol);

					if (converged)
						break;
				}

				if (converged)
					break;
			}

			if (converged)
				break;
    }

    if (!converged)
			return false;

    ss_p->miscibility = true;

    if (xc1 < xc2)
    {
      xb1 = 1 - xc2;
      xb2 = 1 - xc1;
      xc1 = 1 - xb1;
      xc2 = 1 - xb2;
    }
    else
    {
      xb1 = 1 - xc1;
      xb2 = 1 - xc2;
    }

    facb1 = kb * xb1 * exp (xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
    faca1 = kc * xc1 * exp (xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));

    spim1 = log10 (faca1 + facb1);
    xblm1 = (LDBLE)1. / ((LDBLE)1. + faca1 / facb1);

    acrae = facb1 / faca1;
    acrael = log10 (acrae);

    xliapt = log10 (facb1);
    xliapm = log10 (faca1);

    ss_p->tk = t;
    ss_p->xb1 = xb1;
    ss_p->xb2 = xb2;
  }

  xaly = -1.0;
  x = a0 * a0 + 3 * a1 * a1 + 6 * a1 * log (kb / kc);

  if (x > 0)
  {
    if (fabs (x - a0 * a0) >= tol)
    {
      xaly1 = (-(a0 - 3 * a1) + pow ((LDBLE)x, (LDBLE)0.5)) / (6 * a1);
      xaly2 = (-(a0 - 3 * a1) - pow ((LDBLE)x, (LDBLE)0.5)) / (6 * a1);

      if (xaly1 >= 0 && xaly1 <= 1)
      {
				xaly = xaly1;
      }

      if (xaly2 >= 0 && xaly2 <= 1)
      {
				xaly = xaly2;
      }
    }
    else
    {
      xaly = (LDBLE)0.5 + log (kb / kc) / ((LDBLE)2 * a0);
    }

    if (xaly > 0 && xaly < 1)
    {
      faca = kc * (1 - xaly) * exp (xaly * xaly * (a0 - a1 * (3 - 4 * xaly)));
      facb = kb * xaly * exp (((LDBLE)1 - xaly) * ((LDBLE)1 - xaly) * (a0 + a1 * ((LDBLE)4 * xaly - (LDBLE)1.0)));

      spialy = log10 (faca + facb);

      facal = log10 (faca);
      facbl = log10 (facb);
    }
  }

  return true;
}
//-----------------------------------------------------------------------------------------------------------

//===========================================================================================================
// Interface
//===========================================================================================================
bool SSBaseModel::Scan(LDBLE * xx0, LDBLE * xx1)
{
  int i, j, divisions;
  LDBLE fx0, fx1;
  LDBLE x0, x1, diff;

  x0 = *xx0;
  x1 = *xx1;
  diff = x1 - x0;
	
  for (j = 0; j < 3; j++)
  {
    fx0 = FSpinodal(x0);
    divisions = (int) pow ((LDBLE) 10, (LDBLE) j);
		
    for (i = 1; i < divisions; i++)
    {
      x1 = *xx0 + diff * (LDBLE) i / (LDBLE)divisions;
      fx1 = FSpinodal(x1);
			
      if (fx0 * fx1 <= 0)
      {
				*xx0 = x0;
				*xx1 = x1;
				
				return true;
      }
			
      x0 = x1;
      fx0 = fx1;
    }
  }
	
  return false;
}
//-----------------------------------------------------------------------------------------------------------

LDBLE SSBaseModel::FSpinodal(LDBLE x)
{
  return ((LDBLE)-12 * a1 * x * x * x + ((LDBLE)18 * a1 - (LDBLE)2 * a0) * x * x + ((LDBLE)2 * a0 - (LDBLE)6 * a1) * x - (LDBLE)1.0);
}

//-----------------------------------------------------------------------------------------------------------
LDBLE SSBaseModel::Halve(LDBLE x0, LDBLE x1, LDBLE tol)
{
  int i;
  LDBLE x, y, y0, dx;

  y0 = FSpinodal(x0);
  dx = (x1 - x0);
	
  for (i = 0; i < 100; i++)
  {
    dx *= 0.5;
    x = x0 + dx;
    y = FSpinodal(x);
		
    if (dx < tol || y == 0)
    {
      break;
    }
		
    if (y0 * y >= 0)
    {
      x0 = x;
      y0 = y;
    }
  }
	
  return (x0 + dx);
}

//-----------------------------------------------------------------------------------------------------------
bool SSBaseModel::SolveMisc(LDBLE * xxc1, LDBLE * xxc2, LDBLE tol)
{
  int i, max_iter;
	bool repeat, converged;
  LDBLE x1, x2, xb1, xb2;
  LDBLE xc1, xc1_2, xc1_3, xc2, xc2_2, xc2_3;
  LDBLE lc1, lc2, lb1, lb2;
  LDBLE a[6], d[2];
  LDBLE t;

  xc1 = *xxc1;
  xc2 = *xxc2;
  x1 = 0;
  x2 = 0;

  converged = true;
  max_iter = 25;
	
  for (i = 0; i < max_iter; i++)
  {
    xb1 = 1 - xc1;
    xb2 = 1 - xc2;
    xc1_2 = xc1 * xc1;
    xc1_3 = xc1_2 * xc1;
    xc2_2 = xc2 * xc2;
    xc2_3 = xc2_2 * xc2;
		
    lc1 = exp (xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
    lb1 = exp (xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
    lc2 = exp (xb2 * xb2 * (a0 - a1 * (3 - 4 * xb2)));
    lb2 = exp (xc2 * xc2 * (a0 + a1 * (4 * xb2 - 1)));
		
    a[2] = -(xb1 * lb1 - xb2 * lb2);	
    a[5] = -(xc1 * lc1 - xc2 * lc2);
		
    if (fabs (a[2]) < tol && fabs (a[5]) < tol)
      break;
			
    t = exp (a0 * xc1_2 - 4 * a1 * xc1_3 + 3 * a1 * xc1_2);
    a[0] = (2 * a0 * xc1 + 6 * a1 * xc1 - 2 * a0 * xc1_2 + 12 * a1 * xc1_3 - 18 * a1 * xc1_2 - 1) * t;
		
    t = exp (a0 * xc2_2 - 4 * a1 * xc2_3 + 3 * a1 * xc2_2);
    a[1] = (2 * a0 * xc2_2 - 12 * a1 * xc2_3 - 2 * a0 * xc2 + 18 * a1 * xc2_2 - 6 * a1 * xc2 + 1) * t;
		
    t = exp (a0 * xc1_2 - 2 * a0 * xc1 + a0 - 4 * a1 * xc1_3 + 9 * a1 * xc1_2 - 6 * a1 * xc1 + a1);
    a[3] = (2 * a0 * xc1_2 - 2 * a0 * xc1 - 12 * a1 * xc1_3 + 18 * a1 * xc1_2 - 6 * a1 * xc1 + 1) * t;
		
    t = exp (a0 * xc2_2 - 2 * a0 * xc2 + a0 - 4 * a1 * xc2_3 + 9 * a1 * xc2_2 - 6 * a1 * xc2 + a1);
    a[4] = (-2 * a0 * xc2_2 + 2 * a0 * xc2 + 12 * a1 * xc2_3 - 18 * a1 * xc2_2 + 6 * a1 * xc2 - 1) * t;
		
    if (!SLNQ(2, a, d, 3))
			return false;
		
    repeat = true;
    while (repeat)
    {
      x1 = xc1 + d[0];
      x2 = xc2 + d[1];
			
      if (x1 > 1 || x1 < 0 || x2 > 1 || x2 < 0)
      {
				d[0] *= 0.5;
				d[1] *= 0.5;
      }
      else
      {
				repeat = false;
      }
    };
		
    xc1 = x1;
    xc2 = x2;
		
    if (fabs (xc1 - xc2) < .01)
    {
      converged = false;
      break;
    }
  }
	
  if (i == max_iter)
    converged = false;
		
  *xxc1 = xc1;
  *xxc2 = xc2;
	
  return converged;
}

//-----------------------------------------------------------------------------------------------------------
bool SSBaseModel::SLNQ(int n, LDBLE * a, LDBLE * delta, int ncols)
{
  int i, j, k, m;
  LDBLE b;

	if (n == 0)
    return true;

  if (n == 1)
  {
    if (fabs (a[0]) < ZERO_TOL)
      return false;

    delta[0] = a[1] / a[0];
    return true;
  }

  for (i = 0; i < n - 1; i++)
  {
    b = fabs (a[i * ncols + i]);
    m = i;

    for (j = i + 1; j < n; j++)
    {
      if (fabs (a[j * ncols + i]) > b)
      {
				b = fabs (a[j * ncols + i]);
				m = j;
      }
    }

    if (b < ZERO_TOL)
      return false;

    if (m != i)
    {
      for (j = i; j <= n; j++)
      {
				b = a[i * ncols + j];
				a[i * ncols + j] = a[m * ncols + j];
				a[m * ncols + j] = b;
      }
    }

    for (j = n; j >= i; j--)
      a[i * ncols + j] /= a[i * ncols + i];

    for (j = i + 1; j < n; j++)
    {
      if (a[j * ncols + i] == 0.0)
				continue;

      b = -a[j * ncols + i];

      for (k = i + 1; k <= n; k++)
				a[j * ncols + k] += b * a[i * ncols + k];
    }
  }

  if (fabs (a[(n - 1) * ncols + n - 1]) > ZERO_TOL)
    delta[n - 1] = a[(n - 1) * ncols + n] / a[(n - 1) * ncols + n - 1];
  else
    return false;

  for (i = n - 2; i >= 0; i--)
  {
    delta[i] = a[i * ncols + n];

    for (j = i + 1; j < n; j++)
      delta[i] -= a[i * ncols + j] * delta[j];
  }

  return true;
}

LDBLE SSBaseModel::KCalc(LDBLE * logk, LDBLE tempk)
{
  return (
		      logk[0] - logk[1] * ((LDBLE)298.15 - tempk) / ((LDBLE)298.15 * tempk * md->LOG_10 * R_KJ_DEG_MOL) +
	        logk[2] + logk[3] * tempk + logk[4] / tempk + logk[5] * log10 (tempk) + logk[6] / (tempk * tempk)
				 );
}