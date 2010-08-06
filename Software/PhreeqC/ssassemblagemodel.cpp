#include "ssassemblagemodel.h"
#include <math.h>

//===========================================================================================================
// Constructs
//===========================================================================================================
SSAssemblageModel::SSAssemblageModel(GlobalData *gd, ModelData *md)
: SSBaseModel(gd, md)
{
	this->gd = gd;
	this->md = md;
}
//-----------------------------------------------------------------------------------------------------------
SSAssemblageModel::~SSAssemblageModel()
{
}
//-----------------------------------------------------------------------------------------------------------

//===========================================================================================================
// Interface
//===========================================================================================================
//-----------------------------------------------------------------------------------------------------------

//===========================================================================================================
// Internal functions
//===========================================================================================================
bool SSAssemblageModel::SSCalcA0A1(SS *ss)
{
  int i;
	bool done;
  LDBLE r, rt, *p;
  LDBLE q1, q2, xbq1, xbq2, xb1, xb2, xc1, xc2;
  LDBLE r1, r2, pa1, pb1, pa2, pb2, xsm1, xsm2;
  LDBLE pn9, pn10, c5, c6, pl9, pl10, pj9, pj10;
  LDBLE xc, tc;
  LDBLE spialy, azero, phi1, phi2, test;
  LDBLE dq1, dq2, denom, ratio, dr1, dr2, x21, x22, x61, x62;
  LDBLE a0, a1, ag0, ag1;
  LDBLE wg2, wg1, alpha2, alpha3;
  LDBLE kc, kb;
  LDBLE xaly, xcaly, alpha0, alpha1, fx, fx1;
  LDBLE tol;

  tol = (LDBLE)1e-6;
  rt = ss->tk * R_KJ_DEG_MOL;
	
	SSComp *ssc_p0 = (*ss->comps_list)[0],
				 *ssc_p1 = (*ss->comps_list)[1];

  if (ssc_p0->phase == NULL || ssc_p1->phase == NULL)
		return false;
	
  kc = exp (KCalc (ssc_p0->phase->rxn->logk, ss->tk) * md->LOG_10);
  kb = exp (KCalc (ssc_p1->phase->rxn->logk, ss->tk) * md->LOG_10);

  p = ss->p;

  a0 = 0;
  a1 = 0;
  ag0 = 0;
  ag1 = 0;
  dq2 = 0;
	
  switch (ss->input_case)
  {
  case 0:
    a0 = p[0];
    a1 = p[1];
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  case 1:
    q1 = p[0];
    q2 = p[1];
    xbq1 = p[2];
    xbq2 = p[3];
    done = false;
		
    if (fabs (1 - xbq1) > 0 && q1 > 0)
    {
      dq1 = log (q1) / ((1 - xbq1) * (1 - xbq1));
			
      if (xbq2 <= 0 || xbq2 > 1)
      {
				a0 = dq1;
				a1 = 0;
				done = true;
      }
    }
		
    if (!done)
    {
      if (fabs (xbq2) < 0 || q2 <= 0)
				return false;
    }
		
    if (!done)
    {
      dq2 = log (q2) / (xbq2 * xbq2);
			
      if (xbq1 < 0. || xbq2 > 1.)
      {
				a0 = dq2;
				a1 = 0;
				done = true;
      }
    }
		
    if (!done)
    {
      denom = 4 * (xbq1 - xbq2) + 2;
			
      if (fabs (denom) >= tol)
      {
				if (fabs (1 - xbq1) > 0 && q1 > 0)
				{
					dq1 = log (q1) / ((1 - xbq1) * (1 - xbq1));
					a0 = (dq1 * (3 - 4 * xbq2) + dq2 * (4 * xbq1 - 1)) / denom;
					a1 = (dq1 - dq2) / denom;
					done = true;
				}
      }
    }
		
    if (!done)
			return false;
		
    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  case 2:
    q1 = p[0];
    q2 = p[1];
    xbq1 = p[2];
    xbq2 = p[3];
    ratio = kc / kb;
    dr1 = log (q1 / ratio);
    x21 = 2 * xbq1 - 1;
		
    if (fabs (xbq1 - xbq2) < tol || xbq2 < 0)
    {
      a0 = dr1 / x21;
      a1 = 0;
    }
    else
    {
      dr2 = log (q2 / ratio);
      x22 = 2 * xbq2 - 1;
			
      if (xbq1 < 0.)
      {
				a0 = dr2 / x22;
				a1 = 0;
      }
      else
      {
				x61 = 6 * xbq1 * xbq1 - 6 * xbq1 + 1;
				x62 = 6 * xbq2 * xbq2 - 6 * xbq2 + 1;
				
				if (fabs (x22 * x61 - x21 * x62) < tol)
					return false;
				
				a0 = (x61 * dr2 - x62 * dr1) / (x22 * x61 - x21 * x62);
				a1 = (x21 * dr2 - x22 * dr1) / (x21 * x62 - x22 * x61);
      }
    }

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  case 3:
    q1 = p[0];
    q2 = p[1];
    xb1 = q1;
    xb2 = q2;
    xc1 = 1 - xb1;
    xc2 = 1 - xb2;
    r1 = log (xb1 / xb2);
    r2 = log (xc1 / xc2);
    pa1 = xc2 * xc2 - xc1 * xc1;
    pb1 = 3 * (xc2 * xc2 - xc1 * xc1) - 4 * (xc2 * xc2 * xc2 - xc1 * xc1 * xc1);
    pa2 = xb2 * xb2 - xb1 * xb1;
    pb2 = -(3 * (xb2 * xb2 - xb1 * xb1) - 4 * (xb2 * xb2 * xb2 - xb1 * xb1 * xb1));
    a0 = (r1 - pb1 / pb2 * r2) / (pa1 - pa2 * pb1 / pb2);
    a1 = (r1 - pa1 / pa2 * r2) / (pb1 - pb2 * pa1 / pa2);

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  case 4:
    q1 = p[0];
    q2 = p[1];
    xsm1 = q1;
    xsm2 = q2;
    pn9 = 1 / xsm1;
    pn10 = 1 / xsm2;
    c5 = 1 - xsm1;
    c6 = 1 - xsm2;
    pl9 = 6 * c5 - 12 * c5 * c5;
    pl10 = 6 * c6 - 12 * c6 * c6;
    pj9 = 2 * c5;
    pj10 = 2 * c6;
    a0 = (pn9 - pl9 / pl10 * pn10) / (pj9 - pl9 / pl10 * pj10);
    a1 = (pn9 - pj9 / pj10 * pn10) / (pl9 - pj9 / pj10 * pl10);

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  case 5:
    xc = p[0];
    tc = p[1];
    r = R_KJ_DEG_MOL;
    ag1 = r * tc * (2 * xc - 1) / (12 * xc * xc * (1 - xc) * (1 - xc));
    ag0 = (r * tc / (xc * (1 - xc)) - (12 * xc - 6) * ag1) / 2;

    /* io = 0 */
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
  case 6:
    q1 = p[0];
    q2 = p[1];
    xaly = q1;
    r = log (kb / kc);
    alpha0 = 2 * xaly - 1;
    alpha1 = 6 * xaly * (xaly - 1) + 1;
    spialy = pow ((LDBLE)10., (LDBLE)q2);
    a0 = -999.;
    a1 = -999.;
		
    if (fabs (alpha0) < tol)
			return false;
    else
    {
      azero = 1;
			
      if (fabs (alpha0) > tol)
				azero = r / alpha0;
				
      xcaly = 1 - xaly;
			
      for (i = 0; i < 50; i++)
      {
				phi1 = xcaly * xcaly * (azero + (r - azero * alpha0) * (4 * xaly - 1) / alpha1);
				phi2 = xaly * xaly * (azero + (3 - 4 * xaly) * (azero * alpha0 - r) / alpha1);
				phi1 = xaly * kb * exp (phi1);
				phi2 = xcaly * kc * exp (phi2);
				fx = phi1 + phi2 - spialy;
				fx1 = xcaly * xcaly * (1 - alpha0 * (4 * xaly - 1) / alpha1) * phi1 + xaly * xaly * (1 + alpha0 * (3 - 4 * xaly) / alpha1) * phi2;
				
				if (fabs (fx1) < 1e-10)
					return false;
				
				a0 = azero - fx / fx1;
				test = fabs (a0 - azero) + fabs (fx);
				azero = a0;
				
				if (test < tol)
					break;
      }
			
      if (i == 50)
				return false;
      else
      {
				a1 = (r - a0 * alpha0) / alpha1;

				/* io = 0 */
				ag0 = a0 * rt;
				ag1 = a1 * rt;
      }
    }
    break;
  case 7:
    ag0 = p[0];
    ag1 = p[1];
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
  case 8:
    wg2 = p[0];
    wg1 = p[1];
    ag0 = (wg2 + wg1) / 2;
    ag1 = (wg2 - wg1) / 2;
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
  case 9:
    alpha2 = p[0];
    alpha3 = p[1];
    a0 = alpha2 + 3 * alpha3 / 4;
    a1 = alpha3 / 4;
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  }

  ss->ag0 = ag0;
  ss->ag1 = ag1;
  ss->a0 = a0;
  ss->a1 = a1;
	
  return (OK);
}
//-----------------------------------------------------------------------------------------------------------