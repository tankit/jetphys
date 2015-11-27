#ifndef __ptresolution_h__
#define __ptresolution_h__

#include "TMath.h"

// Switch MC truth or data resolutions
bool _ismcjer = true;
//bool _jerkscale = true;//false;
bool _ak7 = true;

// From Sanmay by email 30 Mar 2015
// (direct recommendations from different eta bins => not ideal?)
const double kpar[6][2] = {
  {1.079, 0.026},
  {1.099, 0.028},
  {1.121, 0.029},
  {1.208, 0.039},
  {1.208, 0.053},
  {1.395, 0.062},
};

// Resolutions with Pythia Z2* + ak5ak7resolution12.C
// (produced on iMac desktop, with ROOT 5.30/00, iterating with AK5+2sigma)
// On Sep 15, 2014, using Winter14_V5 private pre-version (root tuples v12)
// Fit of JER for R=0.5, 8 TeV, 53X
const double vpar5[6][3] =
  {{3.13, 0.897, 0.0337},  // y 0.0-0.5, chi2 21.6/33
   {3.58, 0.868, 0.0376},  // y 0.5-1.0, chi2 12.9/33
   {4.78, 0.885, 0.0438},  // y 1.0-1.5, chi2 26.5/33
   {5.36, 0.745, 0.0265},  // y 1.5-2.0, chi2 14.5/32
   {4.45, 0.680, 0.0253},  // y 2.0-2.5, chi2 10.5/25
   {2.86, 0.874, 0.0000}}; // y 2.5-3.0, chi2 10.8/19

// Fit of JER for R=0.7, 8 TeV, 53X
const double vpar7[6][3] =
  // values from Sanmay by email 30 Mar 2015
  {{5.79356, 0.984061, 0.0290218},
   {6.10575, 0.952320, 0.0328014},
   {6.34643, 0.998053, 0.0370081},
   {6.49438, 0.866373, 0.00}, //8.97208e-07},
   {6.06794, 0.734516, 0.00}, //1.31400e-05},
   {4.57993, 0.853656, 0.00}};//1.30360e-06}};


double ptresolution(double pt, double eta) {

  int iy = min(5, int(fabs(eta) / 0.5 + 0.5));
  double res = 0;
  if (_ak7)
    res = sqrt(pow(vpar7[iy][0]/pt,2) + pow(vpar7[iy][1],2)/pt + 
	       pow(vpar7[iy][2],2));
  else
    res = sqrt(pow(vpar5[iy][0]/pt,2) + pow(vpar5[iy][1],2)/pt + 
	       pow(vpar5[iy][2],2));
  if (!_ismcjer) res *= kpar[iy][0];

  return res;
}

#endif // __ptresolution_h__
