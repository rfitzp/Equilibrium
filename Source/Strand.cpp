// Strand.cpp

#include "Strand.h"

// ###########
// Constructor
// ###########
Strand::Strand ()
{
}

// ##########
// Destructor
// ##########
Strand::~Strand ()
{
}

// ####################
// Set class parameters
// ####################
void Strand::SetParameters (double _R, double _Z, double _It)
{
  R  = _R;
  Z  = _Z;
  It = _It;

  cmu = GetCoshMu (_R, _Z);
  eta = GetEta    (_R, _Z);

  if (cmu < 1. || isnan (cmu) || isnan (eta))
    {
      printf ("Strand: Error - (R, Z) = (%11.4e, %11.4e)  (z, eta/M_PI) = (%11.4e, %11.4e)\n",
	      R, Z, cmu, eta);
      exit (1);
    }
}

// #################################################
// Get weighting factor for horizontal force balance
// #################################################
double Strand::GetHorizontalFactor ()
{
  int    ntor = 1;
  double Hfac;

  Hfac  = ToroidalP (ntor, 0, cmu) - ToroidalP (ntor, 1, cmu) * cos (eta);
  Hfac *= sqrt ((cmu*cmu - 1.) /(cmu - cos (eta)));

  return It * Hfac;
}

// ###############################################
// Get weighting factor for vertical force balance
// ###############################################
double Strand::GetVerticalFactor ()
{
  int    ntor = 1;
  double Vfac;

  Vfac  = ToroidalP (ntor, 1, cmu) * sin (eta);
  Vfac *= sqrt ((cmu*cmu - 1.) /(cmu - cos (eta)));

  return It * Vfac;
}

// ###########################
// Get weighting factor for Hn
// ###########################
double Strand::GetCosFactor (double epsa, int npol)
{
  int    ntor = 1;
  double xpol = double (npol);
  double Cfac; 

  Cfac  =   ToroidalP (ntor, npol, cmu) * cos (xpol * eta) * cos (xpol * M_PI);
  Cfac *=   sqrt ((cmu*cmu - 1.) /(cmu - cos (eta)));
  Cfac *= - gsl_sf_gamma (xpol + 1.5) * pow (epsa, xpol - 1.);
  Cfac /=   sqrt(2.*M_PI) * pow (2., xpol) * gsl_sf_gamma (xpol + 1.) * (xpol*xpol - 0.25);

  return It * Cfac;
}

// ###########################
// Get weighting factor for Vn
// ###########################
double Strand::GetSinFactor (double epsa, int npol)
{
  int    ntor = 1;
  double xpol = double (npol);
  double Sfac;

  Sfac  =   ToroidalP (ntor, npol, cmu) * sin (xpol * eta) * cos (xpol * M_PI);
  Sfac *=   sqrt ((cmu*cmu - 1.) /(cmu - cos (eta)));
  Sfac *=   gsl_sf_gamma (xpol + 1.5) * pow (epsa, xpol - 1.);
  Sfac /=   sqrt(2.*M_PI) * pow (2., xpol) * gsl_sf_gamma (xpol + 1.) * (xpol*xpol - 0.25);

  return It * Sfac;
}

// ################################################################
// Function to return associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 124, Eq. (16)
// ################################################################
double Strand::ToroidalP (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("Error: ToroidalP: z < 1.: m = %3d  n = %3d  z = %11.4e\n", m, n, z);
      exit (1);
    }	

  // Calculate factor multiplying hypergeometric series
  double x  = (z - 1.) / (z + 1.);
  double d  = pow (0.5 + z/2., double (abs (n)) - 0.5);
  d        *= pow (x, double (abs (m)) /2.);
  
  if (m > 0)
    for (int j = 1; j <= m; j++)
      d *= (double (n*n - j*j + j) - 0.25) / double (j);
  else if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (j);
    }

  // Sum hypergeometric series
  double a, c, b, e, f, q, r;
  a = 0.5 - double (abs (n));
  c = 1.  + double (abs (m));
  b = a + c - 1.;
  e = 1.;
  f = 1.;
  q = 0.;
  
  do 
    {
      e *= x * (a + q) * (b + q) / (c + q) / (1. + q);
      f += e;
      q += 1.;
      
      r  = (a + q) / (a + q - 1.);
      r *= (b + q) / (b + q - 1.);
      r *= (c + q - 1.) / (c + q);
      r *= q / (q + 1.);
    }
  while (fabs (e) > 1.e-15 || fabs (r) > 1./x || q < - double (n));
  
  return f * d;
} 
  
// ####################################################################
// Function to return associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 134, Eq. (41)
//
// Added cos(m PI) factor to get agreement with EQUILIBRIUM definition.
// ####################################################################
double Strand::ToroidalQ (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("Error: ToroidalQ: z < 1.: m = %3d  n = %3d  z = %11.4e\n", m, n, z);
      exit (1);
    }

  // Calculate factor multiplying hypergeometric series
  double y  = sqrt ((z*z - 1.) /z/z);
  double d  = M_PI / sqrt (2.*z); 
  int    na = abs (n);
  if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (n*n - j*j + j) - 0.25;
    }
  int ma = abs (m);
  if (na > 0)
    for (int k = 1; k <= na; k++)
      d *= (double (ma + k) - 0.5) /2./z / double (k);
  if (ma > 0)
    for (int l = 1; l <= ma; l++)
      d *= - y * (double (l) - 0.5);

  // Sum hypergeometric series
  double a, b, c, e, f, q, r;
  a = 0.5 * (1.5 + double (ma + na));
  b = a - 0.5;
  c = double (na + 1);
  e = 1.;
  f = 1.;
  q = 0.;
  
  do
    {
      e *= (a + q) * (b + q) / (c + q) / (1. + q) / z / z;
      f += e;
      q += 1.;
      
      r  = (a + q) / (a + q - 1.);
      r *= (b + q) / (b + q - 1.);
      r *= (c + q - 1.) / (c + q);
      r *= q / (q + 1.);
    }
  while (e > 1.e-15 || r > z*z);
  
  return cos(ma*M_PI) * f * d;
} 

// ##############################################################
// Function to return hyperbolic cosine of toroidal coordinate mu
// ##############################################################
double Strand::GetCoshMu (double _R, double _Z)
{
  double d1 = sqrt ((_R + 1.) * (_R + 1.) + _Z * _Z);
  double d2 = sqrt ((_R - 1.) * (_R - 1.) + _Z * _Z);

  return 0.5 * (d1 /d2 + d2 /d1);
}

// ##########################################
// Function to return toroidal coordinate eta
// ##########################################
double Strand::GetEta (double _R, double _Z)
{
  double d1 = sqrt ((_R + 1.) * (_R + 1.) + _Z * _Z);
  double d2 = sqrt ((_R - 1.) * (_R - 1.) + _Z * _Z);

  if (_Z == 0.)
    {
      if (_R > 1.)
	return 0.;
      else if (_R < 1.)
	return M_PI;
    }
  else if (_Z >= 0.)
    return   acos ((d1*d1 + d2*d2 - 4.) /2./d1/d2);
  else
    return - acos ((d1*d1 + d2*d2 - 4.) /2./d1/d2);
}

