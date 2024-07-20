// Coil.cpp

#include "Coil.h"

// ###########
// Constructor
// ###########
Coil::Coil ()
{
}

// ##########
// Destructor
// ##########
Coil::~Coil ()
{
}

// ####################
// Set class parameters
// ####################
void Coil::SetParameters (double _R, double _Z, double _dR, double _It)
{
  R  = _R;
  Z  = _Z;
  dR = _dR;
  It = _It;

  strand1.SetParameters (_R,          _Z,          _It/5.);
  strand2.SetParameters (_R + _dR/2., _Z + _dR/2., _It/5.);
  strand3.SetParameters (_R - _dR/2., _Z + _dR/2., _It/5.);
  strand4.SetParameters (_R + _dR/2., _Z - _dR/2., _It/5.);
  strand5.SetParameters (_R - _dR/2., _Z - _dR/2., _It/5.);
}

// ######################
// Get strand coordinates
// ######################
void Coil::GetCoordinates (double* RR, double* ZZ)
{
  RR[0] = R;         ZZ[0] = Z;
  RR[1] = R + dR/2.; ZZ[1] = Z + dR/2.;
  RR[2] = R - dR/2.; ZZ[2] = Z + dR/2.;
  RR[3] = R + dR/2.; ZZ[3] = Z - dR/2.;
  RR[4] = R - dR/2.; ZZ[4] = Z - dR/2.;
}

// #################################################
// Get weighting factor for horizontal force balance
// #################################################
double Coil::GetHorizontalFactor ()
{
  return
    strand1.GetHorizontalFactor () +
    strand2.GetHorizontalFactor () +
    strand3.GetHorizontalFactor () +
    strand4.GetHorizontalFactor () +
    strand5.GetHorizontalFactor ();
}

// ###############################################
// Get weighting factor for vertical force balance
// ###############################################
double Coil::GetVerticalFactor ()
{
  return
    strand1.GetVerticalFactor () +
    strand2.GetVerticalFactor () +
    strand3.GetVerticalFactor () +
    strand4.GetVerticalFactor () +
    strand5.GetVerticalFactor ();
}

// ###########################
// Get weighting factor for Hn
// ###########################
double Coil::GetCosFactor (double epsa, int npol)
{
  return
    strand1.GetCosFactor (epsa, npol) +
    strand2.GetCosFactor (epsa, npol) +
    strand3.GetCosFactor (epsa, npol) +
    strand4.GetCosFactor (epsa, npol) +
    strand5.GetCosFactor (epsa, npol);
}

// ###########################
// Get weighting factor for Vn
// ###########################
double Coil::GetSinFactor (double epsa, int npol)
{
  return
    strand1.GetSinFactor (epsa, npol) +
    strand2.GetSinFactor (epsa, npol) +
    strand3.GetSinFactor (epsa, npol) +
    strand4.GetSinFactor (epsa, npol) +
    strand5.GetSinFactor (epsa, npol); 
}

