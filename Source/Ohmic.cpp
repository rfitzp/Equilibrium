// Ohmic.cpp

#include "Ohmic.h"

// ###########
// Constructor
// ###########
Ohmic::Ohmic ()
{
}

// ##########
// Destructor
// ##########
Ohmic::~Ohmic ()
{
}

// ####################
// Set class parameters
// ####################
void Ohmic::SetParameters (double _R, double _Z, double _dR, double _It)
{
  R  = _R;
  Z  = _Z;
  dR = _dR;
  It = _It;

  strand1 .SetParameters (_R - _dR/2.,   1.*_Z/9., _It/29.);
  strand2 .SetParameters (_R - _dR/2.,   3.*_Z/9., _It/29.);
  strand3 .SetParameters (_R - _dR/2.,   5.*_Z/9., _It/29.);
  strand4 .SetParameters (_R - _dR/2.,   7.*_Z/9., _It/29.);
  strand5 .SetParameters (_R - _dR/2.,   9.*_Z/9., _It/29.);
  strand6 .SetParameters (_R - _dR/2., - 1.*_Z/9., _It/29.);
  strand7 .SetParameters (_R - _dR/2., - 3.*_Z/9., _It/29.);
  strand8 .SetParameters (_R - _dR/2., - 5.*_Z/9., _It/29.);
  strand9 .SetParameters (_R - _dR/2., - 7.*_Z/9., _It/29.);
  strand10.SetParameters (_R - _dR/2., - 9.*_Z/9., _It/29.);
  strand11.SetParameters (_R + _dR/2.,   1.*_Z/9., _It/29.);
  strand12.SetParameters (_R + _dR/2.,   3.*_Z/9., _It/29.);
  strand13.SetParameters (_R + _dR/2.,   5.*_Z/9., _It/29.);
  strand14.SetParameters (_R + _dR/2.,   7.*_Z/9., _It/29.);
  strand15.SetParameters (_R + _dR/2.,   9.*_Z/9., _It/29.);
  strand16.SetParameters (_R + _dR/2., - 1.*_Z/9., _It/29.);
  strand17.SetParameters (_R + _dR/2., - 3.*_Z/9., _It/29.);
  strand18.SetParameters (_R + _dR/2., - 5.*_Z/9., _It/29.);
  strand19.SetParameters (_R + _dR/2., - 7.*_Z/9., _It/29.);
  strand20.SetParameters (_R + _dR/2., - 9.*_Z/9., _It/29.);
  strand21.SetParameters (_R,            0.*_Z/9., _It/29.);
  strand22.SetParameters (_R,            2.*_Z/9., _It/29.);
  strand23.SetParameters (_R,            4.*_Z/9., _It/29.);
  strand24.SetParameters (_R,            6.*_Z/9., _It/29.);
  strand25.SetParameters (_R,            8.*_Z/9., _It/29.);
  strand26.SetParameters (_R,          - 2.*_Z/9., _It/29.);
  strand27.SetParameters (_R,          - 4.*_Z/9., _It/29.);
  strand28.SetParameters (_R,          - 6.*_Z/9., _It/29.);
  strand29.SetParameters (_R,          - 8.*_Z/9., _It/29.);
}

// ######################
// Get strand coordinates
// ######################
void Ohmic::GetCoordinates (double* RR, double* ZZ)
{
  RR[0]  = R - dR/2.; ZZ[0]  =   1.*Z/9.;
  RR[1]  = R - dR/2.; ZZ[1]  =   3.*Z/9.;
  RR[2]  = R - dR/2.; ZZ[2]  =   5.*Z/9.;
  RR[3]  = R - dR/2.; ZZ[3]  =   7.*Z/9.;
  RR[4]  = R - dR/2.; ZZ[4]  =   9.*Z/9.;
  RR[5]  = R - dR/2.; ZZ[5]  = - 1.*Z/9.;
  RR[6]  = R - dR/2.; ZZ[6]  = - 3.*Z/9.;
  RR[7]  = R - dR/2.; ZZ[7]  = - 5.*Z/9.;
  RR[8]  = R - dR/2.; ZZ[8]  = - 7.*Z/9.;
  RR[9]  = R - dR/2.; ZZ[9]  = - 9.*Z/9.;
  RR[10] = R + dR/2.; ZZ[10] =   1.*Z/9.;
  RR[11] = R + dR/2.; ZZ[11] =   3.*Z/9.;
  RR[12] = R + dR/2.; ZZ[12] =   5.*Z/9.;
  RR[13] = R + dR/2.; ZZ[13] =   7.*Z/9.;
  RR[14] = R + dR/2.; ZZ[14] =   9.*Z/9.;
  RR[15] = R + dR/2.; ZZ[15] = - 1.*Z/9.;
  RR[16] = R + dR/2.; ZZ[16] = - 3.*Z/9.;
  RR[17] = R + dR/2.; ZZ[17] = - 5.*Z/9.;
  RR[18] = R + dR/2.; ZZ[18] = - 7.*Z/9.;
  RR[19] = R + dR/2.; ZZ[19] = - 9.*Z/9.;
  RR[20] = R;         ZZ[20] =   0.*Z/9.;
  RR[21] = R;         ZZ[21] =   2.*Z/9.;
  RR[22] = R;         ZZ[22] =   4.*Z/9.;
  RR[23] = R;         ZZ[23] =   6.*Z/9.;
  RR[24] = R;         ZZ[24] =   8.*Z/9.;
  RR[25] = R;         ZZ[25] = - 2.*Z/9.;
  RR[26] = R;         ZZ[26] = - 4.*Z/9.;
  RR[27] = R;         ZZ[27] = - 6.*Z/9.;
  RR[28] = R;         ZZ[28] = - 8.*Z/9.;

}

// #################################################
// Get weighting factor for horizontal force balance
// #################################################
double Ohmic::GetHorizontalFactor ()
{
  return
    strand1. GetHorizontalFactor () +
    strand2. GetHorizontalFactor () +
    strand3. GetHorizontalFactor () +
    strand4. GetHorizontalFactor () +
    strand5. GetHorizontalFactor () +
    strand6. GetHorizontalFactor () +
    strand7. GetHorizontalFactor () +
    strand8. GetHorizontalFactor () +
    strand9. GetHorizontalFactor () +
    strand10.GetHorizontalFactor () +
    strand11.GetHorizontalFactor () +
    strand12.GetHorizontalFactor () +
    strand13.GetHorizontalFactor () +
    strand14.GetHorizontalFactor () +
    strand15.GetHorizontalFactor () +
    strand16.GetHorizontalFactor () +
    strand17.GetHorizontalFactor () +
    strand18.GetHorizontalFactor () +
    strand19.GetHorizontalFactor () +
    strand20.GetHorizontalFactor () +
    strand21.GetHorizontalFactor () +
    strand22.GetHorizontalFactor () +
    strand23.GetHorizontalFactor () +
    strand24.GetHorizontalFactor () +
    strand25.GetHorizontalFactor () +
    strand26.GetHorizontalFactor () +
    strand27.GetHorizontalFactor () +
    strand28.GetHorizontalFactor () +
    strand29.GetHorizontalFactor ();
}

// ###############################################
// Get weighting factor for vertical force balance
// ###############################################
double Ohmic::GetVerticalFactor ()
{
  return
    strand1. GetVerticalFactor () +
    strand2. GetVerticalFactor () +
    strand3. GetVerticalFactor () +
    strand4. GetVerticalFactor () +
    strand5. GetVerticalFactor () +
    strand6. GetVerticalFactor () +
    strand7. GetVerticalFactor () +
    strand8. GetVerticalFactor () +
    strand9. GetVerticalFactor () +
    strand10.GetVerticalFactor () +
    strand11.GetVerticalFactor () +
    strand12.GetVerticalFactor () +
    strand13.GetVerticalFactor () +
    strand14.GetVerticalFactor () +
    strand15.GetVerticalFactor () +
    strand16.GetVerticalFactor () +
    strand17.GetVerticalFactor () +
    strand18.GetVerticalFactor () +
    strand19.GetVerticalFactor () +
    strand20.GetVerticalFactor () +
    strand21.GetVerticalFactor () +
    strand22.GetVerticalFactor () +
    strand23.GetVerticalFactor () +
    strand24.GetVerticalFactor () +
    strand25.GetVerticalFactor () +
    strand26.GetVerticalFactor () +
    strand27.GetVerticalFactor () +
    strand28.GetVerticalFactor () +
    strand29.GetVerticalFactor ();
}

// ###########################
// Get weighting factor for Hn
// ###########################
double Ohmic::GetCosFactor (double epsa, int npol)
{
  return
    strand1. GetCosFactor (epsa, npol) +
    strand2. GetCosFactor (epsa, npol) +
    strand3. GetCosFactor (epsa, npol) +
    strand4. GetCosFactor (epsa, npol) +
    strand5. GetCosFactor (epsa, npol) +
    strand6. GetCosFactor (epsa, npol) +
    strand7. GetCosFactor (epsa, npol) +
    strand8. GetCosFactor (epsa, npol) +
    strand9. GetCosFactor (epsa, npol) +
    strand10.GetCosFactor (epsa, npol) +
    strand11.GetCosFactor (epsa, npol) +
    strand12.GetCosFactor (epsa, npol) +
    strand13.GetCosFactor (epsa, npol) +
    strand14.GetCosFactor (epsa, npol) +
    strand15.GetCosFactor (epsa, npol) +
    strand16.GetCosFactor (epsa, npol) +
    strand17.GetCosFactor (epsa, npol) +
    strand18.GetCosFactor (epsa, npol) +
    strand19.GetCosFactor (epsa, npol) +
    strand20.GetCosFactor (epsa, npol) +
    strand21.GetCosFactor (epsa, npol) +
    strand22.GetCosFactor (epsa, npol) +
    strand23.GetCosFactor (epsa, npol) +
    strand24.GetCosFactor (epsa, npol) +
    strand25.GetCosFactor (epsa, npol) +
    strand26.GetCosFactor (epsa, npol) +
    strand27.GetCosFactor (epsa, npol) +
    strand28.GetCosFactor (epsa, npol) +
    strand29.GetCosFactor (epsa, npol);
  }

// ###########################
// Get weighting factor for Vn
// ###########################
double Ohmic::GetSinFactor (double epsa, int npol)
{
  return
    strand1. GetSinFactor (epsa, npol) +
    strand2. GetSinFactor (epsa, npol) +
    strand3. GetSinFactor (epsa, npol) +
    strand4. GetSinFactor (epsa, npol) +
    strand5. GetSinFactor (epsa, npol) +
    strand6. GetSinFactor (epsa, npol) +
    strand7. GetSinFactor (epsa, npol) +
    strand8. GetSinFactor (epsa, npol) +
    strand9. GetSinFactor (epsa, npol) +
    strand10.GetSinFactor (epsa, npol) +
    strand11.GetSinFactor (epsa, npol) +
    strand12.GetSinFactor (epsa, npol) +
    strand13.GetSinFactor (epsa, npol) +
    strand14.GetSinFactor (epsa, npol) +
    strand15.GetSinFactor (epsa, npol) +
    strand16.GetSinFactor (epsa, npol) +
    strand17.GetSinFactor (epsa, npol) +
    strand18.GetSinFactor (epsa, npol) +
    strand19.GetSinFactor (epsa, npol) +
    strand20.GetSinFactor (epsa, npol) +
    strand21.GetSinFactor (epsa, npol) +
    strand22.GetSinFactor (epsa, npol) +
    strand23.GetSinFactor (epsa, npol) +
    strand24.GetSinFactor (epsa, npol) +
    strand25.GetSinFactor (epsa, npol) +
    strand26.GetSinFactor (epsa, npol) +
    strand27.GetSinFactor (epsa, npol) +
    strand28.GetSinFactor (epsa, npol) +
    strand29.GetSinFactor (epsa, npol);
}

