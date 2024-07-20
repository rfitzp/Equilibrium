// Strand.h

// #########################################################################
// Class to represent single toroidal strand of poloidal magnetic field-coil
// #########################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <gsl/gsl_sf_gamma.h>

// ############
// Class header
// ############
class Strand
{
 private:

  // -----------------
  // Strand parameters
  // -----------------
  double R;   // R coordinate of strand
  double Z;   // Z coordinate of strand
  double It;  // Toroidal current flowing in strand

  double cmu; // Cosh(mu), where mu is toroidal coordinate
  double eta; // Toroidal angular coordinate

 public:

  // Constructor
  Strand ();
  // Destructor
  ~Strand ();
  
  // Set strand properties
  void SetParameters (double _R, double _Z, double _It);

  // Get weighting factor for horizontal force balance
  double GetHorizontalFactor ();
  // Get weighting factor for vertical force balance
  double GetVerticalFactor ();
  // Get weighting factor for Hn
  double GetCosFactor (double epsa, int npol);
  // Get weighting factor for Vn
  double GetSinFactor (double epsa, int npol);

 private:

  // Calculate toroidal P function
  double ToroidalP  (int m, int n, double z);
  // Calculate toroidal Q function
  double ToroidalQ  (int m, int n, double z);
  // Calculate Coshmu
  double GetCoshMu (double _R, double _Z);
  // Calculate eta
  double GetEta (double _R, double _Z);
};

