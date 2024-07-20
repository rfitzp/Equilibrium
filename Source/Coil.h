// Coil.h

// ################################################################################
// Class to represent poloidal magnetic field-coil made up of five toroidal strands
// ################################################################################

#ifndef COIL
#define COIL

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "Strand.h"

// ############
// Class header
// ############
class Coil
{
 private:

  // ---------------
  // Coil parameters
  // ---------------
  double R;       // R coordinate of coil centroid
  double Z;       // Z coordinate of coil centroid
  double dR;      // Spatial extent of coil
  double It;      // Net toroidal current flowing in coil

  Strand strand1; // First strand
  Strand strand2; // Second strand
  Strand strand3; // Third strand
  Strand strand4; // Fourth strand
  Strand strand5; // Fifth strand

 public:

  // Constructor
  Coil ();
  // Destructor
  ~Coil ();
  
  // Set coil properties
  void SetParameters (double _R, double _Z, double _dR, double _It);

  // Get strand coordinates
  void GetCoordinates (double* RR, double* ZZ);
  // Get weighting factor for horizontal force balance
  double GetHorizontalFactor ();
  // Get weighting factor for vertical force balance
  double GetVerticalFactor ();
  // Get weighting factor for Hn
  double GetCosFactor (double epsa, int npol);
  // Get weighting factor for Vn
  double GetSinFactor (double epsa, int npol);

 private:
};

#endif //COIL
