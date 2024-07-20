// Ohmic.h

// ###########################################################################
// Class to represent poloidal ohmic field-coil made up of 29 toroidal strands
// ###########################################################################

#ifndef OHMIC
#define OHMIC

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "Strand.h"

// ############
// Class header
// ############
class Ohmic
{
 private:

  // ---------------
  // Coil parameters
  // ---------------
  double R;        // R coordinate of coil centroid
  double Z;        // Z coordinate of coil top
  double dR;       // Spatial extent of coil
  double It;       // Net toroidal current flowing in coil

  Strand strand1;  // Strand 1
  Strand strand2;  // Strand 2
  Strand strand3;  // Strand 3
  Strand strand4;  // Strand 4
  Strand strand5;  // Strand 5
  Strand strand6;  // Strand 6
  Strand strand7;  // Strand 7
  Strand strand8;  // Strand 8
  Strand strand9;  // Strand 9
  Strand strand10; // strand 10
  Strand strand11; // Strand 11
  Strand strand12; // Strand 12
  Strand strand13; // Strand 13
  Strand strand14; // Strand 14
  Strand strand15; // Strand 15
  Strand strand16; // Strand 16
  Strand strand17; // Strand 17
  Strand strand18; // Strand 18
  Strand strand19; // Strand 19
  Strand strand20; // strand 20
  Strand strand21; // Strand 21
  Strand strand22; // Strand 22
  Strand strand23; // Strand 23
  Strand strand24; // Strand 24
  Strand strand25; // Strand 25
  Strand strand26; // Strand 26
  Strand strand27; // Strand 27
  Strand strand28; // strand 28
  Strand strand29; // strand 29

 public:

  // Constructor
  Ohmic ();
  // Destructor
  ~Ohmic ();
  
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

#endif //OHMIC
