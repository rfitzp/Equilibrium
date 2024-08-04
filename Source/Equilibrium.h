// Equilibrium.h

// ########################################################################################
// Class to solve inverse aspect-ratio expanded tokamak equilibrium problem

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field).
// Radial coordinate, r, normalized to epsa * R_0, where eps_a is inverse-aspect ratio.
// So r = 0 is magnetic axis, and r = 1 is plasma/vacuum interface.

// Flux-surfaces:

// R(r,w) = 1 - epsa r cosw + epsa^2 H1(r) + epsa^2 sum_{n=2,Ns} [Hn(r) cos(n-1)w + Vn(r) sin(n-1)w]
// Z(r,w) =     epsa r sinw                + epsa^2 sum_{n=2,Ns} [Hn(r) sin(n-1)w - Vn(r) cos(n-1)w]
//
// Here, R, phi, Z are cylindrical polar coordinates while r, w, phi are flux coordinates

// Edge shaping: Hna = Hn(1), Vna = Vn(1)

// Equilibrium profiles:
//
//  Lowest order (i.e., cylindrical) safety factor profile is q0(r) = r^2 /f1(r)
//  Pressure profile is P(r) = epsa^2 p2(r)
//
// f1 = (1/nu/qc) [1 - (1 - r^2)^nu]
//
// p2 = pc (1 - r^2)^mu
//
// qc is safety-factor on magnetic axis.
// nu * qc is lowest-order safety-factor at plasma/vacuum interface.

// Inputs:
//  Inputs/Equilibrium.json - JSON file

// Outputs:
//  Plots/Equilibrium.nc

// Ploting scripts:
//  Plots/*.py

// Class uses following external libraries:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Netcdf-c library       (https://github.com/Unidata/netcdf-c)

// Author:
// Richard Fitzpatrick,
// Institute of Fusion Studies,
// Department of Physics,
// University of Texas at Austin
// rfitzp@utexas.edu

// Source: https://github.com/rfitzp/Equilibrium/

// Documentation: ../Documentation/Equilibrium.pdf

// ########################################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <nlohmann/json.hpp>
#include <netcdf.h>

#include "Coil.h"
#include "Ohmic.h"

using namespace blitz;
using           json = nlohmann::json;

// ############
// Class header
// ############
class Equilibrium
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  double qc;          // Safety-factor on magnetic axis (read from JSON file)
  double nu;          // Current peaking parameter (read from JSON file)
  double pc;          // Normalized plasma pressure on magnetic axis (read from JSON file)
  double mu;          // Pressure peaking parameter (read from JSON file)
  double epsa;        // Inverse aspect-ratio of plasma (read from JSON file)
  vector<double> Hna; // H2(1), H3(1), etc (read from JSON file)
  vector<double> Vna; // V2(1), V3(1), etc (read from JSON file)
  
  double Pped;   // Normalized height of pedestal (read from JSON file)
  double wped;   // Normalized width of pedestal (read from JSON file)
  double rped;   // Normalized radius of pedestal (read from JSON file)

  int    CFLAG;  // Flag for inclusion of external coils (read from JSON file)
                 //  If not set then shaping data determined from Hna and Vna values in JSON file

  double Hshift; // Horizontal shift of magnetic axis (read from JSON file)
  double Vshift; // Search for vertical shift of magnetic axis
                 //  takes place for shifts in range -Vshift to +Vshift (read from JSON file)

  Ohmic  ohmic;  // Ohmic heating coil
  double Ro;     // R coordinate of ohmic heating coil (read from JSON file)
  double Zo;     // Height of ohmic heating coil (read from JSON file)
  double dRo;    // Width of ohmic heating coil (read from JSON file)
  double Wo;     // Relative current in ohmic heating coil (read from JSON file)
  double Ito;    // Actual current in ohmic heating coil

  Coil   coil1;  // Poloidal field-coil 1
  double R1;     // R coordinate of poloidal field-coil 1 (read from JSON file)
  double Z1;     // Z coordinate of poloidal field-coil 1 (read from JSON file)
  double dR1;    // Size of poloidal field-coil 1 (read from JSON file)
  double W1;     // Relative current in poloidal field-coil 1 (read from JSON file)
  double It1;    // Actual current in poloidal field-coil 1 

  Coil   coil2;  // Poloidal field-coil 2
  double R2;     // R coordinate of poloidal field-coil 2 (read from JSON file)
  double Z2;     // Z coordinate of poloidal field-coil 2 (read from JSON file)
  double dR2;    // Size of poloidal field-coil 2 (read from JSON file)
  double W2;     // Relative current in poloidal field-coil 2 (read from JSON file)
  double It2;    // Actual current in poloidal field-coil 2

  Coil   coil3;  // Poloidal field-coil 3
  double R3;     // R coordinate of poloidal field-coil 3 (read from JSON file)
  double Z3;     // Z coordinate of poloidal field-coil 3 (read from JSON file)
  double dR3;    // Size of poloidal field-coil 3 (read from JSON file)
  double W3;     // Relative current in poloidal field-coil 3 (read from JSON file)
  double It3;    // Actual current in poloidal field-coil 3

  Coil   coil4;  // Poloidal field-coil 4
  double R4;     // R coordinate of poloidal field-coil 4 (read from JSON file)
  double Z4;     // Z coordinate of poloidal field-coil 4 (read from JSON file)
  double dR4;    // Size of poloidal field-coil 4 (read from JSON file)
  double W4;     // Relative current in poloidal field-coil 4 (read from JSON file)
  double It4;    // Actual current in poloidal field-coil 4

  Coil   coil5;  // Poloidal field-coil 5
  double R5;     // R coordinate of poloidal field-coil 5 (read from JSON file)
  double Z5;     // Z coordinate of poloidal field-coil 5 (read from JSON file)
  double dR5;    // Size of poloidal field-coil 5 (read from JSON file)
  double W5;     // Relative current in poloidal field-coil 5 (read from JSON file)
  double It5;    // Actual current in poloidal field-coil 5

  Coil   coil6;  // Poloidal field-coil 6
  double R6;     // R coordinate of poloidal field-coil 6 (read from JSON file)
  double Z6;     // Z coordinate of poloidal field-coil 6 (read from JSON file)
  double dR6;    // Size of poloidal field-coil 6 (read from JSON file)
  double W6;     // Relative current in poloidal field-coil 6 (read from JSON file)
  double It6;    // Actual current in poloidal field-coil 6
  
  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis
  int    Ns;     // Number of shaping harmonics
  int    Nr;     // Number of radial grid points
  int    Nf;     // Number of magnetic flux-surfaces
  int    Nw;     // Number of angular points on magnetic flux-surfaces

  // ----------------
  // Calculation data
  // ----------------
  double* rr;               // Radial grid points

  double* p2;               // Plasma pressure
  double* f1;               // Lowest-order poloidal flux function
  double* f3;               // Higher-order poloidal flux function
  double* g2;               // Toroidal flux function
  double* q0;               // Lowest-order safety-factor
  double* q2;               // Higher-order safety-factor
  double* It;               // Toroidal plasma current
  double* Ip;               // Poloidal plasma current
  double* Jt;               // Radial derivative of toroidal plasma current
  double* Jp;               // Radial derivative of poloidal plasma current

  double* pp2;              // Pedestal contribution to pressure
  double* ff1;              // Pedestal contribution to f1
  
  Array<double,2> HHfunc;   // Horizontal shaping functions
  Array<double,2> VVfunc;   // Vertical shaping functions
  Array<double,2> HPfunc;   // Radial derivatives of horizontal shaping functions
  Array<double,2> VPfunc;   // Radial derivatives of vertical shaping functions

  gsl_spline* Itspline;     // Interpolated It function
  gsl_spline* Ipspline;     // Interpolated Ip function

  gsl_interp_accel* Itacc;  // Accelerator for interpolated It function
  gsl_interp_accel* Ipacc;  // Accelerator for interpolated Ip function

  gsl_spline* pp2spline;    // Interpolated pp2 function
  gsl_spline* ff1spline;    // Interpolated ff1 function

  gsl_interp_accel* pp2acc; // Accelerator for interpolated pp2 function
  gsl_interp_accel* ff1acc; // Accelerator for interpolated ff1 function

  gsl_spline*  g2spline;    // Interpolated g2 function
  gsl_spline** HHspline;    // Interpolated horizontal shaping functions
  gsl_spline** VVspline;    // Interpolated vertical shaping functions
  gsl_spline** HPspline;    // Interpolated radial derivatives of horizontal shaping functions
  gsl_spline** VPspline;    // Interpolated radial derivatives of vertical shaping functions

  gsl_interp_accel*  g2acc; // Accelerator for interpolated g2 function
  gsl_interp_accel** HHacc; // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc; // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc; // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc; // Accelerator for interpolated radial derivatives of vertical shaping functions

  Array<double,2> RR;       // R coodinates of magnetic flux-surfaces
  Array<double,2> ZZ;       // Z coodinates of magnetic flux-surfaces
  
  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;     // Integration accuracy (read from JSON file)
  double h0;      // Initial step-length (read from JSON file)
  double hmin;    // Minimum step-length (read from JSON file)
  double hmax;    // Maximum step-length (read from JSON file)
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Integration error calcualation flag

  // -----------------------
  // Root finding parameters
  // -----------------------
  int    nint;    // Number of search intervals
  double Eta;     // Min. magnitude of f at root f(x) = 0
  int    Maxiter; // Maximum number of iterations
  
  // ------------------
  // RK4/RK5 parameters
  // ------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int count, rhs_chooser, flg;

 public:

  // Constructor
  Equilibrium ();
  // Destructor
  ~Equilibrium ();

  // Solve problem
  void Solve ();

private:

  // Return f1(r)
  double Getf1 (double r);
  // Return f1'(r)
  double Getf1p (double r);
  // Return p2(r)
  double Getp2 (double r);
  // Return p2'(r)
  double Getp2p (double r);

  // Evaluate right-hand sides of differential equations
  void Rhs (double x, double* y, double* dydx);

  // Adaptive step length RK4/RK5 integration routine
  void RK4RK5Adaptive (int neqns, double& x, double* y, double& h, 
		       double& t_err, double acc, double S, double T, int& rept,
		       int maxrept, double h_min, double h_max, int flag, 
		       int diag, FILE* file);
  // Fixed step length RK4/RK5 integration routine
  void RK4RK5Fixed (int neqns, double& x, double* y, double* err, double h);

  // Target functions for zero finding
  double Feval (double);
  // Zero finding routine
  double RootFind (double, double);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double, double, double, double, double&);

  // Read JSON file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);
};
