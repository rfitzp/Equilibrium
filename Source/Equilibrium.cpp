// Equilibrium.cpp

#include "Equilibrium.h"

// ###########
// Constructor
// ###########
Equilibrium::Equilibrium ()
{
  // -----------------------------------
  // Set adaptive integration parameters
  // -----------------------------------
  maxrept = 50;
  flag    = 2;

  // ---------------------------
  // Set root finding parameters
  // ---------------------------
  nint    = 1;
  Eta     = 1.e-12;
  Maxiter = 60;
  
  // ----------------------
  // Set RK4/RK5 parameters
  // ----------------------
  aa1  = 0.;
  aa2  = 1./5.;
  aa3  = 3./10.;
  aa4  = 3./5.;
  aa5  = 1.;
  aa6  = 7./8.;

  cc1  =  37./378.;
  cc3  = 250./621.;
  cc4  = 125./594.;
  cc6  = 512./1771.;

  ca1  = cc1 -  2825./27648.;
  ca3  = cc3 - 18575./48384.;
  ca4  = cc4 - 13525./55296.;
  ca5  =     -   277./14336.;
  ca6  = cc6 -     1./4.;

  bb21 = 1./5.;

  bb31 = 3./40.;
  bb32 = 9./40.;

  bb41 =   3./10.;
  bb42 = - 9./10.;
  bb43 =   6./5.;

  bb51 = - 11./54.;
  bb52 =    5./2.;
  bb53 = - 70./27.;
  bb54 =   35./27.;

  bb61 =  1631./55296.;
  bb62 =   175./512.;
  bb63 =   575./13824.;
  bb64 = 44275./110592.;
  bb65 =   253./4096.;

  // -----------------------------------------
  // Read namelist file Inputs/Equilibrium.nml
  // -----------------------------------------
  NameListRead (&qc, &nu, &pc, &mu, &epsa,
		&CFLAG, 
		&Hshift, &Vshift,
		&Ro, &Zo, &dRo, &Wo,
		&R1, &Z1, &dR1, &W1,
		&R2, &Z2, &dR2, &W2,
		&R3, &Z3, &dR3, &W3,
		&R4, &Z4, &dR4, &W4,
		&R5, &Z5, &dR5, &W5,
		&R6, &Z6, &dR6, &W6,
		&Pped, &wped, &rped,
		&eps, &Ns, &Nr, &Nf, &Nw,
		&acc, &h0, &hmin, &hmax);

  // ............
  // Sanity check
  // ............
  if (qc <= 0.)
    {
      printf ("Equilibrium:: Error - qc must be positive\n");
      exit (1);
    }
   if (nu <= 0.)
    {
      printf ("Equilibrium:: Error -  nu must be positive\n");
      exit (1);
    }
   if (pc < 0.)
    {
      printf ("Equilibrium:: Error - pc cannot be negative\n");
      exit (1);
    }
   if (epsa <= 0.)
    {
      printf ("Equilibrium:: Error - epsa must be positive\n");
      exit (1);
    }
   if (Ro <= 0.)
    {
      printf ("Equilibrium:: Error - Ro must be positive\n");
      exit (1);
    }
   if (dRo <= 0.)
    {
      printf ("Equilibrium:: Error - dRo must be positive\n");
      exit (1);
    }
    if (R1 <= 0.)
    {
      printf ("Equilibrium:: Error - R1 must be positive\n");
      exit (1);
    }
   if (dR1 <= 0.)
    {
      printf ("Equilibrium:: Error - dR1 must be positive\n");
      exit (1);
    }
    if (R2 <= 0.)
    {
      printf ("Equilibrium:: Error - R2 must be positive\n");
      exit (1);
    }
   if (dR2 <= 0.)
    {
      printf ("Equilibrium:: Error - dR2 must be positive\n");
      exit (1);
    }
    if (R3 <= 0.)
    {
      printf ("Equilibrium:: Error - R3 must be positive\n");
      exit (1);
    }
   if (dR3 <= 0.)
    {
      printf ("Equilibrium:: Error - dR3 must be positive\n");
      exit (1);
    }
    if (R4 <= 0.)
    {
      printf ("Equilibrium:: Error - R4 must be positive\n");
      exit (1);
    }
   if (dR4 <= 0.)
    {
      printf ("Equilibrium:: Error - dR4 must be positive\n");
      exit (1);
    }
    if (R5 <= 0.)
    {
      printf ("Equilibrium:: Error - R5 must be positive\n");
      exit (1);
    }
   if (dR5 <= 0.)
    {
      printf ("Equilibrium:: Error - dR5 must be positive\n");
      exit (1);
    }
    if (R6 <= 0.)
    {
      printf ("Equilibrium:: Error - R6 must be positive\n");
      exit (1);
    }
   if (dR6 <= 0.)
    {
      printf ("Equilibrium:: Error - dR6 must be positive\n");
      exit (1);
    }
   if (Pped < 0.)
    {
      printf ("Equilibrium:: Error - Pped cannot be negative\n");
      exit (1);
    }
   if (wped <= 0.)
     {
       printf ("Equilibrium:: Error - wped must be positive\n");
       exit (1);
    }
   if (rped <= 0.)
     {
       printf ("Equilibrium:: Error - rped must be positive\n");
       exit (1);
    }
   if (eps <= 0.)
     {
       printf ("Equilibrium:: Error - eps must be positive\n");
       exit (1);
    }
   if (Ns <= 0)
     {
       printf ("Equilibrium:: Error - eps must be positive\n");
       exit (1);
    }
   if (Ns < 1)
    {
      printf ("Equilibrium:: Error - Ns cannot be less than unity\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("Equilibrium:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (Nf < 2)
    {
      printf ("Equilibrium:: Error - Nf cannot be less than two\n");
      exit (1);
    }
  if (Nw < 2)
    {
      printf ("Equilibrium:: Error - Nw cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("Equilibrium:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("Equilibrium:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("Equilibrium:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("Equilibrium:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("Equilibrium:: Error - hmax must exceed hmin\n");
      exit (1);
    }
}

// ###########
// Destructor
// ###########
Equilibrium::~Equilibrium ()
{
}

// #########################
// Function to solve problem
// #########################
void Equilibrium::Solve ()
{
  // ...............
  // Allocate memory
  // ...............
  rr  = new double[Nr+1];

  p2  = new double[Nr+1];
  f1  = new double[Nr+1];
  f3  = new double[Nr+1];
  g2  = new double[Nr+1];
  q0  = new double[Nr+1];
  q2  = new double[Nr+1];
  It  = new double[Nr+1];
  Ip  = new double[Nr+1];
  Jt  = new double[Nr+1];
  Jp  = new double[Nr+1];

  pp2 = new double[Nr+1];
  ff1 = new double[Nr+1];

  pp2spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  ff1spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

  pp2acc    = gsl_interp_accel_alloc ();
  ff1acc    = gsl_interp_accel_alloc ();

  HHfunc.resize (Ns+1, Nr+1);
  VVfunc.resize (Ns+1, Nr+1);
  HPfunc.resize (Ns+1, Nr+1);
  VPfunc.resize (Ns+1, Nr+1);

  g2spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Itspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Ipspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

  g2acc    = gsl_interp_accel_alloc ();
  Itacc    = gsl_interp_accel_alloc ();
  Ipacc    = gsl_interp_accel_alloc ();

  HHspline = new gsl_spline* [Ns+1];
  VVspline = new gsl_spline* [Ns+1];
  HPspline = new gsl_spline* [Ns+1];
  VPspline = new gsl_spline* [Ns+1];

  HHacc    = new gsl_interp_accel* [Ns+1];
  VVacc    = new gsl_interp_accel* [Ns+1];
  HPacc    = new gsl_interp_accel* [Ns+1];
  VPacc    = new gsl_interp_accel* [Ns+1];

  for (int n = 0; n <= Ns; n++)
    {
      HHspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VVspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      HPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

      HHacc[n]    = gsl_interp_accel_alloc ();
      VVacc[n]    = gsl_interp_accel_alloc ();
      HPacc[n]    = gsl_interp_accel_alloc ();
      VPacc[n]    = gsl_interp_accel_alloc ();
    }

  RR.resize (Nf, Nr+1);
  ZZ.resize (Nf, Nr+1);
  
  // ..................
  // Set up radial grid
  // ..................
  flg = 0;
  for (int i = 0; i <= Nr; i++)
    {
      rr[i] = double (i) /double (Nr);
      p2[i] = Getp2 (rr[i]);
      f1[i] = Getf1 (rr[i]);
    }

  // ...............................
  // Calculate pp2 and ff1 functions
  // ...............................
  double  r, h, t_err;
  int     rept, neqns = 2;
  double* y0   = new double[neqns];
  double* err0 = new double[neqns];
  rhs_chooser  = 2;

  pp2[0] = 0.;
  ff1[0] = 0.;

  count  = 0;
  h      = h0;
  r      = eps;
  y0[0]  = 0.;
  y0[1]  = 0.;

  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  RK4RK5Adaptive(neqns, r, y0, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      RK4RK5Fixed(neqns, r, y0, err0, rr[i] - r);

      pp2[i] = y0[0];
      ff1[i] = y0[1];
    }

  delete[] y0; delete[] err0;

  double pp2a = pp2[Nr];
  for (int i = 0; i <= Nr; i++)
    {
      pp2[i]  = pp2[i] - pp2a;
      p2 [i] += pp2[i];
      f1 [i]  = sqrt(f1[i]*f1[i] + ff1[i]);
    }
  
  gsl_spline_init (pp2spline, rr, p2,  Nr+1);
  gsl_spline_init (ff1spline, rr, f1,  Nr+1);
  flg = 1;

  // ....................................
  // Integrate shaping function equations
  // ....................................
  neqns        = 1 + 2 + 4 * (Ns - 1);
  double* y    = new double[neqns];
  double* err  = new double[neqns];
  rhs_chooser  = 0;

  count        = 0;
  h            = h0;
  r            = rr[0];
  g2[0]        = 0.;
  HHfunc(1, 0) = 0.;
  HPfunc(1, 0) = 0.;
  int j        = 3;
  for (int n = 2; n <= Ns; n++)
    {
      if (n == 2)
	{
	  HHfunc(n, 0) = 0.;
	  HPfunc(n, 0) = 1.;
	  VVfunc(n, 0) = 0.;
	  VPfunc(n, 0) = 1.;
	}
      else
	{
	  HHfunc(n, 0) = 0.;
	  HPfunc(n, 0) = 0.;
	  VVfunc(n, 0) = 0.;
	  VPfunc(n, 0) = 0.;
	}
    }

  double f1c   = 1./qc;
  double p2ppc = - 2.*mu*pc;

  r    = eps;
  y[0] = - (f1c*f1c + p2ppc/2.)  *r*r;
  y[1] = (2.*p2ppc/f1c/f1c - 1.) *r*r/8.;
  y[2] = (2.*p2ppc/f1c/f1c - 1.) *r/4.;
  j    = 3;
  for (int n = 2; n <= Ns; n++)
    {
      y[j] =                pow (r, double (n-1)); j++;
      y[j] = double (n-1) * pow (r, double (n-2)); j++;
      y[j] =                pow (r, double (n-1)); j++;
      y[j] = double (n-1) * pow (r, double (n-2)); j++;
    }

  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  RK4RK5Adaptive(neqns, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      RK4RK5Fixed(neqns, r, y, err, rr[i] - r);

      g2[i]        = y[0];
      HHfunc(1, i) = y[1];
      HPfunc(1, i) = y[2];
      int j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  HHfunc(n, i) = y[j]; j++;
	  HPfunc(n, i) = y[j]; j++;
	  VVfunc(n, i) = y[j]; j++;
	  VPfunc(n, i) = y[j]; j++;
	}
    }

  delete[] y; delete [] err;

  if (CFLAG)
    {
      // .............................
      // Find vertical shift of plasma
      // .............................
      double vshift = RootFind (- Vshift, Vshift);
      
      // ....................................
      // Set up external magnetic field-coils
      // ....................................
      ohmic.SetParameters (Ro - Hshift, Zo - vshift, dRo, Wo);
      coil1.SetParameters (R1 - Hshift, Z1 - vshift, dR1, W1);
      coil2.SetParameters (R2 - Hshift, Z2 - vshift, dR2, W2);
      coil3.SetParameters (R3 - Hshift, Z3 - vshift, dR3, W3);
      coil4.SetParameters (R4 - Hshift, Z4 - vshift, dR4, W4);
      coil5.SetParameters (R5 - Hshift, Z5 - vshift, dR5, W5);
      coil6.SetParameters (R6 - Hshift, Z6 - vshift, dR6, W6);
      
      // ...........................
      // Set vertical magnetic field
      // ...........................
      double Hfaco = ohmic.GetHorizontalFactor (); 
      double Hfac1 = coil1.GetHorizontalFactor ();
      double Hfac2 = coil2.GetHorizontalFactor ();
      double Hfac3 = coil3.GetHorizontalFactor ();
      double Hfac4 = coil4.GetHorizontalFactor ();
      double Hfac5 = coil5.GetHorizontalFactor ();
      double Hfac6 = coil6.GetHorizontalFactor ();
      double Vfaco = ohmic.GetVerticalFactor   ();
      double Vfac1 = coil1.GetVerticalFactor   ();
      double Vfac2 = coil2.GetVerticalFactor   ();
      double Vfac3 = coil3.GetVerticalFactor   ();
      double Vfac4 = coil4.GetVerticalFactor   ();
      double Vfac5 = coil5.GetVerticalFactor   ();
      double Vfac6 = coil6.GetVerticalFactor   ();
      double f1a   = f1[Nr];
      double H1a   = HPfunc(1, Nr);
      double Bv    = sqrt(2.) * f1a * (log (8./epsa) - 1.5 - H1a);
      
      double Ic    = Bv /(Hfaco + Hfac1 + Hfac2 + Hfac3 + Hfac4 + Hfac5 + Hfac6);
      Ito          = Wo * Ic;
      It1          = W1 * Ic;
      It2          = W2 * Ic;
      It3          = W3 * Ic;
      It4          = W4 * Ic;
      It5          = W5 * Ic;
      It6          = W6 * Ic;
      
      double Fsin  = Vfaco + Vfac1 + Vfac2 + Vfac3 + Vfac4 + Vfac5 + Vfac6;
      
      printf ("\n");
      printf ("B_v = %11.4e  Fsin = %11.4e\n\n", - epsa * (f1a/2.) * (log (8./epsa) - 1.5 - H1a), Fsin);
      printf ("Ohmic  :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      Ro + Hshift, Zo - vshift, Hfaco, Vfaco, Ito);
      printf ("Coil 1 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R1 + Hshift, Z1 - vshift, Hfac1, Vfac1, It1);
      printf ("Coil 2 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R2 + Hshift, Z2 - vshift, Hfac2, Vfac2, It2);
      printf ("Coil 3 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R3 + Hshift, Z3 - vshift, Hfac3, Vfac3, It3);
      printf ("Coil 4 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R4 + Hshift, Z4 - vshift, Hfac4, Vfac4, It4);
      printf ("Coil 5 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R5 + Hshift, Z5 - vshift, Hfac5, Vfac5, It5);
      printf ("Coil 6 :  R = %11.4e  Z = %11.4e  Hfac = %11.4e  Vfac = %11.4e  It  = %11.4e\n",
	      R6 + Hshift, Z6 - vshift, Hfac6, Vfac6, It6);
      printf ("\n");
  
      // .........................
      // Rescale shaping functions
      // .........................
      for (int n = 2; n <= Ns; n++)
	{
	  double Hnfac = Ic * (    coil1.GetCosFactor (epsa, n) + coil2.GetCosFactor (epsa, n)
				 + coil3.GetCosFactor (epsa, n) + coil4.GetCosFactor (epsa, n)
				 + coil5.GetCosFactor (epsa, n) + coil6.GetCosFactor (epsa, n)
				 + ohmic.GetCosFactor (epsa, n));
	  double Vnfac = Ic * (    coil1.GetSinFactor (epsa, n) + coil2.GetSinFactor (epsa, n)
				 + coil3.GetSinFactor (epsa, n) + coil4.GetSinFactor (epsa, n)
				 + coil5.GetSinFactor (epsa, n) + coil6.GetSinFactor (epsa, n)
				 + ohmic.GetSinFactor (epsa, n));
	  
	  double Hnam  = (HPfunc(n, Nr) + double (n - 1) * HHfunc(n, Nr)) /double (2*n);
	  double Vnam  = (VPfunc(n, Nr) + double (n - 1) * VVfunc(n, Nr)) /double (2*n);
	  
	  double Hnfc  = Hnfac /Hnam /f1a;
	  double Vnfc  = Vnfac /Vnam /f1a;
	  
	  for (int i = 0; i <= Nr; i++)
	    {
	      HHfunc(n, i) = HHfunc(n, i) * Hnfc; 
	      HPfunc(n, i) = HPfunc(n, i) * Hnfc; 
	      VVfunc(n, i) = VVfunc(n, i) * Vnfc; 
	      VPfunc(n, i) = VPfunc(n, i) * Vnfc; 
	    }
	  
	  printf ("n = %3d:  Hfac = %11.4e  Vfac = %11.4e  Hna = %11.4e  Vna = %11.4e\n",
		  n, Hnfac, Vnfac, HHfunc(n, Nr), VVfunc(n, Nr));
	}
    }
  else
    {
      // ............................................
      // Read shaping data from file Inputs/Shape.txt
      // ............................................
      int     nshape;
      double  hnax, vnax;
      FILE*   file = OpenFiler ("Inputs/Shape.txt");
      double* hna  = new double[nshape+2];
      double* vna  = new double[nshape+2];

      fscanf (file, "%d", &nshape);
      for (int i = 0; i < nshape; i++)
	{
	  fscanf (file, "%lf %lf", &hnax, &vnax);
	  hna[i+2] = hnax;
	  vna[i+2] = vnax;
	}
      fclose (file);

      if (nshape > Ns)
	nshape = Ns;

      double f1a = f1[Nr];
      double H1a = HPfunc(1, Nr);
		 
      printf ("\n");
      printf ("B_v = %11.4e\n\n", - epsa * (f1a/2.) * (log (8./epsa) - 1.5 - H1a));
  
      // .........................
      // Rescale shaping functions
      // .........................
      for (int n = 2; n <= Ns; n++)
	{
	  double Hnam = (HPfunc(n, Nr) + double (n - 1) * HHfunc(n, Nr)) /double (2*n);
	  double Vnam = (VPfunc(n, Nr) + double (n - 1) * VVfunc(n, Nr)) /double (2*n);


	  double Hnfc, Vnfc;
	  if (n <= nshape+1)
	    {
	      Hnfc = hna[n];
	      Vnfc = vna[n];
	    }
	  else
	    {
	      Hnfc = 0.;
	      Vnfc = 0.;
	    }

	  double qnc =   f1a * Hnam * cos (double (n) * M_PI) * Hnfc /HPfunc(n, Nr);
	  double qns = - f1a * Vnam * cos (double (n) * M_PI) * Vnfc /VPfunc(n, Nr);
	  	  
	  for (int i = 0; i <= Nr; i++)
	    {
	   
	      HHfunc(n, i) = HHfunc(n, i) * Hnfc /HPfunc(n, Nr); 
	      HPfunc(n, i) = HPfunc(n, i) * Hnfc /HPfunc(n, Nr); 
	      VVfunc(n, i) = VVfunc(n, i) * Vnfc /VPfunc(n, Nr); 
	      VPfunc(n, i) = VPfunc(n, i) * Vnfc /VPfunc(n, Nr); 
	    }

	  printf ("n = %3d:  Hna = %11.4e  Vna = %11.4e  qnc = %11.4e  qns = %11.4e\n",
		  n, Hnfc, Vnfc, qnc, qns);
	}

      delete[] hna; delete[] vna;
    }
  
  // .............................
  // Interpolate shaping functions
  // .............................
  gsl_spline_init (g2spline, rr, g2,  Nr+1);

  double* data = new double[Nr+1];

  for (int i = 0; i <= Nr; i++)
    data[i] = HHfunc(1, i);
  gsl_spline_init (HHspline[1], rr, data,  Nr+1);

  for (int i = 0; i <= Nr; i++)
    data[i] = HPfunc(1, i);
  gsl_spline_init (HPspline[1], rr, data, Nr+1);

  for (int n = 2; n <= Ns; n++)
    {
      for (int i = 0; i <= Nr; i++)
	data[i] = HHfunc(n, i);
      gsl_spline_init (HHspline[n], rr, data,  Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = HPfunc(n, i);
      gsl_spline_init (HPspline[n], rr, data,  Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VVfunc(n, i);
      gsl_spline_init (VVspline[n], rr, data,  Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VPfunc(n, i);
      gsl_spline_init (VPspline[n], rr, data,  Nr+1);
    }

  delete[] data;
  
  // ................................
  // Calculate magnetic flux-surfaces
  // ................................
  for (int i = 1; i <= Nf; i++)
    {
     double rf = double (i) /double (Nf);

     double* hn = new double[Ns+1];
     for (int n = 1; n <= Ns; n++)
       hn[n] = gsl_spline_eval (HHspline[n], rf, HHacc[n]);

     double* vn = new double[Ns+1];
     for (int n = 2; n <= Ns; n++)
       vn[n] = gsl_spline_eval (VVspline[n], rf, VVacc[n]);

     double p = rf*rf*rf/8. - rf*hn[1]/2.;
     for (int n = 2; n <= Ns; n++)
       p += - double (n - 1) * hn[n]*hn[n] /2./rf - double (n - 1) * vn[n]*vn[n] /2./rf;
	 
     for (int j = 0; j <= Nw; j++)
       {
	 double w = double (j) * 2.*M_PI /double (Nw);
	 
	 double R = 1. - epsa*rf*cos(w) + epsa*epsa*epsa*p*cos(w) + epsa*epsa*hn[1];
	 for (int n = 2; n <= Ns; n++)
	   R += epsa*epsa * (hn[n] * cos (double (n - 1) * w) + vn[n] * sin (double (n - 1) * w));
	 
	 double Z = epsa*rf*sin(w) - epsa*epsa*epsa*p*sin(w);
	 for (int n = 2; n <= Ns; n++)
	   Z += epsa*epsa * (hn[n] * sin (double (n - 1) * w) - vn[n] * cos (double (n - 1) * w));
	 
	 RR(i-1, j) = R;
	 ZZ(i-1, j) = Z;
       }

     delete[] hn; delete[] vn;
    }

  // .....................
  // Integrate f3 equation
  // .....................
  double* y1   = new double[1];
  double* err1 = new double[1];
  rhs_chooser  = 1;

  f1c        = 1./qc;
  double H2c = HPfunc(2, 0);
  double V2c = VPfunc(2, 0);

  count = 0;
  h     = h0;
  r     = rr[0];
  f3[0] = 0.;
  q0[0] = qc;
  q2[0] = qc * (1. + eps*eps * (H2c*H2c + V2c*V2c));
  It[0] = 0.;
  Ip[0] = 0.;
  
  r     = eps;
  y1[0] = - f1c*f1c * (H2c*H2c + V2c*V2c) * r*r*r*r;
  
  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  RK4RK5Adaptive(1, r, y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      RK4RK5Fixed(1, r, y1, err1, rr[i] - r);
      
      f3[i] = y1[0] /f1[i];
      q0[i] = rr[i]*rr[i] /f1[i];
      q2[i] = rr[i]*rr[i] * (1. + epsa*epsa*(g2[i] - f3[i]/f1[i])) /f1[i];
 
      double gr2 = 3.*rr[i]*rr[i]/4. - HHfunc(1, i) + HPfunc(1, i)*HPfunc(1, i)/2.;
      for (int n = 2; n <= Ns; n++)
	{
	  gr2 += (HPfunc(n, i)*HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i)*HHfunc(n, i)/r/r)/2.;
	  gr2 += (VPfunc(n, i)*VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i)*VVfunc(n, i)/r/r)/2.;
	}

      It[i] =   2.*M_PI * (f1[i] * (1. + epsa*epsa*gr2) + epsa*epsa*f3[i]);
      Ip[i] = - 2.*M_PI * g2[i];
    }

  delete[] y1; delete[] err1;

  gsl_spline_init (Itspline, rr, It,  Nr+1);
  gsl_spline_init (Ipspline, rr, Ip,  Nr+1);

  for (int i = 0; i <= Nr; i++)
    {
      Jt[i] = gsl_spline_eval_deriv (Itspline, rr[i], Itacc);
      Jp[i] = gsl_spline_eval_deriv (Ipspline, rr[i], Ipacc);
    }

  printf ("\n");
  printf ("q2c = %11.4e  q0a  = %11.4e  q2a  = %11.4e  Ip = %11.4e  It = %11.4e\n",
	  q2[0], q0[Nr], q2[Nr], Ip[Nr], It[Nr]);
  
  // ..........................
  // Output data to netcdf file
  // ..........................
  double* Hna  = new double[Ns+1];
  double* Vna  = new double[Ns+1];
  double* npol = new double[Ns+1];

  double* Rdata = new double[Nf*(Nw+1)];
  double* Zdata = new double[Nf*(Nw+1)];
  
  double Rcoil[59], Zcoil[59];
  double Rc1[5], Zc1[5];
  coil1.GetCoordinates (Rc1, Zc1);
  double Rc2[5], Zc2[5];
  coil2.GetCoordinates (Rc2, Zc2);
  double Rc3[5], Zc3[5];
  coil3.GetCoordinates (Rc3, Zc3);
  double Rc4[5], Zc4[5];
  coil4.GetCoordinates (Rc4, Zc4);
  double Rc5[5], Zc5[5];
  coil5.GetCoordinates (Rc5, Zc5);
  double Rc6[5], Zc6[5];
  coil6.GetCoordinates (Rc6, Zc6);
  double Ro[29], Zo[29];
  ohmic.GetCoordinates (Ro, Zo);
  
  memcpy (Rcoil,      Rc1, 5 * sizeof (double));
  memcpy (Rcoil +  5, Rc2, 5 * sizeof (double));
  memcpy (Rcoil + 10, Rc3, 5 * sizeof (double));
  memcpy (Rcoil + 15, Rc4, 5 * sizeof (double));
  memcpy (Rcoil + 20, Rc5, 5 * sizeof (double));
  memcpy (Rcoil + 25, Rc6, 5 * sizeof (double));
  memcpy (Rcoil + 30, Ro, 29 * sizeof (double));
  
  memcpy (Zcoil,      Zc1, 5 * sizeof (double));
  memcpy (Zcoil +  5, Zc2, 5 * sizeof (double));
  memcpy (Zcoil + 10, Zc3, 5 * sizeof (double));
  memcpy (Zcoil + 15, Zc4, 5 * sizeof (double));
  memcpy (Zcoil + 20, Zc5, 5 * sizeof (double));
  memcpy (Zcoil + 25, Zc6, 5 * sizeof (double));
  memcpy (Zcoil + 30, Zo, 29 * sizeof (double));

  for (int n = 0; n <= Ns; n++)
    npol[n] = double (n);
  
  Hna[0] = 0.;
  Vna[0] = 0.;
  Hna[1] = HHfunc(1, Nr);
  Vna[1] = 0.;
  for (int n = 2; n <= Ns; n++)
    {
      Hna[n] = HHfunc(n, Nr);
      Vna[n] = VVfunc(n, Nr);
    }

  for (int n = 0; n < Nf; n++)
    for (int i = 0; i <= Nw; i++)
      {
	Rdata[i + n*(Nw+1)] = RR(n, i);
	Zdata[i + n*(Nw+1)] = ZZ(n, i);
      }
  
  int nerr = 0, dataFile;
  nerr = nc_create ("Plots/Equilibrium.nc", NC_CLOBBER, &dataFile);
  
  if (nerr != 0)
    {
      printf ("Error opening Plots/Equilibrium.nc\n");
      exit (1);
    }

  int r_d, s_d, f_d, w_d, c_d;
  nerr += nc_def_dim (dataFile, "Nr", Nr+1, &r_d);
  nerr += nc_def_dim (dataFile, "Ns", Ns+1, &s_d);
  nerr += nc_def_dim (dataFile, "Nf", Nf,   &f_d);
  nerr += nc_def_dim (dataFile, "Nw", Nw+1, &w_d);
  if (CFLAG)
    nerr += nc_def_dim (dataFile, "Nc", 59, &c_d);
  else
    nerr += nc_def_dim (dataFile, "Nc", 0,  &c_d);

  int shape_d[2];
  shape_d[0] = s_d;
  shape_d[1] = r_d;

  int flux_d[2];
  flux_d[0]  = f_d;
  flux_d[1]  = w_d;

  int r_x, g2_x, p2_x, Hn_x, Hnp_x, Vn_x, Vnp_x, R_x, Z_x;
  int f1_x, f3_x, q0_x, q2_x, It_x, Ip_x, Jt_x, Jp_x;
  int Rc_x, Zc_x, n_x, Hna_x, Vna_x;
  nerr += nc_def_var (dataFile, "r",   NC_DOUBLE, 1, &r_d,    &r_x);
  nerr += nc_def_var (dataFile, "g_2", NC_DOUBLE, 1, &r_d,    &g2_x);
  nerr += nc_def_var (dataFile, "p_2", NC_DOUBLE, 1, &r_d,    &p2_x);
  nerr += nc_def_var (dataFile, "Hn",  NC_DOUBLE, 2, shape_d, &Hn_x);
  nerr += nc_def_var (dataFile, "Hnp", NC_DOUBLE, 2, shape_d, &Hnp_x);
  nerr += nc_def_var (dataFile, "Vn",  NC_DOUBLE, 2, shape_d, &Vn_x);
  nerr += nc_def_var (dataFile, "Vnp", NC_DOUBLE, 2, shape_d, &Vnp_x);
  nerr += nc_def_var (dataFile, "R",   NC_DOUBLE, 2, flux_d,  &R_x);
  nerr += nc_def_var (dataFile, "Z",   NC_DOUBLE, 2, flux_d,  &Z_x);
  nerr += nc_def_var (dataFile, "f_1", NC_DOUBLE, 1, &r_d,    &f1_x);
  nerr += nc_def_var (dataFile, "f_3", NC_DOUBLE, 1, &r_d,    &f3_x);
  nerr += nc_def_var (dataFile, "q_0", NC_DOUBLE, 1, &r_d,    &q0_x);
  nerr += nc_def_var (dataFile, "q_2", NC_DOUBLE, 1, &r_d,    &q2_x);
  nerr += nc_def_var (dataFile, "I_t", NC_DOUBLE, 1, &r_d,    &It_x);
  nerr += nc_def_var (dataFile, "I_p", NC_DOUBLE, 1, &r_d,    &Ip_x);
  nerr += nc_def_var (dataFile, "J_t", NC_DOUBLE, 1, &r_d,    &Jt_x);
  nerr += nc_def_var (dataFile, "J_p", NC_DOUBLE, 1, &r_d,    &Jp_x);
  nerr += nc_def_var (dataFile, "R_c", NC_DOUBLE, 1, &c_d,    &Rc_x);
  nerr += nc_def_var (dataFile, "Z_c", NC_DOUBLE, 1, &c_d,    &Zc_x);
  nerr += nc_def_var (dataFile, "n",   NC_DOUBLE, 1, &s_d,    &n_x);
  nerr += nc_def_var (dataFile, "Hna", NC_DOUBLE, 1, &s_d,    &Hna_x);
  nerr += nc_def_var (dataFile, "Vna", NC_DOUBLE, 1, &s_d,    &Vna_x);

  nerr += nc_enddef (dataFile);

  if (nerr != 0)
    {
      printf ("Error defining variables in Output/Equilibrium.nc\n");
      exit (1);
    }

  nerr += nc_put_var_double (dataFile, r_x,   rr);
  nerr += nc_put_var_double (dataFile, g2_x,  g2);
  nerr += nc_put_var_double (dataFile, p2_x,  p2);
  nerr += nc_put_var_double (dataFile, Hn_x,  HHfunc.data());
  nerr += nc_put_var_double (dataFile, Hnp_x, HPfunc.data());
  nerr += nc_put_var_double (dataFile, Vn_x,  VVfunc.data());
  nerr += nc_put_var_double (dataFile, Vnp_x, VPfunc.data());
  nerr += nc_put_var_double (dataFile, R_x,   Rdata);
  nerr += nc_put_var_double (dataFile, Z_x,   Zdata);
  nerr += nc_put_var_double (dataFile, f1_x,  f1);
  nerr += nc_put_var_double (dataFile, f3_x,  f3);
  nerr += nc_put_var_double (dataFile, q0_x,  q0);
  nerr += nc_put_var_double (dataFile, q2_x,  q2);
  nerr += nc_put_var_double (dataFile, It_x,  It);
  nerr += nc_put_var_double (dataFile, Ip_x,  Ip);
  nerr += nc_put_var_double (dataFile, Jt_x,  Jt);
  nerr += nc_put_var_double (dataFile, Jp_x,  Jp);
  nerr += nc_put_var_double (dataFile, Rc_x,  Rcoil);
  nerr += nc_put_var_double (dataFile, Zc_x,  Zcoil);
  nerr += nc_put_var_double (dataFile, n_x,   npol);
  nerr += nc_put_var_double (dataFile, Hna_x, Hna);
  nerr += nc_put_var_double (dataFile, Vna_x, Vna);

  if (nerr != 0)
    {
      printf ("Error writing Output/Equilibrium.nc\n");
      exit (1);
    }
  
   nerr += nc_close (dataFile);

  if (nerr != 0)
    {
      printf ("Error closing Output/Equilibrium.nc\n");
      exit (1);
    }
  
  // ........
  // Clean up
  // ........
  delete[] rr;

  delete[] p2; delete[] f1; delete[] f3; delete[] g2; delete[] q0;  delete[] q2;
  delete[] It; delete[] Ip; delete[] Jt; delete[] Jp;

  delete[] pp2; delete[] ff1;

  gsl_spline_free (pp2spline);
  gsl_spline_free (ff1spline);

  gsl_interp_accel_free (pp2acc);
  gsl_interp_accel_free (ff1acc);

  gsl_spline_free (g2spline);
  gsl_spline_free (Itspline);
  gsl_spline_free (Ipspline);

  gsl_interp_accel_free (g2acc);
  gsl_interp_accel_free (Itacc);
  gsl_interp_accel_free (Ipacc);
  
  for (int i = 0; i <= Ns; i++)
    {
      gsl_spline_free (HHspline[i]);
      gsl_spline_free (VVspline[i]);
      gsl_spline_free (HPspline[i]);
      gsl_spline_free (VPspline[i]);

      gsl_interp_accel_free (HHacc[i]);
      gsl_interp_accel_free (VVacc[i]);
      gsl_interp_accel_free (HPacc[i]);
      gsl_interp_accel_free (VPacc[i]);
    }
  
  delete[] HHspline; delete[] VVspline; delete[] HPspline; delete[] VPspline;
  delete[] HHacc;    delete[] VVacc;    delete[] HPacc;    delete[] VPacc;

  delete[] Hna; delete[] Vna; delete[] npol;

  delete[] Rdata; delete[] Zdata;
}

// ########################
// Function to return f1(r)
// ########################
double Equilibrium::Getf1 (double r)
{
  if (flg)
    return gsl_spline_eval (ff1spline, r, ff1acc);
  else    
    return (1. - pow (1. - r*r, nu)) /nu/qc; 
}

// #########################
// Function to return f1'(r)
// #########################
double Equilibrium::Getf1p (double r)
{
  if (flg)
    return gsl_spline_eval_deriv (ff1spline, r, ff1acc);
  else
    return 2. * r * pow (1. - r*r, nu-1.) /qc;
}

// ########################
// Function to return p2(r)
// ########################
double Equilibrium::Getp2 (double r)
{
  if (flg)
    return gsl_spline_eval (pp2spline, r, pp2acc);
  else
    return pc * pow (1. - r*r, mu);
}

// #########################
// Function to return p2'(r)
// #########################
double Equilibrium::Getp2p (double r)
{
  if (flg)
    return gsl_spline_eval_deriv (pp2spline, r, pp2acc);
  else
    return - 2. * pc * mu * r * pow (1. - r*r, mu-1.);
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Equilibrium::Rhs (double r, double* y, double* dydr)
{
  if (rhs_chooser == 0)
    {
      // ...............................................
      // Right-hand sides for shaping function equations
      // ...............................................

      double* Hn   = new double[Ns+1];
      double* Hnp  = new double[Ns+1];
      double* Hnpp = new double[Ns+1];
      double* Vn   = new double[Ns+1];
      double* Vnp  = new double[Ns+1];
      double* Vnpp = new double[Ns+1];
      
      double g2 = y[0];

      Hn [1] = y[1];
      Hnp[1] = y[2];
      int j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  Hn [n] = y[j]; j++;
	  Hnp[n] = y[j]; j++;
	  Vn [n] = y[j]; j++;
	  Vnp[n] = y[j]; j++;
	}

      double f1  = Getf1 (r);
      double f1p = Getf1p(r);
      double p2p = Getp2p(r);
      
      double facf = 2.*f1p/f1 - 1./r;
      double facp = 2.*r*r*r*p2p/f1/f1;
      
      double g2p = - p2p - f1*f1p/r/r;

      Hnpp[1] = - facf * Hnp[1] - 1. + facp;

      for (int n = 2; n <= Ns; n++)
	{
	  Hnpp[n] = - facf * Hnp[n] + double (n*n - 1) * Hn[n]/r/r;
	  Vnpp[n] = - facf * Vnp[n] + double (n*n - 1) * Vn[n]/r/r;
  	}
      
      dydr[0] = g2p;
      dydr[1] = Hnp [1];
      dydr[2] = Hnpp[1];
      j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  dydr[j] = Hnp [n]; j++;
	  dydr[j] = Hnpp[n]; j++;
	  dydr[j] = Vnp [n]; j++;
	  dydr[j] = Vnpp[n]; j++;
	}

      delete[] Hn; delete[] Hnp; delete[] Hnpp; delete[] Vn; delete[] Vnp; delete[] Vnpp;
     }
  else if (rhs_chooser == 1)
    {
      // ................................
      // Right-hand side for f13 equation
      // ................................

      double* Hn  = new double[Ns+1];
      double* Hnp = new double[Ns+1];
      double* Vn  = new double[Ns+1];
      double* Vnp = new double[Ns+1];

      double f1  = Getf1 (r);
      double f1p = Getf1p(r);
      double p2p = Getp2p(r);
      
      double g2 = gsl_spline_eval (g2spline, r, g2acc);

      Hn [1] = gsl_spline_eval (HHspline[1], r, HHacc[1]);
      Hnp[1] = gsl_spline_eval (HPspline[1], r, HPacc[1]);
      for (int n = 2; n <= Ns; n++)
	{
	  Hn [n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
	  Hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
	  Vn [n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
	  Vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
	}

      double fac1 = 0., fac2 = 0.;
      for (int n = 2; n <= Ns; n++)
	{
	  fac1 += Hnp[n]*Hnp[n] + 2. * double (n*n - 1) * Hnp[n]*Hn[n]/r - double (n*n - 1) * Hn[n]*Hn[n]/r/r;
	  fac1 += Vnp[n]*Vnp[n] + 2. * double (n*n - 1) * Vnp[n]*Vn[n]/r - double (n*n - 1) * Vn[n]*Vn[n]/r/r;

	  fac2 += (3.*Hnp[n]*Hnp[n] - double (n*n - 1) * Hn[n]*Hn[n]/r/r) /2.;
	  fac2 += (3.*Vnp[n]*Vnp[n] - double (n*n - 1) * Vn[n]*Vn[n]/r/r) /2.;
	}
      
      double f13p =
	- f1*f1 * (3.*r*r/2. - 2.*r*Hnp[1] + Hnp[1]*Hnp[1] + fac1)/r
	+ f1*f1p * (g2 - 3.*r*r/4 + Hn[1] + 3.*Hnp[1]*Hnp[1]/2. + fac2)
	+ r*r*p2p * (g2 + r*r/2. - 3.*r*Hnp[1] - 2.*Hn[1]);

      dydr[0] = f13p;

      delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;
    }
  else if (rhs_chooser == 2)
    {
      // .......................................
      // Right-hand sides for pedestal equations
      // .......................................

      double f1  = Getf1(r);
      double p2p = Getp2p(r);
      
      dydr[0] = - (Pped/2./wped) * ((1. - r*r) /(1. - rped*rped)) /cosh ((r - rped) /wped) /cosh ((r - rped) /wped);
      dydr[1] = - 2. * 1.46 * sqrt (epsa) * pow (r, 2.5) * dydr[0];
    }
}

// ###################################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step length fourth-order/fifth-order Runge-Kutta scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... safety factor
//     T       ... step length cannot change by more than this factor from step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step length
//     h_max   ... maximum allowable step length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step while attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
// ####################################################################################
void Equilibrium::RK4RK5Adaptive (int neqns, double& x, double* y, double& h, 
				  double& t_err, double acc, double S, double T, int& rept,
				  int maxrept, double h_min, double h_max, int flag, 
				  int diag, FILE* file)
{
  double* y0  = new double[neqns];
  double* Err = new double[neqns];
  double  hin = h;

  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0[i] = y[i];

  // Take RK4/RK5 step 
  RK4RK5Fixed (neqns, x, y, Err, h);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err[i] / y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs (Err[i] / y[i]);
	  err2  = fabs (Err[i]);
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc / t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc / t_err), 0.25);

  // Prevent step length from changing by more than factor T
  if (h_est / h > T)
    h *= T;
  else if (h_est / h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Prevent step length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %11.4e hin = %11.4e err = %11.4e acc = %11.4e hout = %11.4e count = %3d\n", 
	     x, hin, t_err, acc, h, count);


  // Check if truncation error acceptable
  if ((t_err <= acc) || (count >= maxrept))
    {
      // If truncation error acceptable take step 
      rept  = count;
      count = 0;
    }
  else 
    {
      // If truncation error unacceptable repeat step 
      count++;
      x = x0;
      for (int i = 0; i < neqns; i++)
	y[i] = y0[i];
      RK4RK5Adaptive (neqns, x, y, h, t_err, acc, S, T, rept, 
		      maxrept, h_min, h_max, flag, diag, file);
    }

  delete[] y0; delete[] Err;
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step length fourth-order/fifth-order Runge-Kutta scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables 
//     err     ... array of errors
//     h       ... step length
//     
// #####################################################################
void Equilibrium::RK4RK5Fixed (int neqns, double& x, double* y, double* err, double h)
{
  double* dydx = new double[neqns];
  double* k1   = new double[neqns];
  double* k2   = new double[neqns];
  double* k3   = new double[neqns];
  double* k4   = new double[neqns];
  double* k5   = new double[neqns];
  double* k6   = new double[neqns];
  double* f    = new double[neqns];

  // First stage
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  Rhs (x + aa6 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k6[i] = h * dydx[i];
    }

  // Actual step 
  for (int i = 0; i < neqns; i++)
    {
      y  [i] = y[i] + cc1 * k1[i] + cc3 * k3[i] + cc4 * k4[i]               + cc6 * k6[i];
      err[i] =        ca1 * k1[i] + ca3 * k3[i] + ca4 * k4[i] + ca5 * k5[i] + ca6 * k6[i];
    }
  x += h;

  delete[] dydx; delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6; delete[] f;
}

// ################################
// Target function for zero finding
// ################################
double Equilibrium::Feval (double x)
{
  ohmic.SetParameters (Ro - Hshift, Zo - x, dRo, Wo);
  coil1.SetParameters (R1 - Hshift, Z1 - x, dR1, W1);
  coil2.SetParameters (R2 - Hshift, Z2 - x, dR2, W2);
  coil3.SetParameters (R3 - Hshift, Z3 - x, dR3, W3);
  coil4.SetParameters (R4 - Hshift, Z4 - x, dR4, W4);
  coil5.SetParameters (R5 - Hshift, Z5 - x, dR5, W5);
  coil6.SetParameters (R6 - Hshift, Z6 - x, dR6, W6);
 
  double Vfaco = ohmic.GetVerticalFactor ();
  double Vfac1 = coil1.GetVerticalFactor ();
  double Vfac2 = coil2.GetVerticalFactor ();
  double Vfac3 = coil3.GetVerticalFactor ();
  double Vfac4 = coil4.GetVerticalFactor ();
  double Vfac5 = coil5.GetVerticalFactor ();
  double Vfac6 = coil6.GetVerticalFactor ();

  return Vfaco + Vfac1 + Vfac2 + Vfac3 + Vfac4 + Vfac5 + Vfac6;
}

// ###################################################################
// Routine to find approximate root of F(x) = 0 using Ridder's method 
// Search takes place in interval (x1, x2)
// Interval is chopped into nint equal segments
// ###################################################################
double Equilibrium::RootFind (double x1, double x2)
{
  double F1, F2 = 0., root;

  // Chop search interval into nint segments  
  for (int i = 0; i < nint; i++)
    {
      double x1_seg = x1 + (x2 - x1) * double (i)     /double (nint);
      double x2_seg = x1 + (x2 - x1) * double (i + 1) /double (nint);
      
      if (i == 0) 
	F1 = Feval (x1_seg);
      else 
	F1 = F2;
      F2 = Feval (x2_seg);
      //printf ("%e %e %e %e\n", x1_seg, F1, x2_seg, F2);
      
      // Call Ridder's method for segment containing zero
      if (F1 * F2 < 0.)
	{
	  Ridder (x1_seg, x2_seg, F1, F2, root);
	  break;
	}
    }

  return root;
}

// ############################################
// Ridder's method for finding root of F(x) = 0
// ############################################
void Equilibrium::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) /2.;    
      double F3 = Feval (x3);
      
      // Iterate x using Ridder's method 
      xold = x;           
      x = x3 - (x3 - x1) * (F2 - F1) * F3 /
	(sqrt (F3 * F3 - F1 * F2) * fabs (F2 - F1));
      Fx = Feval (x);
       
      // Make new value of x upper/lower bound of refined search interval, as appropriate 
      if (Fx * F1 < 0.) 
	{  
	  x2 = x;           
	  F2 = Fx; 
	}
      else 
	{
	  x1 = x;
	  F1 = Fx; 
	}
      //printf ("%d %e %e\n", iter, x, Fx);
      iter++;
    } 
  // Iterate until absolute change in x falls below Eta
  while (fabs (x - xold) > Eta && fabs(Fx) > Eta && iter < Maxiter); 
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Equilibrium::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* Equilibrium::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

