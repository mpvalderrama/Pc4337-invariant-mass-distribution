#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

/* The masses of the decaying particle and its decay products:
 */
const double M_proton = 938.3;
const double M_JPsi = 3096.9;

const double M_Bs0 = 5366.9;

/* Useful physical and mathematical constants:
 */
const double hbc = 197.33;

const double Pi = 3.141592;

/* A small epsilon, mostly to guarantee that one does not accidentally
   go outside of the kinematical limits of the decays:
 */
const double Epsilon = 1.0e-4;

/* Number of plot points for the continuous distribution shapes:
 */
#define N_PLOT 400

/* Integration points within each bin: it is different than in the fit
   programs, but then this serves as a crosscheck (if the results are
   a bit off, it will be clear that the original N_INT_BIN
   in the fit programs were a bit off).
 */
#define N_INT_BIN 200


/* Here we choose what distribution (12, 12+z) and parametes
   (all or only the combined production/decay couplings) to fit:
 */

#define COMPUTE_T4_12 0 /* Fits to m12, all parameters */
#define COMPUTE_T4_12_PD 1 /* Fits to m12,
			      only production/decay parameters */
#define COMPUTE_T4_12_z_PD 0 /* Fits to m12 and z distributions,
				only production/decay parameters */

/* The m12 invariant mass distribution for meson-baryon rescattering:
 */
double
imd_T4_12 (double m12, double beta, double gz, double phiz,
	   double g0, double phi0, double g1, double phi1, double h1,
	   double d8, double dE, double c8, double c10,
	   double Lambda);

double
imd_T4_12_extra (double m12, double beta,
		 double gz, double phiz,
		 double g0, double phi0,
		 double g1, double phi1, double h1,
		 double d8, double dE, double c8, double c10,
		 double Lambda,
		 double *I_peak,
		 double *I_reflection,
		 double *I_nonres);

/* Proton helicity distribution:
 */
double
imd_T4_z (double z, double beta, double gz, double phiz,
	  double g0, double phi0, double g1, double phi1, double h1,
	  double d8, double dE, double c8, double c10,
	  double Lambda);

/* Programs for fitting the different distributions: currently
   the only ones implemented are 12 (J/Psi proton) and 12+z
   (J/Psi proton and the proton helicity distribution)/
 */
double
imd_fit_T4_12_poisson (double *betaT,
		       double *gz_R, double *phiz_R,
		       double *g0_R, double *phi0_R,
		       double *g1_R, double *phi1_R,
		       double *h1_R,
		       double *d8, double *dE, double *c8, double *c10,
		       double betaT_step,
		       double gz_step, double phiz_step,
		       double g0_step, double phi0_step,
		       double g1_step, double phi1_step,
		       double h1_step,
		       double d8_step, double dE_step,
		       double c8_step, double c10_step,
		       double m12_data[], double events_data[],
		       double background_data[],
		       int n_data, double Lambda);


/* This program fits only the production/decay couplings to the m12
   distribution  with a given meson-baryon interaction profile.
 */
double
imd_fit_T4_12_pd_poisson (double *betaT,
			  double *gz_R, double *phiz_R,
			  double *g0_R, double *phi0_R,
			  double *g1_R, double *phi1_R,
			  double *h1_R,
			  double betaT_step,
			  double gz_step, double phiz_step,
			  double g0_step, double phi0_step,
			  double g1_step, double phi1_step,
			  double h1_step,
			  double m12_data[], double events_data[],
			  double background_data[],
			  int n_data,
			  double d8, double dE, double c8, double c10,
			  double Lambda);

/* This program fits only the production/decay couplings to the m12
   and z distributions with a given meson-baryon interaction profile.
 */
double
imd_fit_T4_12_z_pd_poisson (double *betaT,
			    double *gz_R, double *phiz_R,
			    double *g0_R, double *phi0_R,
			    double *g1_R, double *phi1_R,
			    double *h1_R,
			    double betaT_step,
			    double gz_step, double phiz_step,
			    double g0_step, double phi0_step,
			    double g1_step, double phi1_step,
			    double h1_step,
			    double m12_data[], double events_data[],
			    double background_data[],
			    double events_z_data[],
			    double background_z_data[],
			    int n_data,
			    double d8, double dE, double c8, double c10,
			    double Lambda);

/* Here is the main program: the particular type of fit is chosen
   by the #define COMPUTE_T4_{something} statement
   at the beginning of this file. Choose there.
 */
int
main (void)
{
  double m12;

  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;
  
  double h = (m12_max - m12_min)/((double)(N_PLOT - 1));
  int i;

  double z, z_min = -1.0, z_max = 1.0;
  double h_z = (z_max - z_min)/((double)(N_PLOT - 1));
  
  double imd;

  /* Initial parameters for the T-matrix fit
   */

  /* Meson-baryon couplings: */
  double d8 = -0.3 /hbc/hbc;
  double dE = +0.1 /hbc/hbc;
  double c8 = -0.6 /hbc/hbc ;
  double c10 = +0.9 /hbc/hbc ;

  double d8_step = 0.02 /hbc/hbc;
  double dE_step = 0.01 /hbc/hbc;
  double c8_step = 0.02 /hbc/hbc;
  double c10_step = 0.02 /hbc/hbc;

  /* Production+decay combined couplings: */
  double beta_T = 0.32e-2;

  double gz_R = 1.0e-8;
  double phiz_R = 0.4;
  
  double g0_R = 0.88e-8;
  double phi0_R = -0.6;  

  double g1_R = 0.65e-8;
  double phi1_R = -0.5;  

  double h1_R = 0.54;

  double beta_T_step = 0.1e-2;

  double gz_step = 3.0e-9;
  double phiz_step = 0.5;
  
  double g0_step = 3.0e-9;
  double phi0_step = 0.5;

  double g1_step = 3.0e-9;
  double phi1_step = 0.5;

  double h1_step = 0.3;
  
  double Lambda = 1000.0;
  
  double chisq;

  int ndata = 30;

  double events_comb_data[30], err_comb_data[30],
    background_comb_data[30];

  /* Here I list the data: definitely not the most practical way
     to include them, but it works!
   */
  double m12_data[30] = {4.03175,
			 4.04525,
			 4.05875,
			 4.07225,
			 4.08575,
			 4.09925,
			 4.11275,
			 4.12625,
			 4.13975,
			 4.15325,
			 4.16675,
			 4.18025,
			 4.19375,
			 4.20725,
			 4.22075,
			 4.23425,
			 4.24775,
			 4.26125,
			 4.27475,
			 4.28825,
			 4.30175,
			 4.31525,
			 4.32875,
			 4.34225,
			 4.35575,
			 4.36925,
			 4.38275,
			 4.39625,
			 4.40975,
			 4.42325};

  double events_data[30] = {2,
		     8,
		     9,
		     21,
		     21,
		     25,
		     35,
		     40,
		     26,
		     37,
		     44,
		     42,
		     34,
		     37,
		     41,
		     36,
		     36,
		     38,
		     37,
		     33,
		     36,
		     40,
		     51,
		     54,
		     31,
		     34,
		     21,
		     17,
		     14,
		     6};

  double err_data[30] = {1.414213562373,
		      2.828427124746,
		      3.0,
		      4.582575694956,
		      4.582575694956,
		      5.0,
		      5.9160797831,
		      6.324555320337,
		      5.099019513593,
		      6.082762530298,
		      6.633249580711,
		      6.480740698408,
		      5.830951894845,
		      6.082762530298,
		      6.403124237433,
		      6.0,
		      6.0,
		      6.164414002969,
		      6.082762530298,
		      5.744562646538,
		      6.0,
		      6.324555320337,
		      7.141428428543,
		      7.34846922835,
		      5.56776436283,
		      5.830951894845,
		      4.582575694956,
		      4.123105625618,
		      3.74165738677,
		      2.44948974278};

  double background_data[30] = {0.1816923,
				1.834083,
				2.674426,
				3.285455,
				3.746343,
				3.994186,
				4.386583,
				4.722361,
				5.006309,
				5.309543,
				5.194319,
				5.379957,
				5.522682,
				5.605806,
				5.827727,
				5.616472,
				5.733449,
				5.625611,
				5.634692,
				5.646524,
				5.462607,
				5.418713,
				5.236973,
				5.208698,
				4.935687,
				4.69976,
				4.333629,
				3.797254,
				3.096513,
				1.875946};
  
  double events_bar_data[30] = {2,
				8,
				17,
				20,
				24,
				33,
				18,
				35,
				24,
				33,
				43,
				42,
				37,
				40,
				43,
				33,
				39,
				33,
				40,
				39,
				27,
				39,
				45,
				43,
				39,
				43,
				26,
				26,
				10,
				6};

  double background_bar_data[30] = {0.2086673,
				    1.924011,
				    2.859444,
				    3.487688,
				    3.965139,
				    4.377289,
				    4.687739,
				    4.906383,
				    5.418845,
				    5.519567,
				    5.66469,
				    5.507661,
				    5.667747,
				    5.708237,
				    5.765443,
				    5.631274,
				    5.666646,
				    5.565474,
				    5.547987,
				    5.428842,
				    5.314474,
				    5.178233,
				    5.012619,
				    4.820317,
				    4.56275,
				    4.350043,
				    4.077329,
				    3.57677,
				    2.916644,
				    1.825048};
  
  double err_bar_data[30] = {1.41421356237,
			     2.828427124746,
			     4.123105625618,
			     4.472135955,
			     4.898979485566,
			     5.744562646538,
			     4.242640687119,
			     5.9160797831,
			     4.898979485566,
			     5.744562646538,
			     6.557438524302,
			     6.480740698408,
			     6.082762530298,
			     6.324555320337,
			     6.557438524302,
			     5.744562646538,
			     6.244997998398,
			     5.744562646538,
			     6.324555320337,
			     6.244997998398,
			     5.196152422707,
			     6.244997998398,
			     6.708203932499,
			     6.557438524302,
			     6.244997998398,
			     6.557438524302,
			     5.099019513593,
			     5.099019513593,
			     3.162277660168,
			     2.449489742783};

  double z_data[30] = {-0.966666666666667,
		       -0.9,
		       -0.833333333333333,
		       -0.766666666666667,
		       -0.7,
		       -0.633333333333333,
		       -0.566666666666667,
		       -0.5,
		       -0.433333333333333,
		       -0.366666666666667,
		       -0.3,
		       -0.233333333333333,
		       -0.166666666666667,
		       -0.1,
		       -0.0333333333333333,
		       0.0333333333333333,
		       0.1,
		       0.166666666666667,
		       0.233333333333333,
		       0.3,
		       0.366666666666667,
		       0.433333333333333,
		       0.5,
		       0.566666666666667,
		       0.633333333333333,
		       0.7,
		       0.766666666666667,
		       0.833333333333333,
		       0.9,
		       0.966666666666667};

  double events_z_data[30] = {24,
			      21,
			      35,
			      25,
			      34,
			      48,
			      30,
			      36,
			      31,
			      36,
			      23,
			      31,
			      32,
			      29,
			      23,
			      23,
			      32,
			      30,
			      36,
			      25,
			      25,
			      35,
			      37,
			      40,
			      30,
			      35,
			      26,
			      26,
			      28,
			      21};

  double background_z_data[30] = {5.251675,
				  4.891967,
				  4.802771,
				  4.649889,
				  4.519848,
				  4.359397,
				  4.357543,
				  4.331665,
				  4.231905,
				  4.141592,
				  4.052932,
				  4.068803,
				  4.116404,
				  4.057827,
				  4.092458,
				  4.101066,
				  4.169519,
				  4.172583,
				  4.159311,
				  4.284362,
				  4.338489,
				  4.493214,
				  4.481134,
				  4.592039,
				  4.67464,
				  4.885715,
				  4.961477,
				  5.091943,
				  5.275876,
				  5.534956};
  
  /* Change the m12 to MeV:
   */
  for (i=0; i < ndata; i++)
    {
      m12_data[i] = 1000.0 * m12_data[i];
    }

  /* Combination of the data:
   */
  for (i=0; i < ndata; i++)
    {
      events_comb_data[i] = (events_data[i] + events_bar_data[i])/2.0;
      err_comb_data[i] = sqrt(events_comb_data[i]);
      background_comb_data[i] =
	(background_data[i] + background_bar_data[i])/2.0;
    }

  /* First fit block: fits the 4-channel T-matrix to the m12 distro:
   */
  if (COMPUTE_T4_12)
    {
      /* Here there are some test values from a previous fit: they may
	 be used as the starting point for the next fit, or
	 alternatively one may comment out the fit program
	 and check the results with these parameters.
       */
      beta_T = 0.072995e-2;

      gz_R = 1.415179e-8;
      phiz_R = 2.574588;

      g0_R = 0.636120e-8;
      phi0_R = 4.584395;

      g1_R = 1.341118e-8;
      phi1_R = 2.420769;

      h1_R = 3.059468;

      d8 = -0.087767 /hbc/hbc;
      dE = -0.129388 /hbc/hbc;
      c8 = -0.202079 /hbc/hbc;
      c10 = -0.206854 /hbc/hbc;
      
      chisq = 17.045335;

      /* Now we run the fit:
       */
      chisq = imd_fit_T4_12_poisson (&beta_T,
				     &gz_R, &phiz_R,
				     &g0_R, &phi0_R,
				     &g1_R, &phi1_R, &h1_R,
				     &d8, &dE, &c8, &c10,
				     beta_T_step,
				     gz_step, phiz_step,
				     g0_step, phi0_step,
				     g1_step, phi1_step,
				     h1_step,
				     d8_step, dE_step, c8_step, c10_step,
				     m12_data, events_data,
				     background_data,
				     ndata, Lambda);
      
      /* Finally, we plot our results:

	 First, the J/Psi proton mass distribution (M12):
       */
      FILE *out = fopen ("M12-fit12-T4.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h*i;

	  imd = imd_T4_12 (m12, beta_T, gz_R, phiz_R,
			   g0_R, phi0_R, g1_R, phi1_R, h1_R,
			   d8, dE, c8, c10, Lambda);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);
            
      /* Second, the proton helicity mass distribution (Mz):
       */
      out = fopen ("Mz-fit12-Tz.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  z = z_min + h_z*i;

	  imd = imd_T4_z (z, beta_T, gz_R, phiz_R,
			  g0_R, phi0_R, g1_R, phi1_R, h1_R,
			  d8, dE, c8, c10, Lambda);

	  /* I have to normalize for the differences in bin sizes:
	   */
	  imd = imd / 13.5 / 15.0;
	  
	  printf ("%f , \t%f\n", z, imd);
	  fprintf (out, "%f , \t%f\n", z, imd);
	}

      fclose (out);

      /* Third, the J/Psi proton mass distribution as a histogram:
       */
      out = fopen ("M12-fit12-T4-histo.plot", "w");

      fprintf (out,
	       "# beta = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
            
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;
	  
	  double m12_min = M_JPsi + M_proton;
	  double m12_max = M_Bs0 - M_proton;
	  
	  double m12_a, m12_b, bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;

	  double I_bin_nonres = 0.0;
	  double inte_bin_nonres_pre, inte_bin_nonres;
	  
	  double events_th, nth;
	  double events_reflect, events_peak, events_nonres;
	  
	  m12 = m12_data[i];

	  m12_a = m12 - bin/2.0;
	  m12_b = m12 + bin/2.0;

	  if (m12_a < m12_min)
	    m12_a = m12_min + Epsilon;
	  if (m12_b > m12_max)
	    m12_b = m12_max - Epsilon;

	  h = (m12_b - m12_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  m12 = m12_a;
	  inte_bin_pre = imd_T4_12_extra (m12,
					  beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R,
					  g1_R, phi1_R, h1_R,
					  d8, dE, c8, c10, Lambda,
					  &inte_bin_peak_pre,
					  &inte_bin_reflect_pre,
					  &inte_bin_nonres_pre);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      m12 += h;

	      inte_bin = imd_T4_12_extra (m12,
					  beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R,
					  g1_R, phi1_R, h1_R,
					  d8, dE, c8, c10, Lambda,
					  &inte_bin_peak_pre,
					  &inte_bin_reflect_pre,
					  &inte_bin_nonres_pre);
     
	      I_bin += (inte_bin + inte_bin_pre) * h / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      I_bin_peak +=
		(inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;
	      I_bin_nonres +=
		(inte_bin_nonres + inte_bin_nonres_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      inte_bin_nonres_pre = inte_bin_nonres;
	    }

	  events_th = I_bin / bin;

	  events_peak = I_bin_peak / bin;
	  events_reflect = I_bin_reflect / bin;
	  events_nonres = I_bin_nonres / bin;
      
	  nth = events_th + background_data[i];

	  printf ("%f , \t%f\n", m12_data[i], nth);
	  
	  fprintf (out, "%f , \t%f, \t%f, \t%f, \t%f, \t%f\n",
		   m12_data[i]/1000.0, nth,
		   events_peak, events_reflect, events_nonres,
		   background_data[i]/events_nonres);
	}

      fclose (out);

      /* Fourth, the histogram for the proton helicity distribution
       */
   
      out = fopen ("Mz-fit12-T4-histo.plot", "w");

      fprintf (out,
	       "# beta = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;

	  double z_a, z_b, z_i;
	  double z_bin = 2.0 / ((double) (ndata)), bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;
	  
	  double events_z_th, nth_z;
	  double events_z_reflect, events_z_peak;

	  z_min = -1.0;
	  z_max = +1.0;

	  /* Here the bins are evenly spaced between z_min and z_max, and
	     it might be better to calculate them directly like this:
	   */
	  z_i = z_min + ((double) (2*i+1)) / ((double) (ndata));

	  z_a = z_i - z_bin/2.0;
	  z_b = z_i + z_bin/2.0;

	  /* Check that we do not under/overflow z:
	   */
	  if (z_a < z_min)
	    z_a = z_min;
	  if (z_b > z_max)
	    z_b = z_max;

	  h_z = (z_b - z_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  z = z_a;
	  inte_bin_pre = imd_T4_z (z,
				   beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R, h1_R,
				   d8, dE, c8, c10, Lambda);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      z += h_z;

	      inte_bin = imd_T4_z (z, beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R,
				   h1_R, d8, dE, c8, c10, Lambda);
     
	      I_bin += (inte_bin + inte_bin_pre) * h_z / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      /* Not implemented yet:
	      I_bin_peak +=
	      (inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      */
	    }

	  /* The events are the average over the bin, but then we have to
	     normalize them by the factor z_bib/bin so they give
	     the correct counts:
	   */
	  events_z_th = (I_bin / z_bin) * z_bin / bin;

	  /*
	  events_z_peak = I_bin_peak / bin;
	  events_z_reflect = I_bin_reflect / bin;
	  */
	  
	  nth_z = events_z_th + background_z_data[i];
      
	  printf ("%f , \t%f\n", z_i, nth_z);
	  fprintf (out, "%f , \t%f\n", z_i, nth_z);
	}

      fclose (out);

      /* Finally, we print the parameters of the fit:

	 I use a format that can be directly copied back into
	 the program for later checks.
       */

      printf ("\n#############\n\n");
	    
      printf ("beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n\n",
	      1e2 * beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      printf ("d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      printf ("chisq = %f\n;", chisq);
      

    }  

  /* Fit to m12, but only involving the production/decay combined
     couplings (the meson-baryon couplings are fixed):
   */
  if (COMPUTE_T4_12_PD)
    {
      /* Here I begin with a few scenarios, with their starting points
	 for the simplex minimization.

	 Put the one you want in the last position. 
       */
      
      /* Meson-baryon couplings (a):

	 (uses the default starting point)
       */
      d8 = -0.343 /hbc/hbc;
      dE = +0.0 /hbc/hbc;
      c8 = -0.0 /hbc/hbc;
      c10 = -0.481 /hbc/hbc;

      beta_T = 0.32e-2;
      gz_R = 1.0e-8;
      phiz_R = 0.4;
      g0_R = 0.88e-8;
      phi0_R = -0.6;  
      g1_R = 0.65e-8;
      phi1_R = -0.5;  
      h1_R = 0.54;
      
      /* Meson-baryon couplings (b):

	 (here I also include initial points for production/decay)
       */
      d8 = -0.372 /hbc/hbc;
      dE = +0.148 /hbc/hbc;
      c8 = -0.802 /hbc/hbc;
      c10 = -0.457 /hbc/hbc;

      beta_T = 0.151803e-2;
      gz_R = 0.874827e-8;
      phiz_R = 1.476878;
      g0_R = 3.9144643e-8;
      phi0_R = 1.280329;
      g1_R = 1.389569e-8;
      phi1_R = -2.884242;
      h1_R = 1.103743;


      /* Meson-baryon couplings (c):

	 (here I also include initial points for production/decay)
       */
      d8 = +0.139 /hbc/hbc;
      dE = +0.073 /hbc/hbc;
      c8 = -0.165 /hbc/hbc;
      c10 = -0.900 /hbc/hbc;

      beta_T = -0.617909e-2;
      gz_R = 1.593361e-8;
      phiz_R = 2.054721;
      g0_R = 5.268872e-8;
      phi0_R = 2.762407;
      g1_R = 50.213544e-8;
      phi1_R = 2.403764;
      h1_R = 0.182316;

      /* Meson-baryon couplings (d):

	 (here I also include initial points for production/decay)
       */
      d8 = -0.088 /hbc/hbc;
      dE = -0.129 /hbc/hbc;
      c8 = -0.202 /hbc/hbc;
      c10 = -0.207 /hbc/hbc;

      beta_T = 0.072995e-2;
      gz_R = 1.415179e-8;
      phiz_R = 2.574588;
      g0_R = 0.636120e-8;
      phi0_R = 4.584395;
      g1_R = 1.341118e-8;
      phi1_R = 2.420769;
      h1_R = 3.059468;

      
      chisq = imd_fit_T4_12_pd_poisson (&beta_T,
					&gz_R, &phiz_R,
					&g0_R, &phi0_R,
					&g1_R, &phi1_R, &h1_R,
					beta_T_step,
					gz_step, phiz_step,
					g0_step, phi0_step,
					g1_step, phi1_step,
					h1_step,
					m12_data, events_data,
					background_data,
					ndata, d8, dE, c8, c10,
					Lambda);
      
      /* Finally, we plot our results:

	 First, the J/Psi p mass distribution:
       */
      FILE *out = fopen ("M12-fit12-T4-PD.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h*i;

	  imd = imd_T4_12 (m12, beta_T, gz_R, phiz_R,
			   g0_R, phi0_R, g1_R, phi1_R, h1_R,
			   d8, dE, c8, c10, Lambda);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

            
      /* Second, the proton helicity mass distribution (Mz):
       */
      out = fopen ("Mz-fit12-T4-PD.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  z = z_min + h_z*i;

	  imd = imd_T4_z (z, beta_T, gz_R, phiz_R,
			  g0_R, phi0_R, g1_R, phi1_R, h1_R,
			  d8, dE, c8, c10, Lambda);

	  /* I have to normalize for the differences in bin sizes:
	   */
	  imd = imd / 13.5 / 15.0;
	  
	  printf ("%f , \t%f\n", z, imd);
	  fprintf (out, "%f , \t%f\n", z, imd);
	}

      fclose (out);

      /* Third, the histogram for the J/Psi p mass distribution
       */
      out = fopen ("M12-fit12-T4-PD-histo.plot", "w");

      fprintf (out,
	       "# beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;
	  
	  double m12_min = M_JPsi + M_proton;
	  double m12_max = M_Bs0 - M_proton;
	  
	  double m12_a, m12_b, bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;

	  double I_bin_nonres = 0.0;
	  double inte_bin_nonres_pre, inte_bin_nonres;
	  
	  double events_th, nth;
	  double events_reflect, events_peak, events_nonres;
	  
	  m12 = m12_data[i];

	  m12_a = m12 - bin/2.0;
	  m12_b = m12 + bin/2.0;

	  if (m12_a < m12_min)
	    m12_a = m12_min + Epsilon;
	  if (m12_b > m12_max)
	    m12_b = m12_max - Epsilon;

	  h = (m12_b - m12_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  m12 = m12_a;
	  inte_bin_pre = imd_T4_12_extra (m12,
					  beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R,
					  g1_R, phi1_R, h1_R,
					  d8, dE, c8, c10, Lambda,
					  &inte_bin_peak_pre,
					  &inte_bin_reflect_pre,
					  &inte_bin_nonres_pre);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      m12 += h;

	      inte_bin = imd_T4_12_extra (m12, beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R, g1_R, phi1_R,
					  h1_R, d8, dE, c8, c10, Lambda,
					  &inte_bin_peak,
					  &inte_bin_reflect,
					  &inte_bin_nonres);
     
	      I_bin += (inte_bin + inte_bin_pre) * h / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      I_bin_peak +=
		(inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;
	      I_bin_nonres +=
		(inte_bin_nonres + inte_bin_nonres_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      inte_bin_nonres_pre = inte_bin_nonres;
	    }

	  events_th = I_bin / bin;

	  events_peak = I_bin_peak / bin;
	  events_reflect = I_bin_reflect / bin;
	  events_nonres = I_bin_nonres / bin;
      
	  nth = events_th + background_data[i];
      
	  printf ("%f , \t%f\n", m12_data[i], nth);
	  
	  fprintf (out, "%f , \t%f, \t%f, \t%f, \t%f, \t%f\n",
		   m12_data[i]/1000.0, nth,
		   events_peak, events_reflect, events_nonres,
		   background_data[i]/events_nonres);
	}

      fclose (out);
      
      /* Fourth, the histogram for the proton helicity distribution
       */
   
      out = fopen ("Mz-fit12-T4-PD-histo.plot", "w");

      fprintf (out,
	       "# beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;

	  double z_a, z_b, z_i;
	  double z_bin = 2.0 / ((double) (ndata)), bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;
	  
	  double events_z_th, nth_z;
	  double events_z_reflect, events_z_peak;

	  z_min = -1.0;
	  z_max = +1.0;

	  /* Here the bins are evenly spaced between z_min and z_max, and
	     it might be better to calculate them directly like this:
	   */
	  z_i = z_min + ((double) (2*i+1)) / ((double) (ndata));

	  z_a = z_i - z_bin/2.0;
	  z_b = z_i + z_bin/2.0;

	  /* Check that we do not under/overflow z:
	   */
	  if (z_a < z_min)
	    z_a = z_min;
	  if (z_b > z_max)
	    z_b = z_max;

	  h_z = (z_b - z_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  z = z_a;
	  inte_bin_pre = imd_T4_z (z,
				   beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R, h1_R,
				   d8, dE, c8, c10, Lambda);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      z += h_z;

	      inte_bin = imd_T4_z (z, beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R,
				   h1_R, d8, dE, c8, c10, Lambda);
     
	      I_bin += (inte_bin + inte_bin_pre) * h_z / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      /*
	      I_bin_peak += (inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      */
	    }

	  /* The events are the average over the bin, but then we have to
	     normalize them by the factor z_bib/bin so they give
	     the correct counts:
	   */
	  events_z_th = (I_bin / z_bin) * z_bin / bin;

	  /*
	  events_z_peak = I_bin_peak / bin;
	  events_z_reflect = I_bin_reflect / bin;
	  */
	  
	  nth_z = events_z_th + background_z_data[i];
      
	  printf ("%f , \t%f\n", z_i, nth_z);
	  fprintf (out, "%f , \t%f\n", z_i, nth_z);
	}

      fclose (out);

      /* Finally, we print the production/coupling parameters
	 as well as the contact-range couplings:

	 
	 I use a format that can be directly copied back into
	 the program for later checks.
       */
      printf ("\n#############\n\n");
	    
      printf ("beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n\n",
	      1e2 * beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      printf ("d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      printf ("chisq = %f\n;", chisq);      
    }  


  
  if (COMPUTE_T4_12_z_PD)
    {
      /* Meson-baryon couplings:
       */
      d8 = -0.372 /hbc/hbc;
      dE = +0.148 /hbc/hbc;
      c8 = -0.802 /hbc/hbc;
      c10 = - 0.457 /hbc/hbc;
      
      /* Test initial values: change them according to your desired
	 starting point.
       */
      beta_T = 0.151803e-2;

      gz_R = 0.874827e-8;
      phiz_R = 1.476878;
      
      g0_R = 3.914464e-8;
      phi0_R = 1.180329;

      g1_R = 1.389569e-8;
      phi1_R = -2.884242;

      h1_R = 1.103743;
      
      chisq = imd_fit_T4_12_z_pd_poisson (&beta_T,
					  &gz_R, &phiz_R,
					  &g0_R, &phi0_R,
					  &g1_R, &phi1_R, &h1_R,
					  beta_T_step,
					  gz_step, phiz_step,
					  g0_step, phi0_step,
					  g1_step, phi1_step,
					  h1_step,
					  m12_data, events_data,
					  background_data,
					  events_z_data,
					  background_z_data,
					  ndata, d8, dE, c8, c10,
					  Lambda);
            
      /* Finally, we plot our results:

	 First, the J/Psi p mass distribution:
       */
      FILE *out = fopen ("M12-fit12z-T4-PD.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h*i;

	  imd = imd_T4_12 (m12, beta_T, gz_R, phiz_R,
			   g0_R, phi0_R, g1_R, phi1_R, h1_R,
			   d8, dE, c8, c10, Lambda);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

            
      /* Second, the proton helicity mass distribution (Mz):
       */
      out = fopen ("Mz-fit12z-T4-PD.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  z = z_min + h_z*i;

	  imd = imd_T4_z (z, beta_T, gz_R, phiz_R,
			  g0_R, phi0_R, g1_R, phi1_R, h1_R,
			  d8, dE, c8, c10, Lambda);

	  /* I have to normalize for the differences in bin sizes:
	   */
	  imd = imd / 13.5 / 15.0;
	  
	  printf ("%f , \t%f\n", z, imd);
	  fprintf (out, "%f , \t%f\n", z, imd);
	}

      fclose (out);

      /* Third, the histogram for the J/Psi p mass distribution
       */
   
      out = fopen ("M12-fit12z-T4-PD-histo.plot", "w");
      
      fprintf (out,
	       "# beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;
	  
	  double m12_min = M_JPsi + M_proton;
	  double m12_max = M_Bs0 - M_proton;
	  
	  double m12_a, m12_b, bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;

	  double I_bin_nonres = 0.0;
	  double inte_bin_nonres_pre, inte_bin_nonres;
	  
	  double events_th, nth;
	  double events_reflect, events_peak, events_nonres;
	  
	  m12 = m12_data[i];

	  m12_a = m12 - bin/2.0;
	  m12_b = m12 + bin/2.0;

	  if (m12_a < m12_min)
	    m12_a = m12_min + Epsilon;
	  if (m12_b > m12_max)
	    m12_b = m12_max - Epsilon;

	  h = (m12_b - m12_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  m12 = m12_a;
	  inte_bin_pre = imd_T4_12_extra (m12,
					  beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R,
					  g1_R, phi1_R, h1_R,
					  d8, dE, c8, c10, Lambda,
					  &inte_bin_peak_pre,
					  &inte_bin_reflect_pre,
					  &inte_bin_nonres_pre);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      m12 += h;

	      inte_bin = imd_T4_12_extra (m12, beta_T,
					  gz_R, phiz_R,
					  g0_R, phi0_R,
					  g1_R, phi1_R,
					  h1_R,
					  d8, dE, c8, c10, Lambda,
					  &inte_bin_peak,
					  &inte_bin_reflect,
					  &inte_bin_nonres);
     
	      I_bin += (inte_bin + inte_bin_pre) * h / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      I_bin_peak +=
		(inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;
	      I_bin_nonres +=
		(inte_bin_nonres + inte_bin_nonres_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      inte_bin_nonres_pre = inte_bin_nonres;
	    }

	  events_th = I_bin / bin;

	  events_peak = I_bin_peak / bin;
	  events_reflect = I_bin_reflect / bin;
	  events_nonres = I_bin_nonres / bin;
      
	  nth = events_th + background_data[i];
      
	  printf ("%f , \t%f\n", m12_data[i], nth);
	  fprintf (out, "%f , \t%f, \t%f, \t%f, \t%f, \t%f\n",
		   m12_data[i]/1000.0, nth,
		   events_peak, events_reflect, events_nonres,
		   background_data[i]/events_nonres);
	}

      fclose (out);
      
      /* Fourth, the histogram for the proton helicity distribution
       */
      out = fopen ("Mz-fit12z-T4-PD-histo.plot", "w");

      fprintf (out,
	       "# beta = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n",
	      1e2*beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      fprintf (out, "# d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);
      
      fprintf (out, "# chisq = %f;\n", chisq);
      
      for (i=0; i < ndata; i++)
	{
	  int j;

	  double z_a, z_b, z_i;
	  double z_bin = 2.0 / ((double) (ndata)), bin = 13.5;

	  double I_bin = 0.0;
	  double inte_bin_pre, inte_bin;

	  double I_bin_peak = 0.0;
	  double inte_bin_peak_pre, inte_bin_peak;

	  double I_bin_reflect = 0.0;
	  double inte_bin_reflect_pre, inte_bin_reflect;
	  
	  double events_z_th, nth_z;
	  double events_z_reflect, events_z_peak;

	  z_min = -1.0;
	  z_max = +1.0;

	  /* Here the bins are evenly spaced between z_min and z_max, and
	     it might be better to calculate them directly like this:
	   */
	  z_i = z_min + ((double) (2*i+1)) / ((double) (ndata));

	  z_a = z_i - z_bin/2.0;
	  z_b = z_i + z_bin/2.0;

	  /* Check that we do not under/overflow z:
	   */
	  if (z_a < z_min)
	    z_a = z_min;
	  if (z_b > z_max)
	    z_b = z_max;

	  h_z = (z_b - z_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
      
	  z = z_a;
	  inte_bin_pre = imd_T4_z (z,
				   beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R, h1_R,
				   d8, dE, c8, c10, Lambda);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      z += h_z;

	      inte_bin = imd_T4_z (z, beta_T, gz_R, phiz_R,
				   g0_R, phi0_R, g1_R, phi1_R,
				   h1_R, d8, dE, c8, c10, Lambda);
     
	      I_bin += (inte_bin + inte_bin_pre) * h_z / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      /*
	      I_bin_peak += (inte_bin_peak + inte_bin_peak_pre) * h / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      */
	    }

	  /* The events are the average over the bin, but then we have to
	     normalize them by the factor z_bib/bin so they give
	     the correct counts:
	   */
	  events_z_th = (I_bin / z_bin) * z_bin / bin;

	  /*
	  events_z_peak = I_bin_peak / bin;
	  events_z_reflect = I_bin_reflect / bin;
	  */
	  
	  nth_z = events_z_th + background_z_data[i];
      
	  printf ("%f , \t%f\n", z_i, nth_z);
	  fprintf (out, "%f , \t%f\n", z_i, nth_z);
	}

      fclose (out);

      /* Finally, we print the production/coupling parameters
	 as well as the contact-range couplings:
       */

      printf ("\n#############\n\n");
	    
      printf ("beta_T = %fe-2; gz_R = %fe-8; phiz_R = %f; g0_R = %fe-8; phi0_R = %f; g1_R = %fe-8; phi1_R = %f; h1_R = %f;\n\n",
	      1e2 * beta_T, 1e8*gz_R, phiz_R, 1e8*g0_R, phi0_R, 1e8*g1_R, phi1_R, h1_R);

      printf ("d8 = %f/hbc/hbc; dE = %f/hbc/hbc; c8 = %f/hbc/hbc; c10 = %f/hbc/hbc;\n\n",
	      hbc*hbc*d8, hbc*hbc*dE, hbc*hbc*c8, hbc*hbc*c10);

      printf ("chisq = %f\n;", chisq);
      
    }
  
  return 0;
}
