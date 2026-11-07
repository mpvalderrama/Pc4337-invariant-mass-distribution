#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

/* The invariant mass distribution data are included in the header
   file below
 */
#include "main-imd.h"

/* The Dalitz plot data are included in the header file below.

   The two most important arrays defined in the header are:
   
   static const double dalitz_counts[18][18];
   static const double dalitz_background[18][18];

   which contain the counts per bin and the background (assuming
   proportionality with the total background: the Dalitz histogram
   has 885 events, which is less than the 907 reported events,
   from which 797 are signal and 110 background)
 */
#include "main-dalitz.h"

/* Parameters for the Dalitz plot histogram: 18x18 grid where
   the ranges for the m12^2 and m23^2 squared invariant masses
   is the same.
   
   (particle 2 is J/Psi, while particles 1 and 3 are a proton
   and antiproton, which have the same mass).

   M2_MIN and M2_MAX indicate the minimum and maximum values of
   m12^2 and m23^2.
 */
#define DALITZ_NBIN 18
#define DALITZ_M2_MIN 16.28274
#define DALITZ_M2_MAX 19.61159

/* Here is where one decides what do you want to calculate:
 */

#define COMPUTE_BW_DALITZ 1  /* Fit to the binned (18x18) Dalitz plot */
#define COMPUTE_BW_GAMMA_DALITZ 0  /* Fit to the binned (18x18) Dalitz
				      plot, but with a fixed
				      resonance width. */
/* The masses of the initial and final particles:
 */
const double M_proton = 938.3;
const double M_JPsi = 3096.9;

const double M_Bs0 = 5366.9;

const double hbc = 197.33;

/* The small Epsilon is used to guarantee that we fall
   within kinematical ranges:
 */
const double Epsilon = 1.0e-4;

const double Pi = 3.141592;

/* N_INT_DALITZ_BIN is the number of integration points per bin
   for the Dalitz histogram.
 */
#define N_INT_DALITZ_BIN 50

/* N_PLOT is the number of plot points for the invariant mass
   distribution as a continuous function of the invariant mass,

   N_INT_BIN is the number of integration points per bin when
   the invariant mass distribution is plotted as a histogram.
 */
#define N_PLOT 400
#define N_INT_BIN 200

/* The Breit-Wigner amplitude:
 */
double complex amp_bw (double m12, double M_R, double Gamma_R);

/* The invariant mass distribution for a Breit-Wigner
 */
double imd_bw_s (double m12, double beta, double g, double phi,
		 double M_R, double Gamma_R);

/* Same as before, but now includes the explicit contributions from
   the peak, the reflection and the non-resonant amplitude.
 */
double imd_bw_extra (double m12, double beta, double g, double phi,
		     double M_R, double Gamma_R,
		     double *I_peak, double *I_reflection,
		     double *I_nonres);

/* The z distribution for a Breit-Wigner:
 */ 
double
imd_bw_z (double z, double beta, double g, double phi,
	  double M_R, double Gamma_R);

double
imd_bw_z_extra (double z, double beta, double g, double phi,
		double M_R, double Gamma_R, double *I_peak,
		double *I_reflection, double *I_nonres);


/* This program fits the BW parameters to the binned Dalitz plot,
   assuming the counts are Poisson distributed:
 */
double
dalitz_fit_bw_poisson (double *beta, double *g_R, double *phi_R,
		       double *M_R, double *Gamma_R,
		       double beta_step, double g_step, double phi_step,
		       double M_step, double Gamma_step,
		       int n12, int n23,
		       double events_data[n12][n23],
		       double background_data[n12][n23]);

/* This program fits the BW parameters -- except for the width -- to
   the binned Dalitz plot, assuming the counts are Poisson
   distributed.

   The point is to check the sensitivity of the fit to the width.
 */
double
dalitz_fit_bw_gamma_poisson (double *beta, double *g_R,
			     double *phi_R, double *M_R, 
			     double beta_step, double g_step,
			     double phi_step, double M_step,
			     double Gamma_R,
			     int n12, int n23,
			     double events_data[n12][n23],
			     double background_data[n12][n23]);

int
main (void)
{
  /* The m12 -- M(J/Psi p) -- kinematical variable and its extreme
     values (notice that this is identical for m23):
   */
  double m12;

  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  double m23;

  double m23_min = M_proton + M_JPsi + Epsilon;
  double m23_max = M_Bs0 - M_proton - Epsilon;

  /* This are the integration steps for the m12 and mz distributions.
     
     Here z refers to cos(\theta_p), check Fig. 2 in the LHCb paper.
   */
  double h12 = (m12_max - m12_min)/((double)(N_PLOT - 1));
  int i;

  double z, z_min = -1.0, z_max = 1.0;
  double h_z = (z_max - z_min)/((double)(N_PLOT - 1));
  
  double imd;

  /* The number of data for the histograms of the different
     invariant mass distributions.
   */
  int ndata = 30;

  /* Initial parameters for the BW fit:
   */
  double M_R = 4340.9;
  double Gamma_R = 29.2;
  
  double beta_R = 0.793e-2;
  double g_R = 1.2;
  double phi_R = +0.13;  
  
  double M_step = 10.0;
  double Gamma_step = 1.0;

  double beta_step = 0.1e-2;
  double g_step = 1.0;
  double phi_step = 0.1;

  /* The "chi square" of the Poisson fit (so to speak):
   */
  double chisq;

  /* The grid size for the binned Dalitz plot (18 x 18, where
     DALITZ_NBIN is actually equal 18):
   */
  int n_dalitz_data = DALITZ_NBIN;
  
  /* m12_data is the invariant mass for the center of the bin
     in the invariant mass histogram for m12 (J/Psi proton).

     It is defined in the "main-imd.h" header, and it is originally
     written in GeV. Here I change m12 to MeV:
   */
  for (i=0; i < ndata; i++)
    {
      m12_data[i] = 1000.0 * m12_data[i];
    }

  if (COMPUTE_BW_DALITZ)
    {
      chisq = dalitz_fit_bw_poisson (&beta_R, &g_R, &phi_R, &M_R, &Gamma_R,
				     beta_step, g_step, phi_step,
				     M_step, Gamma_step,
				     n_dalitz_data,
				     n_dalitz_data,
				     dalitz_counts,
				     dalitz_background);

      /* After this, we may generate the predicted
	 invariant mass distributions:
       */

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fit-dalitz-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);
	  
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fit-dalitz-BW-histo.plot", "w");

      /* We begin by printing the fitted parameters, for later use: */
      fprintf (out,
	       "# beta_R = %f; g_R = %f; phi_R = %f; M_R = %f; Gamma_R = %f\n",
	       beta_R, g_R, phi_R, M_R, Gamma_R);
  
      fprintf (out, "# chisq = %f;\n", chisq);

      fprintf (out, "# bin(center), calculated counts, peak, reflection, non-resonant, noise/non-resonant \n\n");
      
      /* And now we begin the bin integration: */
      for (i=0; i < ndata; i++)
	{
	  int j;
	  
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
	    m12_a = m12_min;
	  if (m12_b > m12_max)
	    m12_b = m12_max;

	  h12 = (m12_b - m12_a)/((double) (N_INT_BIN -1));

	  I_bin = 0.0;
	  
	  m12 = m12_a;
	  inte_bin_pre = imd_bw_extra (m12, beta_R, g_R, phi_R,
				       M_R, Gamma_R,
				       &inte_bin_peak_pre,
				       &inte_bin_reflect_pre,
				       &inte_bin_nonres_pre);
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      m12 += h12;

	      inte_bin = imd_bw_extra (m12, beta_R, g_R, phi_R,
				       M_R, Gamma_R,
				       &inte_bin_peak,
				       &inte_bin_reflect,
				       &inte_bin_nonres);
	  
	      I_bin += (inte_bin + inte_bin_pre) * h12 / 2.0;
	      
	      inte_bin_pre = inte_bin;

	      I_bin_peak +=
		(inte_bin_peak + inte_bin_peak_pre) * h12 / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h12 / 2.0;
	      I_bin_nonres +=
		(inte_bin_nonres + inte_bin_nonres_pre) * h12 / 2.0;

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

      /* Third, the proton helicity distribution:
       */
      out = fopen ("Mz-fit-dalitz-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  z = z_min + h_z*i;

	  imd = imd_bw_z (z, beta_R, g_R, phi_R,
			  M_R, Gamma_R);
			  
	  /* I have to normalize for the differences in bin sizes:
	   */
	  imd = imd / 13.5 / 15.0;
	  
	  printf ("%f , \t%f\n", z, imd);
	  fprintf (out, "%f , \t%f\n", z, imd);
	}

      fclose (out);


      /* Fourth, the histogram for the proton helicity distribution
       */
   
      out = fopen ("Mz-fit-dalitz-BW-histo.plot", "w");
      
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

	  double I_bin_nonres = 0.0;
	  double inte_bin_nonres_pre, inte_bin_nonres;
	  
	  double events_z_th, nth_z;
	  double events_z_reflect, events_z_peak, events_z_nonres;

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
	  inte_bin_pre = imd_bw_z_extra (z, beta_R, g_R, phi_R,
					 M_R, Gamma_R,			 
					 &inte_bin_peak_pre,
					 &inte_bin_reflect_pre,
					 &inte_bin_nonres_pre);
     
	  
	  for (j=0; j < (N_INT_BIN-1); j++)
	    {
	      z += h_z;

	      inte_bin = imd_bw_z_extra (z,beta_R, g_R, phi_R,
					 M_R, Gamma_R,
					 &inte_bin_peak,
					 &inte_bin_reflect,
					 &inte_bin_nonres);
	  
	      I_bin += (inte_bin + inte_bin_pre) * h_z / 2.0;
	      
	      inte_bin_pre = inte_bin;
	      
	      I_bin_peak +=
		(inte_bin_peak + inte_bin_peak_pre) * h_z / 2.0;
	      I_bin_reflect +=
		(inte_bin_reflect + inte_bin_reflect_pre) * h_z / 2.0;
	      I_bin_nonres +=
		(inte_bin_nonres + inte_bin_nonres_pre) * h_z / 2.0;

	      inte_bin_peak_pre = inte_bin_peak;
	      inte_bin_reflect_pre = inte_bin_reflect;
	      inte_bin_nonres_pre = inte_bin_nonres;
	    }

	  /* The events are the average over the bin, but then we have to
	     normalize them by the factor z_bin/bin so they give
	     the correct counts:

	     (Obviously, the z_bin factors cancel out; I have left
	     them written down explicitly for clarity though.)
	   */
	  events_z_th = (I_bin / z_bin) * z_bin / bin;

	  events_z_peak = (I_bin_peak / z_bin) * z_bin / bin;
	  events_z_reflect = (I_bin_reflect / z_bin) * z_bin / bin;
	  events_z_nonres = (I_bin_nonres / z_bin) * z_bin / bin;
	  
	  nth_z = events_z_th + background_z_data[i];

	  /*
	  printf ("%f , \t%f\n", z_i, nth_z);
	  */
	  fprintf (out, "%f , \t%f, \t%f, \t%f, \t%f, \t%f\n",
		   z_i, nth_z,
		   events_z_peak, events_z_reflect, events_z_nonres,
		   background_z_data[i]/events_z_nonres);
	}

      fclose (out);

      /* Finally, we print the parameters:
       */
      printf ("beta_R = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f;\n",
	      (1.e3)*beta_R, (1.e-2)*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }

  /* For calculating the dependence of the Poisson deviance on
     the resonance width, I find it most useful to start at
     Gamma = 150 MeV and then slowly reducing the width.

     Going from Gamma = 10 MeV upwards generates a discontinuity.
   */
  
  if (COMPUTE_BW_GAMMA_DALITZ)
    {
      double Gamma_a = 150.0;
      double Gamma_b = 10.0;

      double h_Gamma = -1.0;

      FILE *out = fopen ("Poisson-deviance-vs-width.plot", "w");
      
      for (Gamma_R = Gamma_a; Gamma_R > (Gamma_b - h_Gamma/2.0);
	   Gamma_R += h_Gamma)
	{
	  chisq = dalitz_fit_bw_gamma_poisson (&beta_R, &g_R, &phi_R, &M_R, 
					       beta_step, g_step, phi_step,
					       M_step, Gamma_R,
					       n_dalitz_data,
					       n_dalitz_data,
					       dalitz_counts,
					       dalitz_background);

	  printf ("%f, %f ,%f\n", Gamma_R, chisq, chisq / (149-5));
	  fprintf (out, "%f, %f, %f\n", Gamma_R, chisq, chisq / (149-5));
	}

      fclose (out);
    }
  
  return 0;
}
