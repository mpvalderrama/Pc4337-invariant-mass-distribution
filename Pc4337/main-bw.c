#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

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

/*
 */
#define N_PLOT 400
#define N_INT_BIN 200

#define COMPUTE_BW_12 1/* Fits a Breit-Wigner to the m12 distro */
#define COMPUTE_BW_23 0 /* Fits a BW to the m23 distro */
#define COMPUTE_BW_Z 0 /* Fits a BW to the z distro */
#define COMPUTE_BW_12_Z 0 /* Fits a BW to the m12 & z distroz */
#define COMPUTE_BW_12_23_Z 0 /* Fits a BW to the m12, m23 and z distros */

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

/* This program fits the BW parameters to the 12 invariant mass
   distribution, assuming the counts are Poisson distributed:
 */
double
imd_fit_12_bw_poisson (double *beta, double *g_R, double *phi_R,
		       double *M_R, double *Gamma_R,
		       double beta_step, double g_step, double phi_step,
		       double M_step, double Gamma_step,
		       double m12_data[], double events_data[],
		       double background_data[],
		       int n_data);

double
imd_fit_z_bw_poisson (double *beta, double *g_R, double *phi_R,
		      double *M_R, double *Gamma_R,
		      double beta_step, double g_step, double phi_step,
		      double M_step, double Gamma_step,
		      double z_data[], double events_z_data[],
		      double background_z_data[],
		      int n_data);

double
imd_fit_12_z_bw_poisson (double *beta, double *g_R, double *phi_R,
			 double *M_R, double *Gamma_R,
			 double beta_step, double g_step,
			 double phi_step,
			 double M_step, double Gamma_step,
			 double m12_data[], double events_data[],
			 double background_data[],
			 double z_data[], double events_z_data[],
			 double background_z_data[],
			 int n_data);

double
imd_fit_12_23_z_bw_poisson (double *beta, double *g_R, double *phi_R,
			    double *M_R, double *Gamma_R,
			    double beta_step, double g_step,
			    double phi_step,
			    double M_step, double Gamma_step,
			    double m12_data[], double events_data[],
			    double background_data[],
			    double events_bar_data[],
			    double background_bar_data[],
			    double z_data[], double events_z_data[],
			    double background_z_data[],
			    int n_data);

int
main (void)
{
  /* The m12 -- M(J/Psi p) -- kinematical variable and its extreme values:
   */
  double m12;

  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  /* This are the integration steps for the m12 and mz distributions.
     
     Here z refers to cos(\theta_p), check Fig. 2 in the LHCb paper.
   */
  double h12 = (m12_max - m12_min)/((double)(N_PLOT - 1));
  int i;

  double z, z_min = -1.0, z_max = 1.0;
  double h_z = (z_max - z_min)/((double)(N_PLOT - 1));
  
  double imd;

  /* Initial parameters for the BW fit:
   */
  double M_R = 4337.0;
  double Gamma_R = 29.0;
  
  double beta_R = 0.6e-2;
  double g_R = 350;
  double phi_R = 0.0;  
  
  double M_step = 5.0;
  double Gamma_step = 2.0;

  double beta_step = 0.1e-2;
  double g_step = 20.0;
  double phi_step = 0.5;

  /* The "chi square" of the Poisson fit (so to speak):
   */
  double chisq;

  /* The data to which we fit: they are stored directly within the program,
     which may not be ideal but is convenient.
   */
  int ndata = 30;

  double events_comb_data[30], err_comb_data[30], background_comb_data[30];
  
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

  /* Combination of the data: this might be used for a combined fit to m12
     and m23 -- M(J/Psi p) and M(J/Psi pbar) -- though I am not sure whether
     this is a good idea from a statistical point of view.

     Use at your own risk, if you want a combined fit to both of them.
   */
  for (i=0; i < ndata; i++)
    {
      events_comb_data[i] = (events_data[i] + events_bar_data[i])/2.0;
      err_comb_data[i] = sqrt(events_comb_data[i]);
      background_comb_data[i] =
	(background_data[i] + background_bar_data[i])/2.0;
    }

  if (COMPUTE_BW_12)
    {
      chisq = imd_fit_12_bw_poisson (&beta_R, &g_R, &phi_R, &M_R, &Gamma_R,
				     beta_step, g_step, phi_step,
				     M_step, Gamma_step,
				     m12_data, events_data,
				     background_data,
				     ndata);

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fit12-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fit12-BW-histo.plot", "w");

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
      out = fopen ("Mz-fit12-BW.plot", "w");
  
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
   
      out = fopen ("Mz-fit12-BW-histo.plot", "w");
      
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

      
      printf ("beta_R = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f;\n",
	      (1.e3)*beta_R, 1.e-2*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }  

  
  if (COMPUTE_BW_23)
    {
      /* Initial parameters:
       */
      beta_R = 3.937931e-3 + 1.5e-3;
      /*
      g_R = 2.772474e2;
      */
      g_R = 0.0;
      phi_R = -1.139021;
      M_R = 4340.914736;
      Gamma_R = 29.202114;

      
      chisq = imd_fit_12_bw_poisson (&beta_R, &g_R, &phi_R, &M_R, &Gamma_R,
				     beta_step, g_step, phi_step,
				     M_step, Gamma_step,
				     m12_data, events_bar_data,
				     background_bar_data,
				     ndata);

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fit23-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fit23-BW-histo.plot", "w");

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
      out = fopen ("Mz-fit23-BW.plot", "w");
  
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
   
      out = fopen ("Mz-fit23-BW-histo.plot", "w");
      
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

      
      printf ("beta = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f\n",
	      (1.e3)*beta_R, 1.e-2*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }  

  if (COMPUTE_BW_Z)
    {
      chisq = imd_fit_z_bw_poisson (&beta_R, &g_R, &phi_R,
				    &M_R, &Gamma_R,
				    beta_step, g_step, phi_step,
				    M_step, Gamma_step,
				    z_data, events_z_data,
				    background_z_data,
				    ndata);

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fitz-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fitz-BW-histo.plot", "w");

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
      out = fopen ("Mz-fitz-BW.plot", "w");
  
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
   
      out = fopen ("Mz-fitz-BW-histo.plot", "w");
      
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

      
      printf ("beta_R = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f;\n",
	      (1.e3)*beta_R, 1.e-2*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }  

  
  if (COMPUTE_BW_12_Z)
    {
      /* Initial parameters:
       */
      beta_R = 3.937931e-3 + 1.5e-3;
      /*
      g_R = 2.772474e2;
      */
      g_R = 0.0;
      phi_R = -1.139021;
      M_R = 4340.914736;
      Gamma_R = 29.202114;


      chisq = imd_fit_12_z_bw_poisson (&beta_R, &g_R, &phi_R,
				       &M_R, &Gamma_R,
				       beta_step, g_step, phi_step,
				       M_step, Gamma_step,
				       m12_data, events_data,
				       background_data,
				       z_data, events_z_data,
				       background_z_data,
				       ndata);

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fit12z-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fit12z-BW-histo.plot", "w");

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
      out = fopen ("Mz-fit12z-BW.plot", "w");
  
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
   
      out = fopen ("Mz-fit12z-BW-histo.plot", "w");
      
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

      
      printf ("beta_R = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f;\n",
	      (1.e3)*beta_R, 1.e-2*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }  

    
  if (COMPUTE_BW_12_23_Z)
    {
      /* Initial parameters:
       */
      beta_R = 3.937931e-3 + 1.5e-3;
      /*
      g_R = 2.772474e2;
      */
      g_R = 0.0;
      phi_R = -1.139021;
      M_R = 4340.914736;
      Gamma_R = 29.202114;


      chisq = imd_fit_12_23_z_bw_poisson (&beta_R, &g_R, &phi_R,
					  &M_R, &Gamma_R,
					  beta_step, g_step, phi_step,
					  M_step, Gamma_step,
					  m12_data, events_data,
					  background_data,
					  events_bar_data,
					  background_bar_data,
					  z_data, events_z_data,
					  background_z_data,
					  ndata);

      /* Finally, we plot our results:

	 First, the shape of the predicted invariant mass distribution,
	 assuming no noise:
       */
      FILE *out = fopen ("M12-fit12z-BW.plot", "w");
  
      for (i=0; i < N_PLOT; i++)
	{
	  m12 = m12_min + h12*i;

	  imd = imd_bw_s (m12, beta_R, g_R, phi_R, M_R, Gamma_R);

	  printf ("%f , \t%f\n", m12, imd);
	  fprintf (out, "%f , \t%f\n", m12/1000.0, imd);
	}

      fclose (out);

      /* Second, the predicted histogram for the J/Psi p mass distribution,
	 which includes the noise and is identical to the J/Psi pbar
	 distribution
       */   
      out = fopen ("M12-fit12z-BW-histo.plot", "w");

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
      out = fopen ("Mz-fit12z-BW.plot", "w");
  
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
   
      out = fopen ("Mz-fit12z-BW-histo.plot", "w");
      
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

      
      printf ("beta_R = %fe-3; g_R = %fe2; phi_R = %f; M_R = %f; Gamma_R = %f;\n",
	      (1.e3)*beta_R, 1.e-2*g_R, phi_R, M_R, Gamma_R);
  
      printf ("chisq = %f;\n", chisq);
    }  


  
  return 0;
}
