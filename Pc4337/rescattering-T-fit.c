/* Programs for fitting the contact-range couplings: it uses
   the simplex algorithm as implemented in the gsl library.
 */

#include <math.h>

#include <complex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

#define MAX_ITER	5000 	/* Maximum number of iterations. */

#define ACCURACY	5.0e-6

/* Integration points within each bin: increase it or decrease it
   according to your own needs.

   Taking into account that the bin size is 13.5 MeV, twenty is enough
   provided the poles are not too narrow. Later on, as a crosscheck,
   the main program generates the bin counts with a different INT_BIN. 
 */
#define N_INT_BIN 20

/* Masses of the decaying particle and decay products:
 */
extern const double M_proton;
extern const double M_JPsi;

extern const double M_Bs0;

extern const double Epsilon;

extern const double hbc;

/* Data structure for the T (rescattering) fit:
 */
struct min_data_T
{
  int n_data; /* Number of bins,
		 assumed to be the same for all distros. */

  /* J\Psi proton data:
   */
  double *m12_data;
  double *events_data;
  double *err_data;
  double *background_data;

  /* J\Psi antiproton data (no *m23_data : assumed to be
     identical to *m12_data):
   */
  double *events_bar_data;
  double *err_bar_data;
  double *background_bar_data;

  /* Proton helicity angle data:
   */
  double *z_data;
  double *events_z_data;
  double *background_z_data;

  /* Non-resonant amplitude:
   */
  double beta;

  /* Combined production-decay couplings:
   */
  double gz;
  double phi_z;
  
  double g0_R;
  double phi0_R;
  
  double g1_R;
  double phi1_R;

  double h1_R;

  /* Meson-baryon couplings (contact-range potential):
   */
  double d8;
  double dE;

  double c8;
  double c10;

  /* Cutoff for the Gaussian regulator:
   */
  double Lambda;
};


/* J/Psi proton (or antiproton) mass distribution:

   (defined in imd-rescattering.c)
 */
double
imd_T4_12 (double m12, double beta, double gz, double phiz,
	   double g0, double phi0, double g1, double phi1, double h1,
	   double d8, double dE, double c8, double c10,
	   double Lambda);

/* Proton helicity distribution:

   (defined in imd-rescattering.c)
 */
double
imd_T4_z (double z, double beta, double gz, double phiz,
	  double g0, double phi0, double g1, double phi1, double h1,
	  double d8, double dE, double c8, double c10,
	  double Lambda);

/**********************************************************
 ************** Poisson distribution (T4) *****************
 **********************************************************/

/* This calculates D for the m12 invariant mass distribution
   assuming that bins are filled as a Poisson distribution
 */
double poisson_T4_12 (double beta,
		      double gz, double phiz,
		      double g0, double phi0,
		      double g1, double phi1, double h1,
		      double d8, double dE, double c8, double c10,
		      double m12_data[], double events_data[],
		      double background_data[],
		      int n_data, double Lambda);

/* Wrappers for the use of the previous function with the gsl library:
   the default wrapper (poisson_T4_12_wrapper) fits to all parameters,
   including the meson-baryon couplings.

   The second wrapper (poisson_T4_12_pd_wrapper) only fits
   the production/decay combined couplings and
   the non-resonant contribution, while the meson-baryon couplings
   are provided by us.
 */
double poisson_T4_12_wrapper (const gsl_vector *pol, void *params);

double poisson_T4_12_pd_wrapper (const gsl_vector *pol, void *params);

/* Here we write the main program for fitting to the m12 invariant
   mass distribution: the pointers are the initial points (input)
   and the fitted parameters (output), the "_step" variables is
   the initial step for the fit parameters, then we have
   the distribution data and finally the cutoff.
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
		       int n_data, double Lambda)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  /* Multipliers to avoid having a very extreme simplex: this is because
     the natural values of the parameters tend to be orders-of-magnitude
     apart from each other, thus generating a simplex with a very
     extreme shape.

     The solution is to internally rescale by a factor,
     so the parameters are of the same order-of-magnitude:
   */
  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;
  double mul_c = hbc*hbc;

  /* The initial fit point and the step, as a gsl vector (12 components):
   */
  gsl_vector *pol_ini = gsl_vector_alloc (12);
  gsl_vector *pol_step = gsl_vector_alloc (12);

  /* Here we set up the minimizer (a simplex): this is mostly
     standard gsl lingo, one may check their reference
     materials and examples for further details.
   */
  struct min_data_T data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 12);

  gsl_multimin_function F = {&poisson_T4_12_wrapper, 12, NULL};

  /* Here I include a complete history of the simplex exploring
     the function to be minized. It is potentially useful for
     resetting searches and the like:
   */
  FILE *dbg = fopen ("M12-fit-history.txt", "w");
  
  /* Now we can fill in the data structures and complete
     the declaration of the function to minimize.
   */
  data.n_data = n_data;
  
  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

  data.Lambda = Lambda;
  
  F.params = &data;

  /* We set the inital points for the simplex: notice that the order
     is not the one in the parameters of this function, but it is
     slightly altered, the reason being the order in which each
     parameter was included in the previous versions of
     this program.
   */ 
  gsl_vector_set (pol_ini, 0, mul_beta * (*betaT));
  gsl_vector_set (pol_step, 0, mul_beta * (betaT_step));
  
  gsl_vector_set (pol_ini, 1, mul_g * (*g0_R));
  gsl_vector_set (pol_step, 1, mul_g * (g0_step));
  
  gsl_vector_set (pol_ini, 2, *phi0_R);
  gsl_vector_set (pol_step, 2, phi0_step);
  
  gsl_vector_set (pol_ini, 3, mul_g * (*g1_R));
  gsl_vector_set (pol_step, 3, mul_g * g1_step);
  
  gsl_vector_set (pol_ini, 4, *phi1_R);
  gsl_vector_set (pol_step, 4, phi1_step);
  
  gsl_vector_set (pol_ini, 5, *h1_R);
  gsl_vector_set (pol_step, 5, h1_step);
  
  gsl_vector_set (pol_ini, 6, mul_c * (*d8));
  gsl_vector_set (pol_step, 6, mul_c * (d8_step));
  
  gsl_vector_set (pol_ini, 7, mul_c * (*dE));
  gsl_vector_set (pol_step, 7, mul_c * (dE_step));
  
  gsl_vector_set (pol_ini, 8, mul_c * (*c8));
  gsl_vector_set (pol_step, 8, mul_c * (c8_step));
  
  gsl_vector_set (pol_ini, 9, mul_c * (*c10));
  gsl_vector_set (pol_step, 9, mul_c * (c10_step));

  gsl_vector_set (pol_ini, 10, mul_g * (*gz_R));
  gsl_vector_set (pol_step, 10, mul_g * (gz_step));

  gsl_vector_set (pol_ini, 11, *phiz_R);
  gsl_vector_set (pol_step, 11, phiz_step);

  /* Now we set the minimizer:
   */
  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  /* And we begin the minimization loop:
   */
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);

      /* Here the intermediate states of the minimization
	 are written to the debugging file
       */
      fprintf (dbg, "%5d ", iter);

      for (i=0; i < 12; i++)
	{
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      fprintf (dbg, "f() = %.3e size = %.3e\n", s->fval, size);      

      /* And here I choose which values to print on the screen:
	 for this minimization I choose the meson-baryon
	 couplings
       */
      
      printf ("%5d ", iter);

      for (i=6; i < 11; i++)
	{
	  printf ("%.3e ", gsl_vector_get (s->x, i));
	}

      printf ("f() = %.3e size = %.3e\n", s->fval, size);

    }
  while (status == GSL_CONTINUE && iter < max_iter);

  fclose (dbg);

  /* Now, I put the fitted results back in the pointers:
   */
  *betaT = gsl_vector_get (s->x, 0) / mul_beta;

  *g0_R = gsl_vector_get (s->x, 1) / mul_g;
  *phi0_R = gsl_vector_get (s->x, 2);

  *g1_R = gsl_vector_get (s->x, 3) / mul_g;
  *phi1_R = gsl_vector_get (s->x, 4);

  *h1_R = gsl_vector_get (s->x, 5);

  *d8 = gsl_vector_get (s->x, 6) / mul_c;
  *dE = gsl_vector_get (s->x, 7) / mul_c;

  *c8 = gsl_vector_get (s->x, 8) / mul_c;
  *c10 = gsl_vector_get (s->x, 9) / mul_c;

  *gz_R = gsl_vector_get (s->x, 10) / mul_g;
  *phiz_R = gsl_vector_get (s->x, 11);
  
  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  /* And, finally, return the minimum value of D (which is not a chi
     square per se, but I call it like that anyways):
   */
  
  return chi2;
}

/* This is the same as the previous program, but it fits only
   the production/decay couplings for a given meson-baryon
   interaction profile.
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
			  double Lambda)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  /* Multipliers to avoid having a very extreme simplex: this is because
     the natural values of the parameters tend to be orders-of-magnitude
     apart from each other, thus generating a simplex with a very
     extreme shape.

     The solution is to internally rescale by a factor,
     so the parameters are of the same order-of-magnitude:
   */
  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;

  /* The initial fit point and the step, as a gsl vector (8 components):
   */
  gsl_vector *pol_ini = gsl_vector_alloc (8);
  gsl_vector *pol_step = gsl_vector_alloc (8);

  
  /* The initial fit point and the step, as a gsl vector (8 components):
   */
  struct min_data_T data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 8);

  gsl_multimin_function F = {&poisson_T4_12_pd_wrapper, 8, NULL};

  FILE *dbg = fopen ("M12-pd-fit-history.txt", "w");
  
  /* Now we can fill in the data structures and complete
     the declaration of the function to minimize.
   */
  data.n_data = n_data;
  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

  data.d8 = d8;
  data.dE = dE;
  data.c8 = c8;
  data.c10 = c10;
  
  data.Lambda = Lambda;
  
  F.params = &data;

  /* We set the inital points for the simplex:
   */
  gsl_vector_set (pol_ini, 0, mul_beta * (*betaT));
  gsl_vector_set (pol_step, 0, mul_beta * (betaT_step));
  
  gsl_vector_set (pol_ini, 1, mul_g * (*g0_R));
  gsl_vector_set (pol_step, 1, mul_g * (g0_step));
  
  gsl_vector_set (pol_ini, 2, *phi0_R);
  gsl_vector_set (pol_step, 2, phi0_step);
  
  gsl_vector_set (pol_ini, 3, mul_g * (*g1_R));
  gsl_vector_set (pol_step, 3, mul_g * g1_step);
  
  gsl_vector_set (pol_ini, 4, *phi1_R);
  gsl_vector_set (pol_step, 4, phi1_step);
  
  gsl_vector_set (pol_ini, 5, *h1_R);
  gsl_vector_set (pol_step, 5, h1_step);
  
  gsl_vector_set (pol_ini, 6, mul_g * (*gz_R));
  gsl_vector_set (pol_step, 6, mul_g * (gz_step));

  gsl_vector_set (pol_ini, 7, *phiz_R);
  gsl_vector_set (pol_step, 7, phiz_step);

  /* Now we set the minimizer:
   */
  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  /* And we begin the minimization loop:
   */
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);
      
      /* Here the intermediate states of the minimization
	 are written to the debugging file
       */
      fprintf (dbg, "%5d ", iter);

      for (i=0; i < 8; i++)
	{
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      fprintf (dbg, "f() = %.3e size = %.3e\n", s->fval, size);      

      /* And here I choose which values to print on the screen:
	 for this minimization I choose the non-resonant
	 contribution (index 0) and the modulus of the
	 combined production/decay couplings (indices 1,3,6):
       */

      printf ("%5d ", iter);

      printf ("%.3e ", gsl_vector_get (s->x, 0));
      printf ("%.3e ", gsl_vector_get (s->x, 1));
      printf ("%.3e ", gsl_vector_get (s->x, 3));
      printf ("%.3e ", gsl_vector_get (s->x, 6));

      printf ("f() = %.3e size = %.3e\n", s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  fclose (dbg);

  /* Now, I put the fitted results back in the pointers:
   */  
  *betaT = gsl_vector_get (s->x, 0) / mul_beta;

  *g0_R = gsl_vector_get (s->x, 1) / mul_g;
  *phi0_R = gsl_vector_get (s->x, 2);

  *g1_R = gsl_vector_get (s->x, 3) / mul_g;
  *phi1_R = gsl_vector_get (s->x, 4);

  *h1_R = gsl_vector_get (s->x, 5);

  *gz_R = gsl_vector_get (s->x, 6) / mul_g;
  *phiz_R = gsl_vector_get (s->x, 7);
  
  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  /* And, finally, return the minimum value of D (the Poisson deviance,
     which is not a chi square per se, but the variable bears
     that name for historical reasons):
   */

  return chi2;
}

/* This program calculates the Poisson deviance for the m12
   distribution and a given set of parameters.
 */
double
poisson_T4_12 (double beta,
	       double gz, double phiz,
	       double g0, double phi0,
	       double g1, double phi1, double h1,
	       double d8, double dE, double c8, double c10,
	       double m12_data[], double events_data[],
	       double background_data[],
	       int n_data, double Lambda)
{
  int i, j;

  double m12_min = M_JPsi + M_proton;
  double m12_max = M_Bs0 - M_proton;
  
  double m12, events_th, nth, nexp;
  double chisq_val = 0.0;

  double m12_a, m12_b, bin = 13.5;

  double I_bin = 0.0;
  double inte_bin, inte_bin_pre;
  double h;
  
   for (i=0; i < n_data; i++)
    {
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
      inte_bin_pre = imd_T4_12 (m12, beta, gz, phiz, g0, phi0,
				g1, phi1, h1,
				d8, dE, c8, c10, Lambda);
      
      for (j=0; j < (N_INT_BIN-1); j++)
	{
	  m12 += h;

	  inte_bin = imd_T4_12 (m12, beta, gz, phiz, g0, phi0,
				g1, phi1, h1,
				d8, dE, c8, c10, Lambda);
     
	  I_bin += (inte_bin + inte_bin_pre) * h / 2.0;

	  inte_bin_pre = inte_bin;
	}

      events_th = I_bin / bin;
      
      nth = events_th + background_data[i];
      nexp = events_data[i];
      
      chisq_val += 2.0* (nth - nexp + nexp * log(nexp/nth));
    }
  
  return chisq_val;
}

/* gsl wrapper for the Poisson deviance when all parameters
   are fitted.
 */
double
poisson_T4_12_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;
  double mul_c = hbc*hbc;

  double beta = gsl_vector_get (pol, 0) / mul_beta;

  double g0_R = gsl_vector_get (pol, 1) / mul_g;
  double phi0_R = gsl_vector_get (pol, 2);

  double g1_R = gsl_vector_get (pol, 3) / mul_g;
  double phi1_R = gsl_vector_get (pol, 4);

  double h1_R = gsl_vector_get (pol, 5);

  double d8 = gsl_vector_get (pol, 6) / mul_c;
  double dE = gsl_vector_get (pol, 7) / mul_c;
  double c8 = gsl_vector_get (pol, 8) / mul_c;
  double c10 = gsl_vector_get (pol, 9) / mul_c;

  double gz_R = gsl_vector_get (pol, 10) / mul_g;
  double phiz_R = gsl_vector_get (pol, 11);
  
  struct min_data_T *data = (struct min_data_T *) params;
  
  chisq = poisson_T4_12 (beta, gz_R, phiz_R, g0_R, phi0_R,
			 g1_R, phi1_R, h1_R,
			 d8, dE, c8, c10,
			 data->m12_data, data->events_data,
			 data->background_data,
			 data->n_data, data->Lambda);
  
  return chisq;
}

/* gsl wrapper for the Poisson deviance when the meson-baryon couplings
   are excluded from the fitted parameters:
 */
double
poisson_T4_12_pd_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;
  double mul_c = hbc*hbc;

  double beta = gsl_vector_get (pol, 0) / mul_beta;

  double g0_R = gsl_vector_get (pol, 1) / mul_g;
  double phi0_R = gsl_vector_get (pol, 2);

  double g1_R = gsl_vector_get (pol, 3) / mul_g;
  double phi1_R = gsl_vector_get (pol, 4);

  double h1_R = gsl_vector_get (pol, 5);

  double gz_R = gsl_vector_get (pol, 6) / mul_g;
  double phiz_R = gsl_vector_get (pol, 7);
  
  struct min_data_T *data = (struct min_data_T *) params;
  
  chisq = poisson_T4_12 (beta, gz_R, phiz_R, g0_R, phi0_R,
			 g1_R, phi1_R, h1_R,
			 data->d8, data->dE, data->c8, data->c10,
			 data->m12_data, data->events_data,
			 data->background_data,
			 data->n_data, data->Lambda);
  
  return chisq;
}

/*************************************************************
 ************** Poisson distribution (for z) *****************
 ************************************************************/

double poisson_T4_z (double beta,
		   double gz, double phiz,
		   double g0, double phi0,
		   double g1, double phi1, double h1,
		   double d8, double dE, double c8, double c10,
		   double events_z_data[],
		   double background_z_data[],
		   int n_data, double Lambda);

double poisson_T4_12_z_wrapper (const gsl_vector *pol, void *params);

double poisson_T4_12_z_pd_wrapper (const gsl_vector *pol, void *params);


double
imd_fit_T4_12_z_poisson (double *betaT,
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
		    double events_z_data[],
		    double background_z_data[],
		    int n_data, double Lambda)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  /* Multipliers to avoid having a very extreme simplex:
   */
  double mul_beta = 1.0e2;
  double mul_g0 = 1.0e8;
  double mul_g1 = 1.0e8;
  double mul_c0 = hbc*hbc;
  
  gsl_vector *pol_ini = gsl_vector_alloc (12);
  gsl_vector *pol_step = gsl_vector_alloc (12);

  struct min_data_T data;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 12);

  gsl_multimin_function F = {&poisson_T4_12_z_wrapper, 12, NULL};

  FILE *dbg = fopen ("M12-z-fit-history.txt", "w");
  
  /* Now we can fill in the data structures and complete
     the declaration of the function to minimize.
   */
  data.n_data = n_data;
  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;
  
  data.events_z_data = events_z_data;
  data.background_z_data = background_z_data;
  
  data.Lambda = Lambda;
  
  F.params = &data;

  /* We set the inital points for the simplex:
   */
  
  gsl_vector_set (pol_ini, 0, mul_beta * (*betaT));
  gsl_vector_set (pol_step, 0, mul_beta * (betaT_step));
  
  gsl_vector_set (pol_ini, 1, mul_g0 * (*g0_R));
  gsl_vector_set (pol_step, 1, mul_g0 * (g0_step));
  
  gsl_vector_set (pol_ini, 2, *phi0_R);
  gsl_vector_set (pol_step, 2, phi0_step);
  
  gsl_vector_set (pol_ini, 3, mul_g1 * (*g1_R));
  gsl_vector_set (pol_step, 3, mul_g1 * g1_step);
  
  gsl_vector_set (pol_ini, 4, *phi1_R);
  gsl_vector_set (pol_step, 4, phi1_step);
  
  gsl_vector_set (pol_ini, 5, *h1_R);
  gsl_vector_set (pol_step, 5, h1_step);
  
  gsl_vector_set (pol_ini, 6, mul_c0 * (*d8));
  gsl_vector_set (pol_step, 6, mul_c0 * (d8_step));
  
  gsl_vector_set (pol_ini, 7, mul_c0 * (*dE));
  gsl_vector_set (pol_step, 7, mul_c0 * (dE_step));
  
  gsl_vector_set (pol_ini, 8, mul_c0 * (*c8));
  gsl_vector_set (pol_step, 8, mul_c0 * (c8_step));
  
  gsl_vector_set (pol_ini, 9, mul_c0 * (*c10));
  gsl_vector_set (pol_step, 9, mul_c0 * (c10_step));

  gsl_vector_set (pol_ini, 10, mul_g0 * (*gz_R));
  gsl_vector_set (pol_step, 10, mul_g0 * (gz_step));

  gsl_vector_set (pol_ini, 11, *phiz_R);
  gsl_vector_set (pol_step, 11, phiz_step);
  
  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);

      /* Printing everything into the log/debugging file:
       */
      fprintf (dbg, "%5d ", iter);

      for (i=0; i < 8; i++)
	{
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      fprintf (dbg, "f() = %.3e size = %.3e\n", s->fval, size);

      /* Printing part of the output in the terminal:
       */
      printf ("%5d ", iter);
      
      for (i=6; i < 10; i++)
	{
	  printf ("%.3e ", gsl_vector_get (s->x, i));
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      printf ("f() = %.3e size = %.3e\n", s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  fclose (dbg);
  
  *betaT = gsl_vector_get (s->x, 0) / mul_beta;

  *g0_R = gsl_vector_get (s->x, 1) / mul_g0;
  *phi0_R = gsl_vector_get (s->x, 2);

  *g1_R = gsl_vector_get (s->x, 3) / mul_g1;
  *phi1_R = gsl_vector_get (s->x, 4);

  *h1_R = gsl_vector_get (s->x, 5);

  *d8 = gsl_vector_get (s->x, 6) / mul_c0;
  *dE = gsl_vector_get (s->x, 7) / mul_c0;

  *c8 = gsl_vector_get (s->x, 8) / mul_c0;
  *c10 = gsl_vector_get (s->x, 9) / mul_c0;

  *gz_R = gsl_vector_get (s->x, 10) / mul_g0;
  *phiz_R = gsl_vector_get (s->x, 11);
  
  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  return chi2;
}

/* This program fits only the production/decay couplings
   with a given meson-baryon interaction profile.
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
			    double Lambda)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  /* Multipliers to avoid having a very extreme simplex:
   */
  double mul_beta = 1.0e2;
  double mul_g0 = 1.0e8;
  double mul_g1 = 1.0e8;
  double mul_c0 = hbc*hbc;
  
  gsl_vector *pol_ini = gsl_vector_alloc (8);
  gsl_vector *pol_step = gsl_vector_alloc (8);

  struct min_data_T data;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 8);

  gsl_multimin_function F = {&poisson_T4_12_z_pd_wrapper, 8, NULL};

  FILE *dbg = fopen ("M12z-pd-fit-history.txt", "w");
  
  /* Now we can fill in the data structures and complete the declaration of
     the function to minimize.
   */
  data.n_data = n_data;
  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

  data.events_z_data = events_z_data;
  data.background_z_data = background_z_data;
    
  data.d8 = d8;
  data.dE = dE;
  data.c8 = c8;
  data.c10 = c10;
  
  data.Lambda = Lambda;
  
  F.params = &data;

  /* We set the inital points for the simplex:
   */
  
  gsl_vector_set (pol_ini, 0, mul_beta * (*betaT));
  gsl_vector_set (pol_step, 0, mul_beta * (betaT_step));
  
  gsl_vector_set (pol_ini, 1, mul_g0 * (*g0_R));
  gsl_vector_set (pol_step, 1, mul_g0 * (g0_step));
  
  gsl_vector_set (pol_ini, 2, *phi0_R);
  gsl_vector_set (pol_step, 2, phi0_step);
  
  gsl_vector_set (pol_ini, 3, mul_g1 * (*g1_R));
  gsl_vector_set (pol_step, 3, mul_g1 * g1_step);
  
  gsl_vector_set (pol_ini, 4, *phi1_R);
  gsl_vector_set (pol_step, 4, phi1_step);
  
  gsl_vector_set (pol_ini, 5, *h1_R);
  gsl_vector_set (pol_step, 5, h1_step);
  
  gsl_vector_set (pol_ini, 6, mul_g0 * (*gz_R));
  gsl_vector_set (pol_step, 6, mul_g0 * (gz_step));

  gsl_vector_set (pol_ini, 7, *phiz_R);
  gsl_vector_set (pol_step, 7, phiz_step);
  
  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);

      printf ("%5d ", iter);
      fprintf (dbg, "%5d ", iter);

      for (i=0; i < 6; i++)
	{
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      fprintf (dbg, "\n");
      
      for (i=0; i < 4; i++)
	{
	  printf ("%.3e ", gsl_vector_get (s->x, i));
	  fprintf (dbg, "%.3e ", gsl_vector_get (s->x, i));
	}

      printf ("f() = %.3e size = %.3e\n", s->fval, size);
      fprintf (dbg, "f() = %.3e size = %.3e\n", s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  fclose (dbg);
  
  *betaT = gsl_vector_get (s->x, 0) / mul_beta;

  *g0_R = gsl_vector_get (s->x, 1) / mul_g0;
  *phi0_R = gsl_vector_get (s->x, 2);

  *g1_R = gsl_vector_get (s->x, 3) / mul_g1;
  *phi1_R = gsl_vector_get (s->x, 4);

  *h1_R = gsl_vector_get (s->x, 5);

  *gz_R = gsl_vector_get (s->x, 6) / mul_g0;
  *phiz_R = gsl_vector_get (s->x, 7);
  
  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  return chi2;
}

double poisson_T4_z (double beta,
		     double gz, double phiz,
		     double g0, double phi0,
		     double g1, double phi1, double h1,
		     double d8, double dE, double c8, double c10,
		     double events_z_data[],
		     double background_z_data[],
		     int n_data, double Lambda)
{
  int i, j;

  double z;
  double chisq_val = 0.0;

  double nexp_z;

  double h_z;
  double z_min, z_max;
  
  for (i=0; i < n_data; i++)
    {
      int j;

      double z_a, z_b, z_i;
      double z_bin = 2.0 / ((double) (n_data)), bin = 13.5;

      double I_bin = 0.0;
      double inte_bin_pre, inte_bin;
	  
      double events_z_th, nth_z;

      z_min = -1.0;
      z_max = +1.0;

      /* Here the bins are evenly spaced between z_min and z_max, and
	 it might be better to calculate them directly like this:
      */
      z_i = z_min + ((double) (2*i+1)) / ((double) (n_data));

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
			       beta, gz, phiz,
			       g0, phi0, g1, phi1, h1,
			       d8, dE, c8, c10, Lambda);
	  
      for (j=0; j < (N_INT_BIN-1); j++)
	{
	  z += h_z;

	  inte_bin = imd_T4_z (z, beta, gz, phiz,
			       g0, phi0, g1, phi1,
			       h1, d8, dE, c8, c10, Lambda);
     
	  I_bin += (inte_bin + inte_bin_pre) * h_z / 2.0;
	      
	  inte_bin_pre = inte_bin;
	}

      /* The events are the average over the bin, but then we have to
	 normalize them by the factor z_bib/bin so they give
	 the correct counts:
      */
      events_z_th = (I_bin / z_bin) * z_bin / bin;
	  
      nth_z = events_z_th + background_z_data[i];
      nexp_z = events_z_data[i];
      
      chisq_val += 2.0* (nth_z - nexp_z + nexp_z * log(nexp_z/nth_z));
    }

  return chisq_val;
};

double
poisson_T4_12_z_pd_wrapper (const gsl_vector *pol, void *params)
{
  double chisq, chisq_s, chisq_z;

  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;
  double mul_c = hbc*hbc;

  double beta = gsl_vector_get (pol, 0) / mul_beta;

  double g0_R = gsl_vector_get (pol, 1) / mul_g;
  double phi0_R = gsl_vector_get (pol, 2);

  double g1_R = gsl_vector_get (pol, 3) / mul_g;
  double phi1_R = gsl_vector_get (pol, 4);

  double h1_R = gsl_vector_get (pol, 5);

  double gz_R = gsl_vector_get (pol, 6) / mul_g;
  double phiz_R = gsl_vector_get (pol, 7);
  
  struct min_data_T *data = (struct min_data_T *) params;
  
  chisq_s = poisson_T4_12 (beta, gz_R, phiz_R,
			   g0_R, phi0_R, g1_R, phi1_R, h1_R,
			   data->d8, data->dE, data->c8, data->c10,
			   data->m12_data, data->events_data,
			   data->background_data,
			   data->n_data, data->Lambda);
  
  chisq_z = poisson_T4_z (beta, gz_R, phiz_R,
			  g0_R, phi0_R, g1_R, phi1_R, h1_R,
			  data->d8, data->dE, data->c8, data->c10,
			  data->events_z_data,
			  data->background_z_data,
			  data->n_data, data->Lambda);

  chisq = chisq_s + chisq_z;
  
  return chisq;
}



double
poisson_T4_12_z_wrapper (const gsl_vector *pol, void *params)
{
  double chisq, chisq_s, chisq_z;

  double mul_beta = 1.0e2;
  double mul_g = 1.0e8;
  double mul_c = hbc*hbc;

  double beta = gsl_vector_get (pol, 0) / mul_beta;

  double g0_R = gsl_vector_get (pol, 1) / mul_g;
  double phi0_R = gsl_vector_get (pol, 2);

  double g1_R = gsl_vector_get (pol, 3) / mul_g;
  double phi1_R = gsl_vector_get (pol, 4);

  double h1_R = gsl_vector_get (pol, 5);

  double d8 = gsl_vector_get (pol, 6) / mul_c;
  double dE = gsl_vector_get (pol, 7) / mul_c;
  double c8 = gsl_vector_get (pol, 8) / mul_c;
  double c10 = gsl_vector_get (pol, 9) / mul_c;

  double gz_R = gsl_vector_get (pol, 10) / mul_g;
  double phiz_R = gsl_vector_get (pol, 11);
  
  struct min_data_T *data = (struct min_data_T *) params;
  
  chisq_s = poisson_T4_12 (beta, gz_R, phiz_R,
			   g0_R, phi0_R, g1_R, phi1_R, h1_R,
			   d8, dE, c8, c10,
			   data->m12_data, data->events_data,
			   data->background_data,
			   data->n_data, data->Lambda);
  
  chisq_z = poisson_T4_z (beta, gz_R, phiz_R, g0_R,
			  phi0_R, g1_R, phi1_R, h1_R,
			  d8, dE, c8, c10,
			  data->events_z_data,
			  data->background_z_data,
			  data->n_data, data->Lambda);

  chisq = chisq_s + chisq_z;
  
  return chisq;
}
