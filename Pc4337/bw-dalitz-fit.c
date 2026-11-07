/* Programs for fitting the contact-range couplings:
 */

#include <math.h>

#include <complex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

#define MAX_ITER	5000 	/* Maximum number of iterations. */

#define ACCURACY	1.e-5

#define N_INT_BIN 10

extern const double M_proton;
extern const double M_JPsi;

extern const double M_Bs0;

extern const double Epsilon;

/* Rescaling variables used to begin with a more natural-sized simplex:
 */
const double mul_beta = 1.0e2;
const double mul_g = 1.0e-2;

/* Data structure for the BW fit:
 */
struct min_dalitz_data_bw
{
  int n12, n23;
  
  double *events_data;
  double *background_data;
  
  double beta;
  double g_R;
  double phi_R;
  double M_R;
  double Gamma_R;
};

/* Calculates the event count in a square bin of the Dalitz plot:
 */
double
dalitz_bin_bw_s (double m12_sq_a, double m23_sq_a,
		 double m12_sq_b, double m23_sq_b,
		 double beta, double g, double phi,
		 double M_R, double Gamma_R);

/*****************************************************
 ************** Poisson distribution *****************
 *****************************************************/

double poisson_dalitz_bw (double beta, double g_R, double phi_R, 
			  double M_R, double Gamma_R,
			  int n12, int n23,
			  double events_data[n12][n23],
			  double background_data[n12][n23]);

double poisson_dalitz_bw_wrapper (const gsl_vector *pol, void *params);

/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-
 */
double
dalitz_fit_bw_poisson (double *beta, double *g_R, double *phi_R,
		       double *M_R, double *Gamma_R,
		       double beta_step, double g_step,
		       double phi_step,
		       double M_step, double Gamma_step,
		       int n12, int n23,
		       double events_data[n12][n23],
		       double background_data[n12][n23])
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (5);
  gsl_vector *pol_step = gsl_vector_alloc (5);

  struct min_dalitz_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_function F = {&poisson_dalitz_bw_wrapper, 5, NULL};

  /* Now we can fill in the data structures and complete the declaration of
     the function to minimize.
   */
  data.n12 = n12;
  data.n23 = n23;
  data.events_data = &events_data[0][0];
  data.background_data = &background_data[0][0];

  F.params = &data;

  /* We set the inital points for the simplex:

     (note: I usually rescale beta and g as to improve
     the convergence of the simplex search)
   */ 
  gsl_vector_set (pol_ini, 0, (*beta) * mul_beta);
  gsl_vector_set (pol_step, 0, (beta_step) * mul_beta);
  
  gsl_vector_set (pol_ini, 1, (*g_R) * mul_g);
  gsl_vector_set (pol_step, 1, (g_step) * mul_g);
  
  gsl_vector_set (pol_ini, 2, *phi_R);
  gsl_vector_set (pol_step, 2, phi_step);
  
  gsl_vector_set (pol_ini, 3, *M_R);
  gsl_vector_set (pol_step, 3, M_step);
  
  gsl_vector_set (pol_ini, 4, *Gamma_R);
  gsl_vector_set (pol_step, 4, Gamma_step);

  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);

      printf ("%5d ", iter);

      for (i=0; i < 5; i++)
	{
	  printf ("%.3e ", gsl_vector_get (s->x, i));
	}

      printf ("f() = %.3e size = %.3e\n", s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  *beta = gsl_vector_get (s->x, 0) / mul_beta;
  *g_R = gsl_vector_get (s->x, 1) / mul_g;
  *phi_R = gsl_vector_get (s->x, 2);
  *M_R = gsl_vector_get (s->x, 3);
  *Gamma_R = gsl_vector_get (s->x, 4);

  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  return chi2;
}

/* Now we calculate the poisson deviance of the model against
   the binned Dalitz plot:
 */
double
poisson_dalitz_bw (double beta, double g_R, double phi_R, 
		   double M_R, double Gamma_R,
		   int n12, int n23,
		   double events_data[n12][n23],
		   double background_data[n12][n23])
{
  int i, j;

  double m12_min = M_JPsi + M_proton + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  double m23_min = m12_min;
  double m23_max = m12_max;

  double m12_sq_min = m12_min * m12_min;
  double m12_sq_max = m12_max * m12_max;

  double m23_sq_min = m23_min * m23_min;
  double m23_sq_max = m23_max * m23_max;
  
  double events_th, nth, nexp;
  double poisson_deviance = 0.0;

  double bin_12 = (m12_sq_max - m12_sq_min) / ((double) (n12));
  double bin_23 = (m23_sq_max - m23_sq_min) / ((double) (n23));
  
  double m12_sq_a, m12_sq_b;
  double m23_sq_a, m23_sq_b; 

  for (i=0; i < n12; i++)
    for (j=0; j < n23; j++)
      {
	m12_sq_a = m12_sq_min + bin_12 * ((double) (i));
	m12_sq_b = m12_sq_min + bin_12 * ((double) (i+1));
	
	m23_sq_a = m23_sq_min + bin_23 * ((double) (j));
	m23_sq_b = m23_sq_min + bin_23 * ((double) (j+1));

	events_th = dalitz_bin_bw_s (m12_sq_a, m23_sq_a,
				     m12_sq_b, m23_sq_b,
				     beta, g_R, phi_R,
				     M_R, Gamma_R);

	      
	nth = events_th + background_data[i][j];
	nexp = events_data[i][j];

	/* Just in case, put a floor on the theoretical prediction:
	 */
	if (nth < 1.0e-15)
	  nth = 1.0e-15;
	
	if (nexp > 0.0)
	  poisson_deviance += 2.0* (nth - nexp + nexp * log(nexp/nth));
	else
	  poisson_deviance += 2.0* (nth);
      }

  return poisson_deviance;
}

double
poisson_dalitz_bw_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);
  double Gamma_R = gsl_vector_get (pol, 4);

  struct min_dalitz_data_bw *data =
    (struct min_dalitz_data_bw *) params;

  double (*events_data)[data->n23] =
    (double (*)[data->n23]) (data->events_data);
  double (*background_data)[data->n23] =
    (double (*)[data->n23]) (data->background_data);

  chisq = poisson_dalitz_bw (beta, g_R, phi_R, M_R, Gamma_R,
			     data->n12, data->n23,
			     events_data,
			     background_data);

  return chisq;
}


/******************************************************************
 ************** Poisson distribution: fixed width *****************
 ******************************************************************/

double poisson_dalitz_bw_gamma_wrapper (const gsl_vector *pol, void *params);

/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-

   However, in this case, I leave the resonance width as a fixed
   parameter. The idea is to test the sensitivity of the poisson
   deviance with respect to the width.
 */
double
dalitz_fit_bw_gamma_poisson (double *beta, double *g_R,
			     double *phi_R, double *M_R, 
			     double beta_step, double g_step,
			     double phi_step, double M_step,
			     double Gamma_R,
			     int n12, int n23,
			     double events_data[n12][n23],
			     double background_data[n12][n23])
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (4);
  gsl_vector *pol_step = gsl_vector_alloc (4);

  struct min_dalitz_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 4);

  gsl_multimin_function F = {&poisson_dalitz_bw_gamma_wrapper, 4, NULL};

  /* Now we can fill in the data structures and complete the declaration of
     the function to minimize.
   */
  data.n12 = n12;
  data.n23 = n23;
  data.events_data = &events_data[0][0];
  data.background_data = &background_data[0][0];
  data.Gamma_R = Gamma_R;
  
  F.params = &data;

  /* We set the inital points for the simplex:

     (note: I usually rescale beta and g as to improve
     the convergence of the simplex search)
   */ 
  gsl_vector_set (pol_ini, 0, (*beta) * mul_beta);
  gsl_vector_set (pol_step, 0, (beta_step) * mul_beta);
  
  gsl_vector_set (pol_ini, 1, (*g_R) * mul_g);
  gsl_vector_set (pol_step, 1, (g_step) * mul_g);
  
  gsl_vector_set (pol_ini, 2, *phi_R);
  gsl_vector_set (pol_step, 2, phi_step);
  
  gsl_vector_set (pol_ini, 3, *M_R);
  gsl_vector_set (pol_step, 3, M_step);
  
  gsl_multimin_fminimizer_set (s, &F, pol_ini, pol_step);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, ACCURACY);

      printf ("%5d ", iter);

      for (i=0; i < 4; i++)
	{
	  printf ("%.3e ", gsl_vector_get (s->x, i));
	}

      printf ("f() = %.3e size = %.3e\n", s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  *beta = gsl_vector_get (s->x, 0) / mul_beta;
  *g_R = gsl_vector_get (s->x, 1) / mul_g;
  *phi_R = gsl_vector_get (s->x, 2);
  *M_R = gsl_vector_get (s->x, 3);

  gsl_vector_free (pol_ini);
  gsl_vector_free (pol_step);

  chi2 = s->fval;
  
  gsl_multimin_fminimizer_free (s);

  return chi2;
}

double
poisson_dalitz_bw_gamma_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);

  struct min_dalitz_data_bw *data =
    (struct min_dalitz_data_bw *) params;

  double (*events_data)[data->n23] =
    (double (*)[data->n23]) (data->events_data);
  double (*background_data)[data->n23] =
    (double (*)[data->n23]) (data->background_data);

  chisq = poisson_dalitz_bw (beta, g_R, phi_R, M_R, data->Gamma_R,
			     data->n12, data->n23,
			     events_data,
			     background_data);

  return chisq;
}

