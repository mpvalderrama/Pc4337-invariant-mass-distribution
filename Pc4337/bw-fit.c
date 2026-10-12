/* Programs for fitting the contact-range couplings:
 */

#include <math.h>

#include <complex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

#define MAX_ITER	5000 	/* Maximum number of iterations. */

#define ACCURACY	1.e-6

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
struct min_data_bw
{
  int n_data;
  
  double *m12_data;
  double *events_data;
  double *err_data;
  double *background_data;
  
  double *m23_data;
  double *events_bar_data;
  double *err_bar_data;
  double *background_bar_data;
  
  double *z_data;
  double *events_z_data;
  double *err_z_data;
  double *background_z_data;
  
  double beta;
  double g_R;
  double phi_R;
  double M_R;
  double Gamma_R;
};

double imd_bw_s (double m12, double beta, double g, double phi,
		 double M_R, double Gamma_R);


double imd_bw_z (double z, double beta, double g, double phi,
		 double M_R, double Gamma_R);

/*****************************************************
 ************** Poisson distribution *****************
 *****************************************************/


double poisson_bw (double beta, double g_R, double phi_R, 
		   double M_R, double Gamma_R,
		   double m12_data[], double events_data[],
		   double background_data[],
		   int n_data);

double poisson_bw_wrapper (const gsl_vector *pol, void *params);


/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-
 */
double
imd_fit_12_bw_poisson (double *beta, double *g_R, double *phi_R,
		       double *M_R, double *Gamma_R,
		       double beta_step, double g_step, double phi_step,
		       double M_step, double Gamma_step,
		       double m12_data[], double events_data[],
		       double background_data[],
		       int n_data)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (5);
  gsl_vector *pol_step = gsl_vector_alloc (5);

  struct min_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_function F = {&poisson_bw_wrapper, 5, NULL};

  /* Now we can fill in the data structures and complete the declaration of
     the function to minimize.
   */
  data.n_data = n_data;
  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

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

double
poisson_bw (double beta, double g_R, double phi_R, double M_R,
	    double Gamma_R,
	    double m12_data[], double events_data[],
	    double background_data[],
	    int n_data)
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
      inte_bin_pre = imd_bw_s (m12, beta, g_R, phi_R, M_R, Gamma_R);
      
      for (j=0; j < (N_INT_BIN-1); j++)
	{
	  m12 += h;

	  inte_bin = imd_bw_s (m12, beta, g_R, phi_R, M_R, Gamma_R);
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


double
poisson_bw_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);
  double Gamma_R = gsl_vector_get (pol, 4);

  struct min_data_bw *data = (struct min_data_bw *) params;
  
  chisq = poisson_bw (beta, g_R, phi_R, M_R, Gamma_R,
		      data->m12_data, data->events_data,
		      data->background_data,
		      data->n_data);

  return chisq;
}


/*************************************************************
 ************** Poisson distribution (for z) *****************
 ************************************************************/


double poisson_bw_z (double beta, double g_R, double phi_R, 
		     double M_R, double Gamma_R,
		     double z_data[], double events_z_data[],
		     double background_z_data[],
		     int n_data);

double poisson_bw_z_wrapper (const gsl_vector *pol, void *params);


/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-
 */
double
imd_fit_z_bw_poisson (double *beta, double *g_R, double *phi_R,
		      double *M_R, double *Gamma_R,
		      double beta_step, double g_step, double phi_step,
		      double M_step, double Gamma_step,
		      double z_data[], double events_z_data[],
		      double background_z_data[],
		      int n_data)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (5);
  gsl_vector *pol_step = gsl_vector_alloc (5);

  struct min_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_function F = {&poisson_bw_z_wrapper, 5, NULL};

  /* Now we can fill in the data structures and complete the declaration of
     the function to minimize.
   */
  data.n_data = n_data;
  data.z_data = z_data;
  data.events_z_data = events_z_data;
  data.background_z_data = background_z_data;

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


double poisson_bw_z (double beta, double g_R, double phi_R, 
		     double M_R, double Gamma_R,
		     double z_data[], double events_z_data[],
		     double background_z_data[],
		     int n_data)
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
      inte_bin_pre = imd_bw_z (z, beta, g_R, phi_R, M_R, Gamma_R);
	  
      for (j=0; j < (N_INT_BIN-1); j++)
	{
	  z += h_z;

	  inte_bin = imd_bw_z (z, beta, g_R, phi_R, M_R, Gamma_R);
     
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
poisson_bw_z_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);
  double Gamma_R = gsl_vector_get (pol, 4);

  struct min_data_bw *data = (struct min_data_bw *) params;
  
  chisq = poisson_bw_z (beta, g_R, phi_R, M_R, Gamma_R,
			data->z_data, data->events_z_data,
			data->background_z_data,
			data->n_data);

  return chisq;
}



/******************************************************************
 ************** Compound fits (for 12+z, 12+23+z) *****************
 ******************************************************************/

double poisson_bw_12_z_wrapper (const gsl_vector *pol, void *params);
double poisson_bw_12_23_z_wrapper (const gsl_vector *pol, void *params);

/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-
 */
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
			 int n_data)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (5);
  gsl_vector *pol_step = gsl_vector_alloc (5);

  struct min_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_function F = {&poisson_bw_12_z_wrapper, 5, NULL};

  /* Now we can fill in the data structures and complete
     the declaration of the function to minimize.
   */
  data.n_data = n_data;

  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

  data.z_data = z_data;
  data.events_z_data = events_z_data;
  data.background_z_data = background_z_data;

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


/* Here we write the main program for fitting to a Breit-Wigner
   resonance with J^p = 1/2^-
 */
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
			    int n_data)
{
  int status;
  int i, iter = 0, max_iter = MAX_ITER;

  double size;

  double chi2;

  double k_cm, d;

  gsl_vector *pol_ini = gsl_vector_alloc (5);
  gsl_vector *pol_step = gsl_vector_alloc (5);

  struct min_data_bw data;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_function F = {&poisson_bw_12_23_z_wrapper, 5, NULL};

  /* Now we can fill in the data structures and complete
     the declaration of the function to minimize.
   */
  data.n_data = n_data;

  data.m12_data = m12_data;
  data.events_data = events_data;
  data.background_data = background_data;

  data.m23_data = m12_data;
  data.events_bar_data = events_bar_data;
  data.background_bar_data = background_bar_data;

  data.z_data = z_data;
  data.events_z_data = events_z_data;
  data.background_z_data = background_z_data;

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


double
poisson_bw_12_z_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double chisq_12, chisq_z;
  
  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);
  double Gamma_R = gsl_vector_get (pol, 4);

  struct min_data_bw *data = (struct min_data_bw *) params;
  
  chisq_12 = poisson_bw (beta, g_R, phi_R, M_R, Gamma_R,
			 data->m12_data, data->events_data,
			 data->background_data,
			 data->n_data);

  chisq_z = poisson_bw_z (beta, g_R, phi_R, M_R, Gamma_R,
			  data->z_data, data->events_z_data,
			  data->background_z_data,
			  data->n_data);


  chisq = chisq_12 + chisq_z;
  
  return chisq;
}

double
poisson_bw_12_23_z_wrapper (const gsl_vector *pol, void *params)
{
  double chisq;

  double chisq_12, chisq_23, chisq_z;
  
  double beta = gsl_vector_get (pol, 0) / mul_beta;
  double g_R = gsl_vector_get (pol, 1) / mul_g;
  double phi_R = gsl_vector_get (pol, 2);
  double M_R = gsl_vector_get (pol, 3);
  double Gamma_R = gsl_vector_get (pol, 4);

  struct min_data_bw *data = (struct min_data_bw *) params;
  
  chisq_12 = poisson_bw (beta, g_R, phi_R, M_R, Gamma_R,
			 data->m12_data, data->events_data,
			 data->background_data,
			 data->n_data);

  chisq_23 = poisson_bw (beta, g_R, phi_R, M_R, Gamma_R,
			 data->m12_data, data->events_bar_data,
			 data->background_bar_data,
			 data->n_data);

  chisq_z = poisson_bw_z (beta, g_R, phi_R, M_R, Gamma_R,
			  data->z_data, data->events_z_data,
			  data->background_z_data,
			  data->n_data);

  chisq = chisq_12 + chisq_23 + chisq_z;
  
  return chisq;
}
