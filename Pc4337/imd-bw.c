#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

extern const double M_proton;
extern const double M_JPsi;

extern const double M_Bs0;

extern const double Epsilon;

#define N_INT 50

/* The BW amplitude for a negative parity, spin-1/2 pentaquark:
 */
double complex amp_bw (double s, double M_R, double Gamma_R);

/* The Kallen function:
 */
double
lambda_kallen (double a, double b, double c)
{
  double val = a*a + b*b + c*c - 2.0*a*b -2.0*a*c - 2.0*b*c;

  return val;
}


/* The invariant mass distribution for a Breit-Wigner for a spin 1/2^-
   resonance, which implies that both production and decay
   happen in S-wave.
 */
double
imd_bw_s (double m12, double beta, double g, double phi,
	   double M_R, double Gamma_R)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;

  /* Then, I determine m23_sq_min and m23_sq_max
     (m23_sq is the square of m23; for details one may check
     the PDG review on kinematics):
   */
  double E2 = (m12*m12 - m1*m1 + m2*m2) / 2.0 / m12;
  double E3 = (M*M - m12*m12 - m3*m3) / 2.0 / m12;

  double m12_sq = m12*m12;

  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;
  
  double m23_sq;
  double m23_sq_min = (E2+E3)*(E2+E3) -
    pow(sqrt(E2*E2-m2*m2) + sqrt(E3*E3-m3*m3),2);
  double m23_sq_max = (E2+E3)*(E2+E3) -
    pow(sqrt(E2*E2-m2*m2) - sqrt(E3*E3-m3*m3),2);

  /* After this, I prepare the integration over m23_sq:

     i) amp_12, amp_23 are the intermediate amplitudes
     (e.g. a pentaquark modeled as a Breit-Wigner or
     a meson-baryon loop).

     ii) inte_amp is the amplide to be integrated, before squaring it.

     iii) inte_amp_sq, inte_amp_sq_pre are the integrands
     (I am using trapezoidal)

     iv) I_amp is the integral ($\int |A|^2 d m_{23}^2$
     in LaTeX notation)
   */
  int i;
  
  double complex amp_12, amp_23;

  double I_amp = 0.0;
  double complex inte_amp;
  double inte_amp_sq, inte_amp_sq_pre;
  double h;

  m23_sq = m23_sq_min;
  h = (m23_sq_max - m23_sq_min) / ((double) (N_INT-1));

  amp_12 = amp_bw (m12_sq, M_R, Gamma_R);
  amp_23 = amp_bw (m23_sq, M_R, Gamma_R);
  
  inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);

  inte_amp_sq_pre = pow(cabs(inte_amp),2);

  /* I integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m23_sq += h;

      amp_12 = amp_bw (m12_sq, M_R, Gamma_R);
      amp_23 = amp_bw (m23_sq, M_R, Gamma_R);

      
      inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);

      inte_amp_sq = pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;
    }

  /* I still have to multiply the result by m12 because
     what I am calculating is d \Gamma / d m_12.

     I include the averaging over the range of available m12, but
     this is simply an arbitrary normalization choice.
   */
  I_amp = I_amp * (2.0 * m12 / (m12_max + m12_min));
  
  return I_amp;
}

/* The invariant mass distribution for a Breit-Wigner: this version of
   the program returns 
 */
double
imd_bw_extra (double m12, double beta, double g, double phi,
	      double M_R, double Gamma_R, double *I_peak,
	      double *I_reflection, double *I_nonres)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;

  /* Then, I determine m23_sq_min and m23_sq_max
     (m23_sq is the square of m23; for details one may check
     the PDG review on kinematics):
   */
  double E2 = (m12*m12 - m1*m1 + m2*m2) / 2.0 / m12;
  double E3 = (M*M - m12*m12 - m3*m3) / 2.0 / m12;

  double m12_sq = m12*m12;

  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;
  
  double m23_sq;
  double m23_sq_min = (E2+E3)*(E2+E3) -
    pow(sqrt(E2*E2-m2*m2) + sqrt(E3*E3-m3*m3),2);
  double m23_sq_max = (E2+E3)*(E2+E3) -
    pow(sqrt(E2*E2-m2*m2) - sqrt(E3*E3-m3*m3),2);

  /* After this, I prepare the integration over m23_sq:

     i) amp_12, amp_23 are the intermediate amplitudes (e.g. a pentaquark
        modeled as a Breit-Wigner or a meson-baryon loop).

     ii) inte_amp is the amplide to be integrated, before squaring it

     iii) inte_amp_sq, inte_amp_sq_pre are the integrands
     (I am using trapezoidal)

     iv) I_amp is the integral (i.e. $\int |A|^2 d m_{23}^2$
     in LaTeX notation)
   */
  int i;
  
  double complex amp_12, amp_23;

  double I_amp = 0.0;
  double complex inte_amp;
  double inte_amp_sq, inte_amp_sq_pre;
  double h;

  double I_amp_peak = 0.0;
  double complex inte_amp_peak;
  double inte_amp_peak_sq, inte_amp_peak_sq_pre;
  
  double I_amp_reflect = 0.0;
  double complex inte_amp_reflect;
  double inte_amp_reflect_sq, inte_amp_reflect_sq_pre;
  
  double I_amp_nonres = 0.0;
  double complex inte_amp_nonres;
  double inte_amp_nonres_sq, inte_amp_nonres_sq_pre;
  
  m23_sq = m23_sq_min;
  h = (m23_sq_max - m23_sq_min) / ((double) (N_INT-1));

  amp_12 = amp_bw (m12_sq, M_R, Gamma_R);
  amp_23 = amp_bw (m23_sq, M_R, Gamma_R);
  
  inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);

  inte_amp_sq_pre = pow(cabs(inte_amp),2);

  /* Integrands for the peak, reflection and non-resonant contribs
     (by themselves, in isolation):
   */ 
  inte_amp_peak_sq_pre = g*g * pow(cabs(amp_12), 2);
  inte_amp_reflect_sq_pre = g*g * pow(cabs(amp_23), 2);
  inte_amp_nonres_sq_pre = pow(beta, 2);
 
  /* I integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m23_sq += h;

      amp_12 = amp_bw (m12_sq, M_R, Gamma_R);
      amp_23 = amp_bw (m23_sq, M_R, Gamma_R);
      
      inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);

      inte_amp_sq = pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;

      /* Integrands and integrals for the peak, reflection and
	 non-resonant contribs:
      */ 
      inte_amp_peak_sq = g*g * pow(cabs(amp_12), 2);
      inte_amp_reflect_sq = g*g * pow(cabs(amp_23), 2);
      inte_amp_nonres_sq = pow(beta, 2);

      I_amp_peak += (inte_amp_peak_sq + inte_amp_peak_sq_pre)*h/2.0;
      I_amp_reflect +=
	(inte_amp_reflect_sq + inte_amp_reflect_sq_pre)*h/2.0;
      I_amp_nonres +=
	(inte_amp_nonres_sq + inte_amp_nonres_sq_pre)*h/2.0;

      inte_amp_peak_sq_pre = inte_amp_peak_sq;
      inte_amp_reflect_sq_pre = inte_amp_reflect_sq;
      inte_amp_nonres_sq_pre = inte_amp_nonres_sq;
    }

  /* I still have to multiply the result by m12 because
     what I am calculating is d \Gamma / d m_12.

     I include the averaging over the range of available m12, but
     this is simply an arbitrary normalization choice.
   */
  I_amp = I_amp * (2.0 * m12 / (m12_max + m12_min));

  /* For the extra integrals (peak, reflection, non-resonant) I do
     the same change of normalization:
   */
  *I_peak = I_amp_peak * (2.0 * m12 / (m12_max + m12_min));
  *I_reflection = I_amp_reflect * (2.0 * m12 / (m12_max + m12_min));
  *I_nonres = I_amp_nonres * (2.0 * m12 / (m12_max + m12_min));
  
  return I_amp;
}

/* The z (cosinus of theta_p) mass distribution for a Breit-Wigner for a
   spin 1/2^- resonance.
 */
double
imd_bw_z (double z, double beta, double g, double phi,
	  double M_R, double Gamma_R)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;


  /* For this distribution, one has to integrate in m13 (the proton-antiproton
     invariant mass).
   */
  double m13, m13_sq;

  double m13_min = 2.0*M_proton + Epsilon;
  double m13_max = M_Bs0 - M_JPsi - Epsilon;

  double m13_sq_min = m13_min * m13_min;
  double m13_sq_max = m13_max * m13_max;
  
  /* I still need the min and max values of m12 (the J/Psi-proton invariant mass)
     for the relative normalization with respect to the J/Psi-proton invariant
     mass distribution (otherwise it will be incorrect).
   */
  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  /* Then, I have to define the variables for the integration:

     q_u, k_u, s_uz, t_uz (here u indicates dependence wrt m13, z wrt to
     the cosine of the proton helicity angle, while s and t represent
     m12 and m23).
   */
  double q_u, k_u, s_uz, t_uz;
 
  /* After this, I prepare the integration over m13_sq (i.e. u):

     i) amp_12, amp_23 are the intermediate amplitudes (e.g. a pentaquark
        modeled as a Breit-Wigner or a meson-baryon loop).

     ii) inte_amp is the amplide to be integrated, before squaring it

     iii) inte_amp_sq, inte_amp_sq_pre are the integrands
          (I am using trapezoidal)

     iv) I_amp is the integral ($\int q(m12) k(m12) |A|^2 d m_{12}^2$
     in LaTeX notation)
   */
  int i;
  
  double complex amp_12, amp_23;

  double I_amp = 0.0;
  double complex inte_amp;
  double inte_amp_sq, inte_amp_sq_pre;
  double h;

  
  
  m13_sq = m13_sq_min;
  h = (m13_sq_max - m13_sq_min) / ((double) (N_INT-1));

  q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
  k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) / sqrt(m13_sq) / 2.0;
  
  s_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
  t_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;
  
  amp_12 = amp_bw (s_uz, M_R, Gamma_R);
  amp_23 = amp_bw (t_uz, M_R, Gamma_R);
  inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);
  inte_amp_sq_pre = (q_u*k_u) * pow(cabs(inte_amp),2);

  /* I integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m13_sq += h;

      q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
      k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) / sqrt(m13_sq) / 2.0;
  
      s_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
      t_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;

      amp_12 = amp_bw (s_uz, M_R, Gamma_R);
      amp_23 = amp_bw (t_uz, M_R, Gamma_R);
      inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);
  
      inte_amp_sq = (q_u*k_u) * pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;
    }

  /* I still have to divide the result by (m12_max + m12_min)/2.0
     because I want the normalization to be consistent
     with the one I already used for d \Gamma / d m_12.

     Check the program for the J/Psi proton invariant mass distribution
     for further details.
   */
  I_amp = I_amp * (2.0 / (m12_max + m12_min));
    
  return I_amp;
}


/* The z (cosinus of theta_p) mass distribution for a Breit-Wigner for a
   spin 1/2^- resonance.
 */
double
imd_bw_z_extra (double z, double beta, double g, double phi,
		double M_R, double Gamma_R, double *I_peak,
		double *I_reflection, double *I_nonres)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;


  /* For this distribution, one has to integrate in m13 (the proton-antiproton
     invariant mass).
   */
  double m13, m13_sq;

  double m13_min = 2.0*M_proton + Epsilon;
  double m13_max = M_Bs0 - M_JPsi - Epsilon;

  double m13_sq_min = m13_min * m13_min;
  double m13_sq_max = m13_max * m13_max;
  
  /* I still need the min and max values of m12 (the J/Psi-proton invariant mass)
     for the relative normalization with respect to the J/Psi-proton invariant
     mass distribution (otherwise it will be incorrect).
   */
  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  /* Then, I have to define the variables for the integration:

     q_u, k_u, s_uz, t_uz (here u indicates dependence wrt m13, z wrt to
     the cosine of the proton helicity angle, while s and t represent
     m12 and m23).
   */
  double q_u, k_u, s_uz, t_uz;
 
  /* After this, I prepare the integration over m13_sq (i.e. u):

     i) amp_12, amp_23 are the intermediate amplitudes (e.g. a pentaquark
        modeled as a Breit-Wigner or a meson-baryon loop).

     ii) inte_amp is the amplide to be integrated, before squaring it

     iii) inte_amp_sq, inte_amp_sq_pre are the integrands
          (I am using trapezoidal)

     iv) I_amp is the integral ($\int q(m12) k(m12) |A|^2 d m_{12}^2$
     in LaTeX notation)
   */
  int i;
  
  double complex amp_12, amp_23;

  double I_amp = 0.0;
  double complex inte_amp;
  double inte_amp_sq, inte_amp_sq_pre;
  double h;

  /* Now the integrals and integrands for the peak, reflection and non-resonant
     contributions considered in isolation:
   */
  double I_amp_peak = 0.0;
  double complex inte_amp_peak;
  double inte_amp_peak_sq, inte_amp_peak_sq_pre;
  
  double I_amp_reflect = 0.0;
  double complex inte_amp_reflect;
  double inte_amp_reflect_sq, inte_amp_reflect_sq_pre;
  
  double I_amp_nonres = 0.0;
  double complex inte_amp_nonres;
  double inte_amp_nonres_sq, inte_amp_nonres_sq_pre;

  /* We initialize the integration:
   */
  m13_sq = m13_sq_min;
  h = (m13_sq_max - m13_sq_min) / ((double) (N_INT-1));

  q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
  k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) / sqrt(m13_sq) / 2.0;
  
  s_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
  t_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;
  
  amp_12 = amp_bw (s_uz, M_R, Gamma_R);
  amp_23 = amp_bw (t_uz, M_R, Gamma_R);
  inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);
  inte_amp_sq_pre = (q_u*k_u) * pow(cabs(inte_amp),2);

  /* We consider now the peak/reflection/non-resonant contributions
     by themselves:
   */
  inte_amp_peak_sq_pre = (q_u*k_u) * g*g * pow(cabs(amp_12), 2);
  inte_amp_reflect_sq_pre = (q_u*k_u) * g*g * pow(cabs(amp_23), 2);
  inte_amp_nonres_sq_pre = (q_u*k_u) * pow(beta, 2);

  /* We integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m13_sq += h;

      q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
      k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) / sqrt(m13_sq) / 2.0;
  
      s_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
      t_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;

      amp_12 = amp_bw (s_uz, M_R, Gamma_R);
      amp_23 = amp_bw (t_uz, M_R, Gamma_R);
      inte_amp = beta + g*(cos(phi) + I*sin(phi))*(amp_12 + amp_23);
  
      inte_amp_sq = (q_u*k_u) * pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;

      
      /* Integrands and integrals for the peak, reflection and
	 non-resonant contribs:
      */ 
      inte_amp_peak_sq = (q_u*k_u) * g*g * pow(cabs(amp_12), 2);
      inte_amp_reflect_sq = (q_u*k_u) * g*g * pow(cabs(amp_23), 2);
      inte_amp_nonres_sq = (q_u*k_u) * pow(beta, 2);

      I_amp_peak += (inte_amp_peak_sq + inte_amp_peak_sq_pre)*h/2.0;
      I_amp_reflect +=
	(inte_amp_reflect_sq + inte_amp_reflect_sq_pre)*h/2.0;
      I_amp_nonres +=
	(inte_amp_nonres_sq + inte_amp_nonres_sq_pre)*h/2.0;

      inte_amp_peak_sq_pre = inte_amp_peak_sq;
      inte_amp_reflect_sq_pre = inte_amp_reflect_sq;
      inte_amp_nonres_sq_pre = inte_amp_nonres_sq;
    }

  /* I still have to divide the result by (m12_max + m12_min)/2.0
     because I want the normalization to be consistent
     with the one I already used for d \Gamma / d m_12.

     Check the program for the J/Psi proton invariant mass distribution
     for further details.
   */
  I_amp = I_amp * (2.0 / (m12_max + m12_min));

  /* For the extra integrals (peak, reflection, non-resonant) I do
     the same change of normalization:
   */
  *I_peak = I_amp_peak * (2.0 / (m12_max + m12_min));
  *I_reflection = I_amp_reflect * (2.0 / (m12_max + m12_min));
  *I_nonres = I_amp_nonres * (2.0 / (m12_max + m12_min));

  return I_amp;
}


double complex
amp_bw (double s, double M_R, double Gamma_R)
{
  double complex amp_val;

  double m12 = sqrt(s);
  double Gamma_m;

  /* This are q(m) and q(M_R) as described
     in the appendix of the LHCb paper:
   */
  double q_m = sqrt((s - (M_JPsi + M_proton)*(M_JPsi + M_proton))*
		    (s - (M_JPsi - M_proton)*(M_JPsi - M_proton)))
    / 2.0 / m12;
  double q_R = sqrt((M_R*M_R - (M_JPsi + M_proton)*(M_JPsi + M_proton))*
		    (M_R*M_R - (M_JPsi - M_proton)*(M_JPsi - M_proton)))
    / 2.0 / M_R;

  /* Now I define Gamma(m):
   */
  Gamma_m = Gamma_R * (M_R / m12) * (q_m / q_R);

  /* And with this I have the amplitude:
   */
  amp_val = 1.0 / (M_R*M_R - s - I*M_R*Gamma_m);

  return amp_val;
}
