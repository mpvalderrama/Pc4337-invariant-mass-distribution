#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <cerf.h>

extern const double M_proton;
extern const double M_JPsi;

extern const double M_Bs0;

extern const double Epsilon;

extern const double Pi;

/* Masses of the particles involved in the meson-baryon loops:
 */
const double M_D0 = 1864.84;
const double M_DC = 1869.66;
const double M_DD0 = 2006.85;

const double M_Lambda_c = 2286.46;
const double M_Sigma_c_p = 2452.65;
const double M_Sigma_c_pp = 2453.97;

/* Widths of the particles involved in the meson-baryon loops:
 */
const double Gamma_DD0 = 0.054;
const double Gamma_Sigma_c_p = 2.3;
const double Gamma_Sigma_c_pp = 1.89;

/* Integration steps in the second invariant mass squared (e.g. steps
   in m_23^2 if we are calculation the distribution in m_12).
 */
#define N_INT 100

/* The loop function, using my mathematica convention (i.e. opposite
   sign to the usual convention, one still has to multiply it by
   the reduced mass and divide it by a Pi^2 factor):
 */
double complex IG1st (double complex k, double Lambda);

/* The rescattering amplitude for the four channel meson-baryon loop,
   which involves Lambda_c^+ Dbar0 (uncoupled), Lambda_c^+ Dbar*0,
   Sigma_c^+ Dbar0 and Sigma_c^{++} Dbar0 (three coupled channels).

   Here (gz, phiz) is the combined production-decay coupling for the
   uncoupled channel, while (g0, phi0) is for Lambda_c^+ Dbar*0 and
   (g1, phi1) for Sigma_c^+ Dbar0 - Sigma_c^{++} Dbar0 in the I=1/2
   configuration (there is no coupling to I=3/2 in the B_s^0 decay).

   The meson-baryon couplings are d8 (Lambda_c^+ Dbar0/Dbar0* diagonal),
   dE (Lambda_c^+ Dbar0* to Sigma_c Dbar(I=1/2) non-diagonal),
   c8 (Sigma_c Dbar(I=1/2) diagonal) and c10 (Sigma_c Dbar(I=3/2)
   diagonal, which creeps in thanks to isospin breaking
   in the meson and baryon masses).
 */
double complex amp_T4 (double s,
		       double gz, double phiz,
		       double g0, double phi0,
		       double g1, double phi1, double h1,
		       double d8, double dE, double c8, double c10,
		       double Lambda);


/* The Kallen function:
 */
double
lambda_kallen (double a, double b, double c)
{
  double val = a*a + b*b + c*c - 2.0*a*b -2.0*a*c - 2.0*b*c;

  return val;
}

/* The loop function in the first Riemann sheet, in the normalization
   I use: to obtain the usual G0(E) multiply by -mu/Pi/Pi, with mu
   the reduced mass of the particular two-body channel considered.
 */
double complex
IG1st (double complex k, double Lambda)
{
  double complex I0_val;

  double complex I0_val_a = 2.0 * I * k * Pi * cexp(-2.0*k*k/Lambda/Lambda);
  double complex I0_val_b = sqrt(2.0*Pi) * Lambda;
  double complex I0_val_c = - 2.0 * cexp(-2.0*k*k/Lambda/Lambda) * k * Pi
    * cerfi (sqrt(2.0) * k / Lambda);

  I0_val = (I0_val_a + I0_val_b + I0_val_c) / 4.0;
  
  return I0_val;
}

/* The 12 invariant mass distribution for rescattering, including
   the four relevant meson-baryon channels:
   
   Here (gz, phiz) is the combined production-decay coupling for the
   uncoupled channel, while (g0, phi0) is for Lambda_c^+ Dbar*0 and
   (g1, phi1) for Sigma_c^+ Dbar0 - Sigma_c^{++} Dbar0 in the I=1/2
   configuration (there is no coupling to I=3/2 in the B_s^0 decay).

   The non-resonant contribution is given by beta.
   
   The meson-baryon couplings are d8 (Lambda_c^+ Dbar0/Dbar0* diagonal),
   dE (Lambda_c^+ Dbar0* to Sigma_c Dbar(I=1/2) non-diagonal),
   c8 (Sigma_c Dbar(I=1/2) diagonal) and c10 (Sigma_c Dbar(I=3/2)
   diagonal, which creeps in thanks to isospin breaking
   in the meson and baryon masses).
 */
double
imd_T4_12 (double m12, double beta,
	   double gz, double phiz,
	   double g0, double phi0,
	   double g1, double phi1, double h1,
	   double d8, double dE, double c8, double c10,
	   double Lambda)
{
  /* First, the masses of the decaying particle
     and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;

  /* Then, I determine m23_sq_min and m23_sq_max
     (m23_sq is the square of m23):

     (The relevant kinematics can be found
     in the PDG review on kinematics.)
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

     iv) I_amp is the integral ($\int |A|^2 d m_{23}^2$
     in LaTeX notation)
   */
  int i;
  
  double complex amp_12, amp_23;

  double I_amp = 0.0;
  double complex inte_amp;
  double inte_amp_sq, inte_amp_sq_pre;
  double h;

  double A12;
  double A23;

  /* We set a few of the integration parameters:
   */
  m23_sq = m23_sq_min;
  h = (m23_sq_max - m23_sq_min) / ((double) (N_INT-1));

  /* We prepare the first integration poin:
   */
  amp_12 = amp_T4 (m12_sq, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  amp_23 = amp_T4 (m23_sq, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  inte_amp = beta + (amp_12 + amp_23);
  inte_amp_sq_pre = pow(cabs(inte_amp),2);
  
  /* We integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m23_sq += h;

      amp_12 = amp_T4 (m12_sq, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
      amp_23 = amp_T4 (m23_sq, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);

      inte_amp = beta + (amp_12 + amp_23);
      inte_amp_sq = pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;
    }

  /* I still have to multiply the result by m12 because what I am
     calculating is d \Gamma / d m_12 (instead of d \Gamma / d m_12^2);
     the averaging over (m12_min + m12_max)/2.0 is simply
     an arbitrary normalization choice.
   */
  I_amp = I_amp * (2.0 * m12 / (m12_max + m12_min));
  
  return I_amp;
}


/* The event distribution for z, the cosine of the proton helicity angle.

   The inputs
 */
double
imd_T4_z (double z, double beta, double gz, double phiz,
	  double g0, double phi0, double g1, double phi1, double h1,
	  double d8, double dE, double c8, double c10,
	  double Lambda)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;

  /* For this distribution, one has to integrate in m13
     (the proton-antiproton invariant mass).
   */
  double m13, m13_sq;

  double m13_min = 2.0*M_proton + Epsilon;
  double m13_max = M_Bs0 - M_JPsi - Epsilon;

  double m13_sq_min = m13_min * m13_min;
  double m13_sq_max = m13_max * m13_max;
  
  /* I still need the min and max values of m12
     (the J/Psi-proton invariant mass) for the relative normalization
     with respect to the J/Psi-proton invariant mass distribution
     (otherwise the result will generate an incorrect
     number of counts; the shape will be correct though).
   */
  double m12_min = M_proton + M_JPsi + Epsilon;
  double m12_max = M_Bs0 - M_proton - Epsilon;

  /* Then, I have to define the variables for the integration:

     q_u, k_u, s_uz, t_uz (here u indicates dependence
     with respect to m13, z wrt to the cosine of the proton
     helicity angle, while s and t represent m12 and m23).
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
  
  /* Here I set the integration step and initial point:
   */
  m13_sq = m13_sq_min;
  h = (m13_sq_max - m13_sq_min) / ((double) (N_INT-1));

  q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
  k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) / sqrt(m13_sq) / 2.0;
  
  s_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
  t_uz = (M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;
  
  amp_12 = amp_T4
    (s_uz, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  amp_23 = amp_T4
    (t_uz, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  inte_amp = beta + (amp_12 + amp_23);
  inte_amp_sq_pre = (q_u*k_u) * pow(cabs(inte_amp),2);

  /* Then I integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m13_sq += h;

      q_u = sqrt (m13_sq - (m1+m3)*(m1+m3)) / 2.0;
      k_u = sqrt (lambda_kallen (M*M, m13_sq, m2*m2)) /
	sqrt(m13_sq) / 2.0;
  
      s_uz =
	(M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 + 2.0*q_u*k_u*z;
      t_uz =
	(M*M + m1*m1 + m2*m2 + m3*m3 - m13_sq)/2.0 - 2.0*q_u*k_u*z;

      amp_12 = amp_T4 (s_uz, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
      amp_23 = amp_T4 (t_uz, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);

      inte_amp = beta + (amp_12 + amp_23);
      inte_amp_sq = (q_u*k_u) * pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;
    }

  /* I still have to divide the result by (m12_max + m12_min)/2.0
     because I want the normalization to be consistent
     with the one I already used for d \Gamma / d m_12.

     (Check the program for the J/Psi proton invariant mass
     distribution for further details.)
   */
  I_amp = I_amp * (2.0 / (m12_max + m12_min));
  
  return I_amp;
}


/* The invariant mass distribution for rescattering: "_extra" indicates
   that the program also returns the peak, reflection and non-resonant
   contributions separately, which might be useful for the plots.

   The rest is identical to imd_T4_12.
 */
double
imd_T4_12_extra (double m12, double beta,
		 double gz, double phiz,
		 double g0, double phi0,
		 double g1, double phi1, double h1,
		 double d8, double dE, double c8, double c10,
		 double Lambda,
		 double *I_peak, double *I_reflection, double *I_nonres)
{
  /* First, the masses of the decaying particle and the decaying products:
   */
  double M = M_Bs0;
    
  double m1 = M_proton;
  double m2 = M_JPsi;
  double m3 = M_proton;

  /* Then, I determine m23_sq_min and m23_sq_max
     (m23_sq is the square of m23):
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

     iv) I_amp is the integral ($\int |A|^2 d m_{23}^2$
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
  
  /* I prepare the initialization of the integral: integration step and
     first integration point.
   */

  m23_sq = m23_sq_min;
  h = (m23_sq_max - m23_sq_min) / ((double) (N_INT-1));

  amp_12 = amp_T4 (m12_sq, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  amp_23 = amp_T4 (m23_sq, gz, phiz, g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
  inte_amp = beta + (amp_12 + amp_23);
  inte_amp_sq_pre = pow(cabs(inte_amp),2);

  inte_amp_peak_sq_pre = pow(cabs(amp_12), 2);
  inte_amp_reflect_sq_pre = pow(cabs(amp_23), 2);
  inte_amp_nonres_sq_pre = pow(beta, 2);
  
  /* And now I integrate:
   */
  for (i=0; i < (N_INT-1); i++)
    {
      m23_sq += h;

      amp_12 = amp_T4 (m12_sq, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);
      amp_23 = amp_T4 (m23_sq, gz, phiz,
		       g0, phi0, g1, phi1, h1, d8, dE, c8, c10, Lambda);

      inte_amp = beta + (amp_12 + amp_23);
      inte_amp_sq = pow(cabs(inte_amp),2);
      
      I_amp += (inte_amp_sq + inte_amp_sq_pre)*h/2.0;

      inte_amp_sq_pre = inte_amp_sq;

      inte_amp_peak_sq = pow(cabs(amp_12), 2);
      inte_amp_reflect_sq = pow(cabs(amp_23), 2);
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

  /* I still have to multiply the result by m12 because what I am
     calculating is d \Gamma / d m_12 ; the averaging over
     (m12_min + m_12_max)/2.0 is simply an arbitrary
     normalization choice.
   */
  I_amp = I_amp * (2.0 * m12 / (m12_max + m12_min));

  I_amp_peak = I_amp_peak * (2.0 * m12 / (m12_max + m12_min));
  I_amp_reflect = I_amp_reflect * (2.0 * m12 / (m12_max + m12_min));
  I_amp_nonres = I_amp_nonres * (2.0 * m12 / (m12_max + m12_min));

  *I_peak = I_amp_peak;
  *I_reflection = I_amp_reflect;
  *I_nonres = I_amp_nonres;
  
  return I_amp;
}

/* The rescattering amplitude for the four channel meson-baryon loop,
   which involves Lambda_c^+ Dbar0 (uncoupled), Lambda_c^+ Dbar*0,
   Sigma_c^+ Dbar0 and Sigma_c^{++} Dbar0 (three coupled channels).

   Here (gz, phiz) is the combined production-decay coupling for the
   uncoupled channel, while (g0, phi0) is for Lambda_c^+ Dbar*0 and
   (g1, phi1) for Sigma_c^+ Dbar0 - Sigma_c^{++} Dbar0 in the I=1/2
   configuration (there is no coupling to I=3/2 in the B_s^0 decay).

   The meson-baryon couplings are d8 (Lambda_c^+ Dbar0/Dbar0* diagonal),
   dE (Lambda_c^+ Dbar0* to Sigma_c Dbar(I=1/2) non-diagonal),
   c8 (Sigma_c Dbar(I=1/2) diagonal) and c10 (Sigma_c Dbar(I=3/2)
   diagonal, which creeps in thanks to isospin breaking
   in the meson and baryon masses).
 */
double complex
amp_T4 (double s,
	double gz, double phiz,
	double g0, double phi0, double g1, double phi1, double h1,
	double d8, double dE, double c8, double c10,
	double Lambda)
{
  double complex amp_val = 0.0;

  double m12 = sqrt(s);
  double E12 = m12 - (M_DC + M_Sigma_c_pp);

  double Delta_a = (M_DC + M_Sigma_c_pp) - (M_DD0 + M_Lambda_c);
  double Delta_b = (M_DC + M_Sigma_c_pp) - (M_D0 + M_Sigma_c_p);

  double Delta_z = (M_DC + M_Sigma_c_pp) - (M_D0 + M_Lambda_c);
  
  double mu_a = M_DD0 * M_Lambda_c / (M_DD0 + M_Lambda_c);
  double mu_b = M_D0 * M_Sigma_c_p / (M_D0 + M_Sigma_c_p);
  double mu_c = M_DC * M_Sigma_c_pp / (M_DC + M_Sigma_c_pp);

  double mu_z = M_D0 * M_Lambda_c / (M_D0 + M_Lambda_c);
  
  double complex ka = csqrt (2.0*mu_a*(E12 + Delta_a + I*Gamma_DD0/2.0));
  double complex kb = csqrt (2.0*mu_b*(E12 + Delta_b + I*Gamma_Sigma_c_p/2.0));
  double complex kc = csqrt (2.0*mu_c*(E12 + I*Gamma_Sigma_c_pp/2.0));

  double complex kz = csqrt (2.0*mu_z*(E12 + Delta_z + I*Epsilon));
  
  double complex I0_a = (mu_a/Pi/Pi) * IG1st (ka, Lambda);
  double complex I0_b = (mu_b/Pi/Pi) * IG1st (kb, Lambda);
  double complex I0_c = (mu_c/Pi/Pi) * IG1st (kc, Lambda);

  double complex I0_z = (mu_z/Pi/Pi) * IG1st (kz, Lambda);
  
  double caa = + d8;
  double cab = + sqrt(1.0/3.0) * dE;
  double cac = - sqrt(2.0/3.0) * dE;

  double cba = + sqrt(1.0/3.0) * dE;
  double cbb = + c8/3.0 + 2.0*c10/3.0;
  double cbc = -sqrt(2.0)/3 * (c8-c10);

  double cca = - sqrt(2.0/3.0) * dE;
  double ccb = - sqrt(2.0)/3 * (c8-c10);
  double ccc = + 2.0*c8/3.0 + c10/3.0;
  
  double complex IVG0[9] =
    {1.0 + caa * I0_a, cab * I0_b, cac * I0_c,
     cba * I0_a, 1.0 + cbb * I0_b, cbc * I0_c,
     cca * I0_a, ccb * I0_b, 1.0 + ccc * I0_c};

  double complex det_IVG0 =
    (1.0 + caa * I0_a) * (1.0 + cbb * I0_b) * (1.0 + ccc * I0_c) +
    (cca * I0_a) * (cab * I0_b) * (cbc * I0_c) +
    (cba * I0_a) * (ccb * I0_b) * (cac * I0_c) -
    (cca * I0_a) * (1.0 + cbb * I0_b) * (cac * I0_c) -
    (1.0 + caa * I0_a) * (ccb * I0_b) * (cbc * I0_c) -
    (cba * I0_a) * (cab * I0_b) * (1.0 + ccc * I0_c);

  double complex adj_IVG0[9] =
    {IVG0[4]*IVG0[8] - IVG0[5]*IVG0[7],
     IVG0[2]*IVG0[7] - IVG0[1]*IVG0[8],
     IVG0[1]*IVG0[5] - IVG0[2]*IVG0[4],
     IVG0[5]*IVG0[6] - IVG0[3]*IVG0[8],
     IVG0[0]*IVG0[8] - IVG0[2]*IVG0[6],
     IVG0[2]*IVG0[3] - IVG0[0]*IVG0[5],
     IVG0[3]*IVG0[7] - IVG0[4]*IVG0[6],
     IVG0[1]*IVG0[6] - IVG0[0]*IVG0[7],
     IVG0[0]*IVG0[4] - IVG0[1]*IVG0[3]};

  double complex inv_IVG0[9];
  
  double complex T_CC[9];
  
  double complex VC[9] =
    { d8, + sqrt(1.0/3.0) * dE, - sqrt(2.0/3.0) * dE,
     + sqrt(1.0/3.0) * dE, c8/3.0 + 2.0*c10/3.0, -sqrt(2.0)/3 * (c8-c10),
    - sqrt(2.0/3.0) * dE, -sqrt(2.0)/3 * (c8-c10), 2.0*c8/3.0 + c10/3.0};
  
  double complex g[3] = { g0 * (cos(phi0) + I*sin(phi0)),
    + sqrt(1.0/3.0) * g1 * (cos(phi1) + I*sin(phi1)),
    - sqrt(2.0/3.0) * g1 * (cos(phi1) + I*sin(phi1))};
  
  double complex h[3] =
    {1.0, + sqrt(1.0/3.0) * h1, - sqrt(2.0/3.0) * h1};
  
  double complex G0[9] =
    {-I0_a, 0.0, 0.0,
     0.0, -I0_b, 0.0,
     0.0, 0.0, -I0_c};

  double complex Amp_CC[9];
  double complex GT[9], GTG[9];

  double complex T_z = caa / (1.0 + caa * I0_z);
  double complex Amp_z = - I0_z + I0_z * T_z * I0_z;
  
  int i, j, k;

  /* Now I calculate the inverse of I - V G_0 :
   */
  for (i=0; i < 9; i++)
    {
      inv_IVG0[i] = adj_IVG0[i] / det_IVG0;
    }

  /* And the T-matrix as T = (I - V G_0)^{-1} V
   */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      {
	T_CC[3*i + j] = 0.0;

	for (k=0; k < 3; k++)
	  {
	    T_CC[3*i + j] += inv_IVG0[3*i + k] * VC[3*k + j];
	  }
      }

  /* Still, I have to calculate g*(G0 + G0*TCC*G0)*h :

     (i) I define Amp_CC = G0 + G0*T_CC*G0, and I begin
         with the calculation of G0*T_CC*G0
   */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      {
	Amp_CC[3*i + j] = G0[3*i+ j];
	GT[3*i + j] = 0.0;
	
	for (k=0; k < 3; k++)
	  {
	    GT[3*i + j] += G0[3*i + k] * T_CC[3*k + j];
	  }
      }

  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      {
	GTG[3*i + j] = 0.0;
	
	for (k=0; k < 3; k++)
	  {
	    GTG[3*i + j] += GT[3*i + k] * G0[3*k + j];
	  }

	Amp_CC[3*i + j] += GTG[3*i+ j];
      }

  /* Then, I sandwich it between the production and decay vectors:
   */

  amp_val = 0.0; /* Already done, but just to be sure. */
  
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      {
	amp_val += g[i]*Amp_CC[3*i + j]*h[j];
      }

  /* Then, I add the Lambda_c Dbar0 channel:
   */
  
  amp_val += gz * (cos(phiz) + I*sin(phiz)) * Amp_z;
 
  /* And we are done, just return the value:
   */
  
  return amp_val;
}
