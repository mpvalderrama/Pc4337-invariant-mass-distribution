# Pc4337-invariant-mass-distribution
It describes the J/Psi p invariant mass distribution in the Bs decay in terms of meson-baryon loops.

It depends on the gsl and cerf libraries. In Debian 13.7 the relevant packages can be installed with


$ sudo apt-get install gsl-bin libsgl28 \
$ sudo apt-get install libcerf2 libcerf-dev

but the specifics may depend on your distribution.

To build the programs do either $ make imd_bw or $ make imd_rescattering : the first fits the invariant mass distribution(s) to a Breit-Wigner and the second to meson-baryon scattering. Edit main-bw.c or main-rescattering.c to change which distribution you want to fit, which can be done easily by modifying the macros #define COMPUTE_XYZ where the XYZ continuation depends on one's own choice (just read the comments).
