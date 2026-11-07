# Pc4337-invariant-mass-distribution
It describes the J/Psi p invariant mass distribution in the Bs decay in terms of meson-baryon loops.

It depends on the gsl and cerf libraries. In Debian 13.7 the relevant packages can be installed with


$ sudo apt-get install gsl-bin libsgl28 libgsl-dev \
$ sudo apt-get install libcerf2 libcerf-dev

but the specifics may depend on your distribution.

To build the programs do either $ make imd_bw or $ make imd_rescattering : the first fits the invariant mass distribution(s) to a Breit-Wigner and the second to meson-baryon rescattering.

There is also $ make dalitz_bw for fitting the binned Dalitz plot with a Breit-Wigner. However, currently there is no direct fit to the binned Dalitz plot for meson-baryon rescattering.

Edit main-bw.c or main-rescattering.c to change which distribution you want to fit, which can be done easily by modifying the macros #define COMPUTE_XYZ where the XYZ continuation depends on one's own choice (just read the comments). However, for meson-baryon rescattering, not all fit options are currently implemented.

The output of these programs are a series of files in .csv format: for plotting them with gnuplot one has to use:

set datafile separator ',' # For compatibility with the format \
plot [4.03:4.43][0:70] 'M12-fit12-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)" # This one would plot the M(J/Psi p) mass distribution as a histogram
