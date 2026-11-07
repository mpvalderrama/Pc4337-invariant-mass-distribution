# Pc4337-invariant-mass-distribution
It describes the J/Psi p invariant mass distribution in the Bs decay in terms of meson-baryon loops.

It depends on the gsl and cerf libraries. In Debian 13.7 the relevant packages can be installed with

```bash
$ sudo apt-get install gsl-bin libsgl28 libgsl-dev 
$ sudo apt-get install libcerf2 libcerf-dev
```

but the specifics may depend on your distribution.

To build the programs do either $ make imd_bw or $ make imd_rescattering : the first fits the invariant mass distribution(s) to a Breit-Wigner and the second to meson-baryon rescattering.

There is also $ make dalitz_bw for fitting the binned Dalitz plot with a Breit-Wigner. However, currently there is no direct fit to the binned Dalitz plot for meson-baryon rescattering.

Edit main-bw.c or main-rescattering.c to change which distribution you want to fit, which can be done easily by modifying the macros #define COMPUTE_XYZ where the XYZ continuation depends on one's own choice (just read the comments). However, for meson-baryon rescattering, not all fit options are currently implemented.

The Pc4337.tar.gz file contains all the C code and may be the more convenient option if you don't use github regularly.

The output of these programs are a series of files in .csv format. For plotting them with gnuplot one has to use:


```gnuplot
# For compatibility with the .csv format:
set datafile separator ','  
# This one would plot the M(J/Psi p) mass distribution as a histogram: 
plot [4.03:4.43][0:70] 'M12-fit12-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)"
```

The most practical thing to do here is to go to tutorial/gnuplot/figures-tutorial.gnuplot and check there the different examples provided. This directory also contains the invariant mass distribution histograms from the LHCb manuscript as .csv files: Fig2c.csv for the cos(\theta_p) distro, and Fig2e.csv, Fig2f.csv for the M(J/Psi proton) and M(J/Psi antiproton) distros.

The 'Pc4337-tutorial.pdf' file contains a brief tutorial to the physics and the code used here. Its code is in tutorial/ and for compiling it the first step is to generate the gnuplot figures and then compile the .tex file:


```bash
tutorial/gnuplot$ gnuplot figures-tutorial.gnuplot 
tutorial/gnuplot$ cd .. 
tutorial$ pdflatex Pc4337-tutorial.tex
```