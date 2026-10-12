# Building rules for the programs fitting the BW and the meson-baryon
# rescattering loops to the relevant invariant mass distributions
# in the B_s^0 decay where the Pc(4337) was found:

CC = gcc
CFLAGS = -g

o_imd_bw = main-bw.o bw-fit.o imd-bw.o

o_imd_rescattering = main-rescattering.o rescattering-T-fit.o \
	imd-rescattering.o 

imd_bw : $(o_imd_bw)
	$(CC) $(CFLAGS) -o imd_bw $(o_imd_bw) -lm -lcerf -lgsl -lgslcblas

imd_rescattering : $(o_imd_rescattering)
	$(CC) $(CFLAGS) -o imd_rescattering $(o_imd_rescattering) \
	-lm -lcerf -lgsl -lgslcblas
