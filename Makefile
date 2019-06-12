#
# Makefile for writing nc file code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

COMP = ifort

FLAGS = -O3 -assume byterecl -convert big_endian -fpic -mcmodel large -qopenmp
NETCDF_DIR = /scratch/lustre/plgmacmuz/local2/netcdf
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ2  = model_sat_plot.o
OBJ4  = model_sat_coper.o
	
all: \
	model_sat_plot\
	model_sat_coper

model_sat_plot: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o $(EXEDIR)/model_sat_plot

model_sat_plot.o: model_sat_plot.f90
	$(COMP) $(FLAGS) $(INCL) -c model_sat_plot.f90

model_sat_coper: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o $(EXEDIR)/model_sat_coper

model_sat_coper.o: model_sat_coper.f90
	$(COMP) $(FLAGS) $(INCL) -c model_sat_coper.f90

clean: 
	/bin/rm *.o *.mod  model_sat_plot
