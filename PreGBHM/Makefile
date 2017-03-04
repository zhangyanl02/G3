#!/bin/bash
HM          	= .
LHM		= $(HM)
BINPATH 	= $(LHM)
NAME		= SHAWDHM
BINS		= $(BINPATH)/$(NAME)

CFLAGS	= -g
PREFIX	= /opt/bin/

ifeq ($(FC),f77)
	FC=gfortran
endif

ifeq ($(FC),ifort)
	LIB_NETCDF_INC =/opt/netcdf-4.2/include
	LIB_NETCDF_LIB =/opt/netcdf-4.2/lib
	FLAGS          =-fpe0 -traceback -WB #-openmp
	NETCDF_FLAG    =-I$(LIB_NETCDF_INC) -L$(LIB_NETCDF_LIB) -lnetcdff -lnetcdf
endif

ifeq ($(FC),gfortran)
	LIB_NETCDF_INC =/usr/include
	LIB_NETCDF_LIB =/usr/lib
	FLAGS          = -ffpe-trap=invalid,zero,overflow -fbacktrace -fopenmp #-fbounds-check
	NETCDF_FLAG    =-I$(LIB_NETCDF_INC) -L$(LIB_NETCDF_LIB) -lnetcdf -lnetcdff
endif

ifneq (,$(findstring login,$(HOSTNAME)))
	FC=gfortran
	LIB_NETCDF_INC =/public3/home/interior/careeri/huser020/program/include
	LIB_NETCDF_LIB =/public3/home/interior/careeri/huser020/program/lib
	FLAGS          =-fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace -fopenmp
	NETCDF_FLAG    =-I$(LIB_NETCDF_INC) -L$(LIB_NETCDF_LIB) -lnetcdf -lnetcdff
endif


SCRSPATH1	= .


OBJ		= 	$(SCRSPATH1)/gisutil.mod $(SCRSPATH1)/mod_preprocess.mod			
OBJ_O		= 	$(SCRSPATH1)/gisutil.o $(SCRSPATH1)/mod_preprocess.o
			




DEBUG   = -g -Wall
CC	= gcc $(DEBUG)




%.mod: %.F90 
	$(FC) -c $(CFLAGS) $(FLAGS) $(NETCDF_FLAG) $<
#	$(FC) -c -fbounds-check -g -ffpe-trap=invalid,zero,overflow -fbacktrace $<
%.o: %.F90 
	$(FC) -c $(CFLAGS) $(FLAGS) $(NETCDF_FLAG) $<

	
	
all: test

test: $(OBJ) 
	$(FC) -o ./PreGBHM $(OBJ_O) PreGBHM.F90 $(CFLAGS) $(NETCDF_FLAG) $(FLAGS) 
	
clean:
	rm -rf *.o *~ $(OBJ)
	rm ./PreGBHM
	
install:
	sudo cp ./PreGBHM $(PREFIX)
