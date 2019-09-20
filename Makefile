#----------------------------------------------------------#
#                    Makefile for Eddy			   #
#----------------------------------------------------------#

# Machine: MIRA

# paths
VPATH       =./post/
ANNPATH     = /opt/ann_1.1.2
HDF5_PATH   = /home/jose/WORK/hdf5-1.8.12/hdf5
LAPACK_PATH = /home/jose/WORK/lapack-3.6.1
BLAS_PATH   = /opt/BLAS-3.6.0


#Compiler command
F_COMP   = /opt/openmpi/bin/mpif90
CC_COMP  = /opt/openmpi/bin/mpicc


#Libraries
ANNLIB     =  ANN
CLIBS      =  -L$(ANNPATH)/lib -l$(ANNLIB) -L$(HDF5_PATH)/lib -lhdf5_fortran -lhdf5 -lstdc++ -lm -lz 
FLIBS      =  -L$(LAPACK_PATH) -llapack -L$(BLAS_PATH)$ -lblas -lstdc++ -lm -lz
#HDF5LIBS  =  -L$(HDF5_PATH)/lib -lhdf5_fortran.so -lstdc++ -lm -lz
#TECLIBS   =  /opt/tecplot10/lib/tecio64.a -lstdc++
TECLIBS    =  -lm /home/jose/WORK/libtecio.a -lstdc++ #point to tecplot as soon as possible


#Include directories
BLASDIRS    = -I/opt/BLAS-3.6.0
LAPACKDIRS  = -I/home/jose/WORK/lapack-3.6.1/SRC


#C++ Section
CPPFLAGS = $(CC_COMP) -c $(OPT) -I$(ANNPATH)/include


#Debugging options
# DEBUG = -check all -warn all,nodec,interfaces -gen-interfaces -traceback -fpe0 -fp-stack-check
# DEBUG =-debug -traceback          # use -debug to get debugging info from gdb
# DEBUG =-p                         # -p lets you use profiling with gprof
# DEBUG = -warn all -warn nointerfaces
# DEBUG = -check all -debug all -g -heap-arrays 100000 -fpe0 -traceback


#Intel Compiler
EXTRA  = -heap-arrays 100000
OPT    = -O3 #-O3/-O0 HIGH/NO OPTIMIZATION
CFLAGS  = $(F_COMP) -c -r8 -132 $(OPT) $(EXTRA) $(DEBUG) -I$(ANNPATH)/include -I$(HDF5_PATH)/include 
LDFLAGS = $(F_COMP) -r8 -132 $(OPT) $(DEBUG) -o 
DLDFLAGS = -lm /home/jose/WORK/libtecio.a -lstdc++
#point to tecplot as soon as possible
EXTRA_MIRA = -L$(LAPACK_PATH) -llapack -L$(BLAS_PATH) -lblas
#the code was not compiling without this extra line

export prefix=/home/jose/WORK/export


##########################################################################
# Solver

EXEC = hybrid.x

WDIR = ./run

CMD = $(WDIR)/$(EXEC)


OBJ =  modules.o Eddy6.o Planes.o mpi_setup.o boundary.o direct.o grid.o\
 info.o inpall.o Initialize.o io.o rt.o matint.o object.o rhs.o setup.o\
 timeindeg.o timestep.o refreshbc.o numrec.o comf.o blktri.o fftpack.o\
 genbun.o gnbnaux.o flowstat.o vt.o calcforce.o triinter.o momforc.o\
 interp.o rbm.o devel.o stats.o output.o output_cart.o pdc2dn.o rhs_density.o\
 density.o hdf5io.o ANN.o hybrid.o


###########################################################################
# Postprocessing


mean := ./run/post_mean.x
energy := ./run/post_energy.x

OBJ_MEAN   = mean.o
OBJ_ENERGY = energy.o

OBJ_PP =  modules.o Planes_post.o mpi_setup.o boundary.o direct.o grid.o\
 info.o inpall.o inpall_post.o Initialize.o io_post.o rt.o matint.o object.o rhs.o setup.o\
 timeindeg.o timestep.o refreshbc.o numrec.o comf.o blktri.o fftpack.o\
 genbun.o gnbnaux.o flowstat.o vt.o calcforce.o triinter.o momforc.o\
 interp.o rbm.o devel.o stats.o output.o output_cart.o pdc2dn.o rhs_density.o\
 density.o hdf5io.o ANN.o


###########################################################################

.SUFFIXES:  .F .f .h

.cpp.o:
	$(CPPFLAGS) $<

.F.o: 
	$(CFLAGS) $<

.f.o: 
	$(CFLAGS) $< 

 
$(CMD): $(OBJ)
	@echo Makefile: ... compiling $@
	$(LDFLAGS) $(@) $(OBJ) ${BLASDIRS} ${LAPACKDIRS} $(CLIBS) $(FLIBS) $(TECLIBS) $(EXTRA_MIRA)\
	 $(DLDFLAGS) 

dep:
	makedepend -fMakefile.par *.F *.f *.f90 *.h


.PHONY : mean
mean : $(mean)

$(mean):   $(OBJ_MEAN) $(OBJ_PP)
	@echo Makefile: ... compiling $
	$(LDFLAGS) $(@) $(OBJ_MEAN) $(OBJ_PP) ${BLASDIRS} ${LAPACKDIRS} $(CLIBS) $(FLIBS) $(TECLIBS) $(EXTRA_MIRA) $(DLDFLAGS) 

 
.PHONY : energy 
energy: $(energy)

$(energy): $(OBJ_ENERGY) $(OBJ_PP)
	@echo Makefile: ... compiling $
	$(LDFLAGS) $(@) $(OBJ_ENERGY) $(OBJ_PP) ${BLASDIRS} ${LAPACKDIRS} $(CLIBS) $(FLIBS) $(TECLIBS) $(EXTRA_MIRA) $(DLDFLAGS)



clean:
	rm *.o *.mod



#######################################################################################################
# Extra 

#CRAY Compiler
#(comf.f must be compiled with -O1)
#OPT = -O3
#OPT = -e D
#CFLAGS = $(F_COMP) $(OPT) -s real64 -N 132 $(DEBUG) $(EXTRA) $(PROFILE) -c
#LDFLAGS = $(F_COMP) $(EXTRA) -o

#PATHSCALE Compiler
#OPT = -O3
#EXTRA = -intrinsic=PGI
#CFLAGS = $(F_COMP) $(OPT) -r8 -extend-source $(DEBUG) $(EXTRA) $(PROFILE) -c
#LDFLAGS = $(F_COMP) -o

#GNU Compiler
#OPT = -O3
#DEBUG = -g -fbacktrace -fdump-core -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal
#DEBUG = -g -fbacktrace -fdump-core -fbounds-check
#EXTRA = -fmax-stack-var-size=100000
#PROFILE = -pg
#CFLAGS = $(F_COMP) $(OPT) -fdefault-real-8 -ffixed-line-length-none $(DEBUG) $(EXTRA) $(PROFILE) -c
#LDFLAGS = $(F_COMP) -o

#PGI Compiler
#DEBUG = -g -C
#OPT = -fastsse -O3 -Mvect -Minline=levels:10 -tp amd64
#CFLAGS  = $(F_COMP) -r8 $(OPT) -Mextend $(DEBUG) -c 
#LDFLAGS = $(F_COMP) -r8 $(DEBUG) -o



