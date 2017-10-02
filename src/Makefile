
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -check all -fpe0 -warn -traceback -debug extended 
#FCFLAGS = -O3 
#FCFLAGS = -O0 -g -ffast-math -funroll-loops -fbacktrace -fcheck=bounds -ffpe-trap=invalid,overflow,underflow,zero,denormal
FCFLAGS = -O3 -funroll-loops -march=native
#FCFLAGS = -O3 -xHost -opt-prefetch -unroll-aggressive -no-prec-div -fp-model fast=2

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

#INCLUDES = -I${FFTW_ROOT}/include/ -I${H5ROOT}/include/ -m64 -I$(MKLROOT)/include 
#INCLUDES = -I/usr/local/include/
#INCLUDES = -I/usr/local/include/ -I/usr/include/hdf5/serial/
INCLUDES = -I/usr/local/include -I/usr/local/hdf5/include/

OBJECTS = common.o sfunctions.o translations.o mie.o possu.o sparse.o singularity_subtraction.o singularity_subtraction_N.o integrals.o integration_points.o geometry.o sparse_mat.o io.o precorrection.o projection.o build_G.o gmres_module.o rhs.o field.o solver.o transformation_matrices.o T_matrix.o main.o

OBJECTS2 = common.o io.o sfunctions.o translations.o mie.o main_multiT.o

#Required libraries: Lapack, FFTW3, HDF5

#LIBS = -lm -L/usr/local/lib -lfftw3 -L/usr/lib -llapack -lhdf5_fortran -lhdf5 -lblas
LIBS = -lm -L/usr/lib -llapack -L/usr/local/hdf5/lib -lhdf5_fortran -lhdf5 -L/usr/local/lib -lfftw3

# Libs for csc taito gcc
#LIBS = -L${FFTW_ROOT}/lib/ -lfftw3_mpi -lfftw3 -L${H5ROOT}/lib/ -lhdf5_fortran -lhdf5 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

#intel
#LIBS = -L${FFTW_ROOT}/lib/ -lfftw3_mpi -lfftw3 -L${H5ROOT}/lib/ -lhdf5_fortran -lhdf5 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm



NAME = T_matrix
NAME2 = multi_T

all: ${NAME}

${NAME}: ${OBJECTS}
	 ${FC} ${FCFLAGS} ${OBJECTS} ${LIBS} -o ${NAME}

main.o: main.f90
	 ${FC} ${FCFLAGS} -c main.f90 -o main.o

%.o : %.f90
	 ${FC} ${FCFLAGS} ${INCLUDES} -c $*.f90 -o $*.o 

clean:
	 rm -f ${NAME}
	 rm -f ${OBJECTS}
	 rm -f *.mod
	 rm -f *~

veryclean:

	 rm -f ${NAME}
	 rm -f ${OBJECTS}
	 rm -f *.mod
	 rm -f *~
	 rm -f *.h5
	 rm -f *~
	 rm -f *~

all: ${NAME2}

${NAME2}: ${OBJECTS2}
	 ${FC} ${FCFLAGS} ${OBJECTS2} ${LIBS} -o ${NAME2}

main_multiT.o: main_multiT.f90
	 ${FC} ${FCFLAGS} -c main_multiT.f90 -o main_multiT.o

%.o : %.f90
	 ${FC} ${FCFLAGS} ${INCLUDES} -c $*.f90 -o $*.o 
