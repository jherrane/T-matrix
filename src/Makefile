FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -O0 -g -ffast-math -funroll-loops -fbacktrace -fcheck=bounds -ffpe-trap=invalid,overflow,underflow,zero,denormal
FCFLAGS = -O3 -funroll-loops -march=native

H5ROOT = usr/lib/x86_64-linux-gnu/hdf5/serial
OPENBLAS_ROOT = usr
FFTW_ROOT = usr/local

INCLUDES = -I/${FFTW_ROOT}/include/ -I/${H5ROOT}/include/ -I/usr/include -m64  

OBJECTS = common.o sfunctions.o translations.o mie.o possu.o sparse.o singularity_subtraction.o singularity_subtraction_N.o integrals.o integration_points.o geometry.o sparse_mat.o io.o precorrection.o projection.o build_G.o gmres_module.o rhs.o field.o solver.o transformation_matrices.o T_matrix.o main.o

OBJECTS2 = common.o io.o sfunctions.o translations.o mie.o main_multiT.o

LIBS = -lm -L/$(OPENBLAS_ROOT)/lib -llapack -L/$(H5ROOT)/lib -lhdf5_fortran -lhdf5 -L/$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3

NAME = T_matrix

all: ${NAME}

${NAME}: ${OBJECTS}
	 ${FC} ${FCFLAGS} ${OBJECTS} ${LIBS} -o ${NAME}

main.o: main.f90
	 ${FC} ${FCFLAGS} -c main.f90 -o main.o

%.o : %.f90
	 ${FC} ${FCFLAGS} ${INCLUDES} -c $*.f90 -o $*.o 

clean:
	 rm -f ${NAME}
	 rm -f ${NAME2}
	 rm -f ${OBJECTS}
	 rm -f ${OBJECTS2}
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


