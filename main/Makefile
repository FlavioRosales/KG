#FC = h5pfc
#FC = mpif90
FC = mpifort
#FC=gfortran
#FC=/opt/intel/bin/ifort

FFLAGS = -g -O3 # -fbounds-check -fcheck=all -fbacktrace
#FFLAGS = -g

#LNK = mpif90
LNK = mpifort
#LNK = h5pfc
#LNK=gfortran
#LNK=/opt/intel/bin/ifort

LIB=home/flavio/usr/lib
INC=home/flavio/include
#LIB=usr/lib/x86_64-linux-gnu/hdf5/openmpi
#INC=usr/include/hdf5/openmpi

OBJS = main.o initial_conditions.o evolve.o save_gnuplot.o

MODS = \
../mpi/mpi_lib.o \
../hdf5/hdf5_lib.o \
../tools/finite_differences.o ../tools/integral.o ../tools/ODE.o ../tools/strings_lib.o \
../IO/save_data.o ../IO/save_data_lib.o ../IO/read_data.o \
../MoL/mesh_refinement_lib.o ../MoL/geometry.o ../MoL/lines_method_lib.o \
global_numbers.o 
#../hydrodinamics/geometry_lib.o ../hydrodinamics/hydro_base_lib.o ../hydrodinamics/hydro_lib.o ../hydrodinamics/perfect_fluid_lib.o \


$(OBJS):	$(MODS)

main:	$(OBJS) $(MODS)
		$(LNK) $(FFLAGS) -L/$(LIB) -o main $(OBJS) $(MODS) -lhdf5 -lhdf5_fortran
	@ mkdir -p exe
	@ mv main exe
#	@ cp input.par xxx


.PHONY:	clean

clean:
	-rm -f *.o exe/main ../hdf5/*.o ../hydrodinamics/*.o ../schrodinger/*.o ../IO/*.o ../MoL/*.o ../mpi/*.o ../Poisson/*.o ../tools/*.o *.mod *.o 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@ -I/$(INC)
