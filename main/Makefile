# Compilador y flags
FC = mpifort
FFLAGS = -O3
#FFLAGS = -g -O3  -fbounds-check -fcheck=all -fbacktrace 
LNK = mpifort

# Rutas
LIB = /home/flavio/usr/lib
INC = /home/flavio/include

# Objetos
OBJS = main.o initial_conditions.o evolve.o save_gnuplot.o

MODS = \
../mpi/mpi_lib.o \
../hdf5/hdf5_lib.o \
../tools/finite_differences.o ../tools/integral.o ../tools/ODE.o ../tools/strings_lib.o \
../tools/FFT_sph.o  ../tools/fft_lib.o\
../IO/save_data.o ../IO/save_data_lib.o ../IO/read_data.o \
../MoL/mesh_refinement_lib.o ../MoL/geometry.o ../MoL/lines_method_lib.o ../MoL/diagnostic_lib.o ../tools/random.o \
global_numbers.o 

# Dependencias
$(OBJS): $(MODS)

main: $(OBJS) $(MODS)
	$(LNK) $(FFLAGS) -L$(LIB) -I$(INC) -o main $(OBJS) $(MODS) -lhdf5 -lhdf5_fortran -lgsl -lgslcblas -lm
	@ mkdir -p exe
	@ mv main exe

.PHONY: clean

clean:
	-rm -f *.o exe/main ../hdf5/*.o ../hydrodinamics/*.o  ../IO/*.o ../MoL/*.o ../mpi/*.o ../Poisson/*.o ../tools/*.o *.mod *.o 

%.o : %.f90
	$(FC) -c $(FFLAGS) -I$(INC) $< -o $@
