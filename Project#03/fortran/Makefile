# Flags for GNU
FFLAGS = -g -fcheck=all
FC = gfortran

# Suffixes
.SUFFIXES:
	.SUFFIXES: .f90 .F90 .o .mod
		.f90.o:
			$(FC) $(FFLAGS) -c $<
		.f90.mod:
			$(FC) $(FFLAGS) -c $<
		.F90.o:
			$(FC) $(FFLAGS) -c $<
		.F90.mod:
			$(FC) $(FFLAGS) -c $<

# Declare all files
OBJ = const.F90 system.f90 read_in_integrals.f90 main.F90

# Rules for making
els.x: $(OBJ)
	 $(FC) $(FFLAGS) -o $@ $(OBJ)

const.o const.mod: const.F90
	$(FC) $(FFLAGS) -c const.F90

# Clean
clean:
	rm *.o *.mod els.x
