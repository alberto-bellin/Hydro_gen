

F90= ifort 
F90FLAGS = -O3


VPATH=.

# Define all Fortran object files which make up hydrogen2D

OBJECTS = \
	vars_common_hydrogen.o \
	SEED_FROM_URANDOM.o \
	COMB.o \
	BESSIK.o \
	FUNC_COVY.o \
	COVARIAN.o \
	covspaz_2D.o \
	COEFCY.o \
	COEFL2.o \
	coeff_2D.o \
	normal.o \
	LINPACK.o \
	BLAS.o \
	hydrogen2D_input.o \
	main_hydrogen.o \
	gen_field2D.o \
	
        
 

hydrogen2D: $(OBJECTS) 
	-$(F90) $(F90FLAGS) -o hydrogen2D $(OBJECTS) 

# Object codes

.f.o:
	$(F90) $(F90FLAGS) -c $<


vars_common_hydrogen.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/vars_common_hydrogen.f90
SEED_FROM_URANDOM.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/SEED_FROM_URANDOM.f90
COMB.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/COMB.f90
BESSIK.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/BESSIK.f90
FUNC_COVY.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/FUNC_COVY.f90
COVARIAN.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/COVARIAN.f90
covspaz_2D.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/covspaz_2D.f90
COEFCY.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/COEFCY.f90
COEFL2.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/COEFL2.f90
coeff_2D.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/coeff_2D.f90
normal.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/normal.f90
LINPACK.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/LINPACK.f90
BLAS.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/BLAS.f90
hydrogen2D_input.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/hydrogen2D_input.f90
main_hydrogen.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/main_hydrogen.f90
gen_field2D.o:
	$(F90) $(F90FLAGS) -c   $(VPATH)/gen_field2D.f90





#clean:
#	rm $(VPATH)/*.o
#	rm $(VPATH)/*.mod

#  end

