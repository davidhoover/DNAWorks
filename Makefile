# compiler
FC = gfortran

# compile flags
#FCFLAGS = -g -fbounds-check -O2 -static-libgcc -static
FCFLAGS = -g -fbounds-check -O2

# link flags
FLFLAGS = -g

# program name
PROGRAM = dnaworks

# required objects
objects = dnaworks.o dnaworks_data.o dnaworks_test.o \
	control_func.o email_func.o encoding.o input.o misc_func.o \
	mutate.o output.o overlaps.o scores.o str_func.o time_func.o

# required modules
modules = dnaworks_data.mod dnaworks_test.mod

# the main linking step
$(PROGRAM): $(objects)
	$(FC) $(FCFLAGS) -o $(PROGRAM) $(objects)

# specific requirements for each object
$(objects): $(modules)

# compile recipe for modules
%.mod: %.f90
	$(FC) $(FLFLAGS) -c $<

# compile recipe for objects
%.o: %.f90
	$(FC) $(FLFLAGS) -c $<

# extra rules
.PHONY: clean
clean:
	rm -f $(objects) $(modules) $(PROGRAM)
