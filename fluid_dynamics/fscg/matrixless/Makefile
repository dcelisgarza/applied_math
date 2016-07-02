# Defining variables.
target1 = test.out
target2 = fluid.out
objects1 = nrtype.o fluid.o test_fluid.o
objects2 = nrtype.o fluid.o fluid_main.o
f08comp = gfortran-6
FCFLAGS = -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
#FCFLAGS = -Ofast -funsafe-loop-optimizations -fmerge-all-constants
# Makefile
default: $(target2) $(target1)
$(target1): $(objects1)
	$(f08comp) $(objects1) $(FCFLAGS) -o $@
$(target2): $(objects2)
	$(f08comp) $(objects2) $(FCFLAGS) -o $@
nrtype.o: ../../../nrtype.f08
	$(f08comp) $(FCFLAGS) -c ../../../nrtype.f08
fluid.o: fluid.f08
	$(f08comp) $(FCFLAGS) -c fluid.f08
test_fluid.o: test_fluid.f08
	$(f08comp) $(FCFLAGS) -c test_fluid.f08
fluid_main.o: fluid_main.f08
	$(f08comp) $(FCFLAGS) -c fluid_main.f08
# Cleaning everything
clean:
	rm -f *.o *.mod *.dat *.gnu *.gif
clobber: clean
	rm -f $(target)
.PHONY: default clean cobbler
