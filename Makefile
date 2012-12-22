cc = gfortran

ode = dvode_f90_m.o
newdiff = newdiff.o
obj = $(ode) $(newdiff)

default = newdiff

clean = *.o newdiff *.mod 
cleand = *.o newdiff *.mod *.dat

newdiff: $(obj)
	$(cc) -o newdiff $(obj)

%.o: %.f90
	$(cc) -c $*.f90

$(newdiff): $(ode) 

clean:
	rm -rf $(clean)

cleand:
	rm -rf $(cleand)
