INC=-I /usr/lib/x86_64-linux-gnu/openmpi/include
LIB=  #-L/usr/lib/x86_64-linux-gnu -lplplotd -lplplotf95d
INC=
FC=gfortran
FCPARA=mpif90
FCOPT= -fimplicit-none -finit-local-zero  -fopenmp -finit-local-zero  -ffpe-trap=zero,overflow,underflow  -cpp -DDEBUG # -Wunused -Wall -Wargument-mismatch -Wimplicit-interface

#FCOPT=-O2  -stand f03    -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone

#FCOPT=-g  -cpp -DDEBUG
# make FCOPT="-g  -cpp -DDEBUG"

SRC=analysis.f90

OBJ= global.o Atom.o Cell.o Bond.o Configuration.o Element.o   Machine.o IO.o Molecule.o Univers.o tools.o scripts.o
#all: $(OBJ) analysis.f90
#	$(FC) $(FCOPT)  $(OBJ) analysis.f90  -o analysis.x  -lblas -llapack	

#all: $(OBJ) generate_configuration_v1.1.f90
#	$(FC) $(FCOPT)   $(OBJ) generate_configuration_v1.1.f90  -o generate_configuration.x  -lblas -llapack
all: $(OBJ) new.f90
	$(FC) $(FCOPT)   $(OBJ) new.f90  -o generate_configuration.x  -lblas -llapack
Atom.o: Atom.f90 global.o
	$(FC) $(FCOPT) $< -c
Bond.o: Bond.f90 global.o
	$(FC) $(FCOPT) $< -c
Cell.o: Cell.f90  Bond.o global.o tools.o Molecule.o
	$(FC) $(FCOPT) $< -c
Configuration.o: Configuration.f90 global.o Molecule.o
	$(FC) $(FCOPT) $< -c
Element.o: Element.f90 global.o
	$(FC) $(FCOPT) $< -c
global.o: global.f90
	$(FC) $(FCOPT) $< -c
IO.o: IO.f90 Element.o global.o
	$(FC) $(FCOPT) $< -c
Machine.o: Machine.f90 global.o
	$(FC) $(FCOPT) $< -c
Molecule.o: Molecule.f90 Atom.o Element.o global.o
	$(FC) $(FCOPT) $< -c
Univers.o: Univers.f90 global.o 
	$(FC) $(FCOPT) $< -c
tools.o: tools.f90 
	$(FC) $(FCOPT) $< -c
scripts.o: scripts.f90 global.o Molecule.o Univers.o Cell.o IO.o Atom.o
	$(FC) $(FCOPT) $< -c



#serial:  $(SRC)
#	$(FC) $(FCOPT)  $(SRC) -o Hbinitio.x  -lblas -llapack

#analysis.o: analysis.f90#
#	$(FC) $(FCOPT) $< -c
clean:
	rm *.mod *.o
#param.mod: param.f90 global.o
#	$(FC) $(FCOPT) $< -c



#IO.mod: IO.f90 global.o tools.o
#	$(FC) $(FCOPT) $< -c



#	$(FCPARA) $(FCOPT) Hbinitio.f90 -o Hbinitio.x $(INC) -lblas -llapack -lmpi
#pp: ppHbinitio.f90#
#	$(FC) ppHbinitio.f90 -o ppHbinitio.x $(LIB)
#parallel: mpi_Hbinitio.f90
#	$(FCPARA) $(FCOPT) mpi_Hbinitio.f90 -o mpi_Hbinitio.x $(INC) -lblas -llapack -lmpi
#dbg: dbg.f90
#	$(FCPARA) dbg.f90 -o dbg.x  $(INC) -lblas -llapack -lmpi
