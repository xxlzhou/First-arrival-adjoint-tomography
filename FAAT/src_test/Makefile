include ./Makefile.inc

EXEDIR=../bin
	
.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) $(INC) -c $*.f90

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) -c $*.c

#.SUFFIXES: .o .f90
#.f90.o:
#       $(FC) -c $(FFLAGS) $(INC) $<

%.o: %.f90
	$(FC) $(OPTF) $(INC) -c -o $@ $<

all: compile

compile: FAAT

clean:
	rm *.o *.mod

SUB = FAATd.o modforward_rec.o subdata_rec.o subresiduals_rec.o subgradient.o submisfit_rec.o modadjoint_rec.o submodels.o scalegradient.o
FTDET2d : $(SUB)
	$(F90) $(OPTF) $(SUB) $(LIB) -o $(EXEDIR)/FAAT
