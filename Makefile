.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod

BINDIR = ./
FC = mpifort

FFLAGS = -O3 -ffree-line-length-none # -g -Wunused-variable -fbounds-check
FINC   =
LDFLAGS=
LD = $(FC) 

# Object files
OBJS = main.o prec_def.o sbp_6_order.o mms.o 

.PHONY:clean

$(BINDIR)/Maxwell_TDSE.x: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)

.f.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.mod:
	$(FC) -c $(FFLAGS) $(FINC) $<

# Dependencies betwenn modules
main.f90: prec_def.mod mms.mod
mms.f90: prec_def.mod
sbp_6_order.f90: prec_def.mod

clean:
	rm -f *.o *.a *.mod *.txt main.x