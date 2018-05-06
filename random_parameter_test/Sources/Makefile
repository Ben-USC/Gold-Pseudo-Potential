#################################################
#    Makefile for construct pseudopotentials    #
#################################################
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES:	.F90 .o .l

EXTF90  = f90

#IFORT#    LINKER  = ifort
#IFORT#    FFLAGS  = -c
#IFORTC#    LINKER  = ifort
#IFORTC#    FFLAGS  = -c -check all -traceback

#FRT#    LINKER  = frt
#FRT#    FFLAGS  = -c -Am -Eg





FC	= $(LINKER)

SDIR	= ../Sources

PROGRAM	= atm

OBJS	= \
          ftmain.o \
          funcs.o \
          input.o \
          vxc.o \
          ae.o \
          vloc.o \
          pp.o \
          ecut.o \
          trans.o

.F90.o:
	$(FC) $(FFLAGS) $(CPPDEFS) $*.F90

$(PROGRAM):	$(OBJS)
	@echo "Loading $(PROGRAM) ... "
	$(LINKER) $(OBJS) -o $(PROGRAM)
	mv $(PROGRAM) ..
	@echo "done"

###############
# maintenance #
###############
clean:
	@rm -f $(OBJS)
	@rm -f *.o *.mod *.MOD

allclean:
	@rm -f $(OBJS)
	@rm -f *.o *.mod *.MOD
	@rm -f *.F90 *.h *.i


input.o:   ftmain.o
ae.o:      ftmain.o
vloc.o:    ftmain.o
pp.o:      ftmain.o
ecut.o:    ftmain.o
trans.o:   ftmain.o


ftmain.F90:   $(SDIR)/ftmain.f90   ; cp $(SDIR)/ftmain.f90 ./ftmain.F90
funcs.F90:    $(SDIR)/funcs.f90    ; cp $(SDIR)/funcs.f90  ./funcs.F90
input.F90:    $(SDIR)/input.f90    ; cp $(SDIR)/input.f90  ./input.F90
vxc.F90:      $(SDIR)/vxc.f90      ; cp $(SDIR)/vxc.f90    ./vxc.F90
ae.F90:       $(SDIR)/ae.f90       ; cp $(SDIR)/ae.f90     ./ae.F90
vloc.F90:     $(SDIR)/vloc.f90     ; cp $(SDIR)/vloc.f90   ./vloc.F90
pp.F90:       $(SDIR)/pp.f90       ; cp $(SDIR)/pp.f90     ./pp.F90
ecut.F90:     $(SDIR)/ecut.f90     ; cp $(SDIR)/ecut.f90   ./ecut.F90
trans.F90:    $(SDIR)/trans.f90    ; cp $(SDIR)/trans.f90  ./trans.F90

