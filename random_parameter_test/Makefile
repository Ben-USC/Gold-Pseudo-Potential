#################################################
#    Makefile for construct pseudopotentials    #
#################################################
SHELL= /bin/sh

MAKEFILE  = Makefile

SDIR      = Sources

WDIR      = Atm
DDIR      = data

#################################################

PROGATM      = atm


$(PROGATM):	$(WDIR)/$(MAKEFILE)	$(WDIR)/$(PROGATM)


$(WDIR)/$(PROGATM):
	cd $(WDIR);  make $(PROGATM)


$(WDIR)/$(MAKEFILE):	$(SDIR)/$(MAKEFILE)
	@echo Execute: make compiler
	@make help
	@exit 1


#################################################

SOURCEFILE = $(SDIR)/$(MAKEFILE)
WORKFILE   = $(WDIR)/$(MAKEFILE)

ifort:	$(WDIR)	$(DDIR)	$(SOURCEFILE)
	sed "s/^#IFORT#//" $(SOURCEFILE) > $(WORKFILE)

ifortc:	$(WDIR)	$(DDIR)	$(SOURCEFILE)
	sed "s/^#IFORTC#//" $(SOURCEFILE) > $(WORKFILE)

frt:	$(WDIR)	$(DDIR)	$(SOURCEFILE)
	sed "s/^#FRT#//" $(SOURCEFILE) > $(WORKFILE)


$(WDIR):
	mkdir $(WDIR)

$(DDIR):
	mkdir $(DDIR)

$(SDIR)/$(MAKEFILE):
	@echo cannot find $(SDIR)/$(MAKEFILE)
	@exit 1

#################################################
help:
	@echo 'make ifort    : Intel(R) Fortran Compiler'
	@echo 'make frt      : Fujitsu Fortran Compiler'
