# Makefile for Radex

# This file is part of the RADEX software package
# to calculate molecular excitation and radiative
# transfer in a homogeneous medium.
#
# Documentation for the program is posted at
# https://sron.rug.nl/~vdtak/radex/index.shtml  
#
# Although this program has been thoroughly tested, the
# authors do not claim that it is free of errors and
# gives correct results in all situations.
#
# Publications using this program should make a reference
# to our paper: A&A 468, 627 (2007).


# it's modified by bakl 2024/11/28


# Do you want to compile with debug options? 1=Yes, 0=No
IS_DEBUG = 0

#sys := $(shell uname)

# SYSTYPE ?= "ifort"
# SYSTYPE := "lf"
SYSTYPE ?= gfortran

MD := mkdir
MV := mv
DEL = rm -f

########################################################################
# Paths

ifeq ($(RADEXHOME),)
#  PWD       := $(abspath $(shell pwd))
 PWD       = .
 RADEXHOME   := $(PWD)/..
#  STLHOME := $(abspath $(shell pwd))
 export RADEXHOME
endif

SRCDIR    := $(RADEXHOME)/src/
BINDIR    := $(RADEXHOME)/bin/

vpath %.f   $(SRCDIR)
vpath %.f90 $(SRCDIR)
vpath %.o   $(SRCDIR)


VPATH := $(SRCDIR) $(BINDIR)


########################################
## For gfortran fortran
ifeq ($(SYSTYPE),gfortran)
  FC = gfortran
  ifeq ($(IS_DEBUG),1)
    FFLAGS_F77 += -g -O0 -fbounds-check -Wuninitialized -Wall
    FFLAGS_F77 += -ffpe-trap=invalid,zero,overflow -finit-real=nan
    # FFLAGS_F77 += -g  #  debugging
    # FFLAGS_F77 += -ffpe-trap=invalid,zero,overflow  #  debugging
    # FFLAGS_F77 += -finit-real=nan  #  debugging
    FFLAGS_F77 += -Wall  #  debugging
    FFLAGS_F77 += -fcheck=all  #  debugging
    # FFLAGS_F77 += -Warray-temporaries
    FFLAGS_F77 += -Wextra   -Wconversion #  debugging
    FFLAGS_F77 += -pedantic #  debugging
    FFLAGS_F77 += -fbacktrace
    FFLAGS_F77 += -frange-check
    FFLAGS_F77 += -fdump-core
    FFLAGS_F77 += -Wno-do-subscript  -Wno-surprising #  tmp
    FFLAGS_F77 += -Wno-unused-dummy-argument #  tmp
    FFLAGS_F77 += -Wno-maybe-uninitialized # gcc bug 
  else
    FFLAGS_F77 += -O2  # optimize 
  endif
  # FFLAGS_F77 =  -fbounds-check   #  -fno-range-check 
  FFLAGS_F77 += -ffixed-line-length-none
  FFLAGS_F77 += -ffree-line-length-none #-ffixed-line-length-132 
  # FFLAGS_F77 += -fdefault-real-8 
  FFLAGS_F77 += -cpp   # preprocessor
  # FFLAGS_F77 += -O2   # optimize 
  # FFLAGS_F77 += -O0   # valgrind --track-origins=yes --keep-stacktraces=alloc-and-fre --leak-check=full 
  
  FFLAGS_F77 += -finit-local-zero
  # FFLAGS_F77 += -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan  #  debugging
  # FFLAGS_F77 += -Wno-function-elimination   # tmp https://gcc.gnu.org/onlinedocs/gfortran/Evaluation-of-logical-expressions.html
	FFLAGS_F77 += -ftrapping-math
  FFLAGS_F77 += -D_GFORTRAN
  # FFLAGS_F77 += -fopenmp  

  # FFLAGS_F77 += -fmax-stack-var-size=4500000
  
  # FFLAGS_F77 += -std=f2008   #  debugging
  # FFLAGS := -ffree-form
  # FFLAGS := $(FFLAGS_F77) $(FFLAGS)
  FFLAGS := $(FFLAGS_F77) $(FFLAGS) -ffree-form

  # FFLAGS_F77 += -frecursive  #  https://gcc.gnu.org/onlinedocs/gcc-4.6.3/gfortran/Code-Gen-Options.html
  # FFLAGS_F77 += -Wno-unused      # -Wno-unused-parameter -Wunused-dummy-argument
  FFLAGS_F77 +=-Wno-unused-label
#   FFLAGS_TRF := $(FFLAGS_F77) -fno-automatic  #  https://gcc.gnu.org/onlinedocs/gcc-4.6.3/gfortran/Code-Gen-Options.html

  LDFLAGS = 

endif
########################################

# Option 3:
#FC      = ifort
#FFLAGS += -O2
# Option 4:
#FC = g95
#FFLAGS += -Wall -O2
# there may be other options of course ...
FILES =  cla.f90 mRadexInc.f90 \
         io.f90 readdata.f90 matrix.f90 background.f90 slatec.f90 main.f90 
TEMP := $(FILES:.f=.o)
OBJECTS := $(patsubst %.f90,%.o, $(TEMP))

# OBJECTS = mRadexInc.o main.o io.o readdata.o matrix.o background.o slatec.o
EXEC    = radex

.IGNORE:
LDFLAGS += -o

# %.o: %.f radex.inc Makefile
# 	$(COMPILE.F) $*.f -o $@
### Dependences

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@  $(INCL_DIR)
  
%.o: %.f
	$(FC) $(FFLAGS_F77) -c $< -o  $@  $(INCL_DIR)

%.o: %.mod ## if make wants m2c thinking that mod is from Modula language

%.mod : %.o
	@if [! -f $@ ]; then \
	@rm $< \
	$(MAKE) $< \
	fi

# %.o: %.f90 Makefile
# 	$(COMPILE.F) $*.f -o $@

.PHONY:

all: help 
	@echo files: $(FILES) 
	@echo OBJECTS: $(OBJECTS)  

radexinstall: radex install clean  ## Build radex, install radex, clean *-,mod-files
	@ls -l $(BINDIR)/$(EXEC)

radex: $(OBJECTS) $(EXEC)  ## Build radex
	$(FC) $(LDFLAGS) $(EXEC)  $(OBJECTS) 
#	strip $@
	@ls -l $(EXEC)


install:   ## Install radex
	install -m 755 -p -s $(EXEC) $(BINDIR)


help: ## Show this dynamic help
	$(info Available targets)
	@awk '/^[a-zA-Z\-_0-9]+:/ {                                   \
      nb = sub( /^## /, "", helpMsg );                             \
      if(nb == 0) {                                                \
        helpMsg = $$0;                                             \
        nb = sub( /^[^:]*:.* ## /, "", helpMsg );                  \
      }                                                            \
      if (nb)                                                      \
        printf "\033[1;31m%-" width "s\033[0m %s\n", $$1, helpMsg; \
                                                                   \
      helpMsg = $$0 }'                                             \
    width=$$(grep -o '^[a-zA-Z_0-9]\+:' $(MAKEFILE_LIST) |  wc -L)  \
    $(MAKEFILE_LIST)


.PHONY:

clean:  ## Clean obj-files, etc 
	rm *.o *.mod
#	 - rm $(EXEC)
