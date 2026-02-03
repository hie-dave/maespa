##
## Makefile for building a Maespa, could be nicer...
## Martin De Kauwe, 05/12/2013
##

PROG =	maespa.out
WG_LIB = libweathergen.so

SRCS=$(wildcard *.f90)
OBJS =	default_conditions.o switches.o getmet.o maindeclarations.o inout.o \
        maespa.o maestcom.o metcom.o physiol.o radn.o \
	    unstor.o utils.o watbal.o

# Objects for the standalone converter (avoid linking maespa.o)
WG_OBJS = weathergen.o getmet.o maestcom.o metcom.o switches.o radn.o utils.o

F90 = gfortran

RM ?= rm -f

INCLS = -I/usr/include
FFLAGS = -g -Wall -fbounds-check -finit-local-zero -Wuninitialized -ftrapv -ffree-form -ffree-line-length-none -O3 -fPIC

# CDI for CF-time and NetCDF I/O (future integration)
FFLAGS += $(shell pkg-config --cflags cdi_f2003 2>/dev/null)
LIBS =	$(shell pkg-config --libs cdi_f2003 2>/dev/null)

.PHONY: all clean run
.DEFAULT_GOAL := all

all: $(PROG) $(WG_LIB)

$(PROG): $(OBJS)
	$(F90) $(FFLAGS) $(OBJS) $(LIBS) -o $@

$(WG_LIB): $(WG_OBJS)
	$(F90) $(FFLAGS) $(WG_OBJS) $(LIBS) -shared -o $@

clean:
	$(RM) $(PROG) $(WG_LIB) $(OBJS) $(WG_OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(FFLAGS) $(INCLS) -c $< -o $@

default_conditions.o: switches.o
getmet.o: maestcom.o metcom.o switches.o
maindeclarations.o: maestcom.o
inout.o: maestcom.o switches.o
maespa.o: maestcom.o metcom.o switches.o maindeclarations.o
physiol.o: maestcom.o metcom.o
radn.o: maestcom.o
unstor.o: maestcom.o
utils.o: maestcom.o
watbal.o: maestcom.o metcom.o

# Ensure module build order for the standalone converter
weathergen.o: maestcom.o

run:
	./weathergen.jl \
		--file-tmin data/CumberlandPlain.nc --file-tmax data/CumberlandPlain.nc \
		--file-rs data/CumberlandPlain.nc --file-pr data/CumberlandPlain.nc \
		--file-ps data/CumberlandPlain.nc --file-ws data/CumberlandPlain.nc \
		--out-temp output/CumberlandPlain_temp.nc --out-rs output/CumberlandPlain_rs.nc \
		--out-pr output/CumberlandPlain_pr.nc --out-ps output/CumberlandPlain_ps.nc \
		--out-ws output/CumberlandPlain_wind.nc --out-vpd output/CumberlandPlain_vpd.nc \
		--name-tmax tasmax --name-tmin tasmin --show-progress --verbosity 3
