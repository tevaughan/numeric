
EXTRA_DIST = interp_1.gpi
LDADD = ../libnumeric.la

noinst_PROGRAMS = units_1 interp_1

units_1_SOURCES = units_1.cpp
interp_1_SOURCES = interp_1.cpp

dox : interp_1.png

interp_1.png : interp_1.dat interp_1.gpi
	gnuplot interp_1.gpi

interp_1.dat : interp_1
	./interp_1 > $@

