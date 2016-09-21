
# Copyright 2016  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.

EXTRA_DIST = interp_1.gpi interp_2.gpi
LDADD = ../../src/libnumeric.la

noinst_PROGRAMS =\
 units_1\
 units_2\
 integral_1\
 interp_1\
 interp_2

units_1_SOURCES = units_1.cpp
units_2_SOURCES = units_2.cpp
integral_1_SOURCES = integral_1.cpp
interp_1_SOURCES = interp_1.cpp
interp_2_SOURCES = interp_2.cpp

dox : $(EXTRA_DIST:.gpi=.png)

interp_%.png : interp_%.gpi interp_%.dat
	gnuplot $<

interp_%.dat : interp_%
	./$< > $@

