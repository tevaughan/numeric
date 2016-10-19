
# Copyright 2016  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.

EXTRA_DIST = interp_1.gpi interp_2.gpi
LDADD = ../../src/libnumeric.la -lcln -lginac

noinst_PROGRAMS =\
 dense_table_1\
 integral_1\
 interp_1\
 interp_2\
 units_1\
 units_2

dense_table_1_SOURCES = dense_table_1.cpp
integral_1_SOURCES = integral_1.cpp
interp_1_SOURCES = interp_1.cpp
interp_2_SOURCES = interp_2.cpp
units_1_SOURCES = units_1.cpp
units_2_SOURCES = units_2.cpp

dox : $(EXTRA_DIST:.gpi=.png)

interp_1.png : interp_1.gpi interp_1.dat
	gnuplot $<

interp_2.png : interp_2.gpi interp_2.dat
	gnuplot $<

interp_1.dat : interp_1
	./$< > $@

interp_2.dat : interp_2
	./$< > $@

