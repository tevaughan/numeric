
# Copyright 2016-2017  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.


# A compiled library is built for numeric.  At the moment, this library exists
# only to support units for dimval.

# Libtool's versioning system:
#
# current
#
#     The most recent interface number that this library implements.
#
# revision
#
#     The implementation number of the current interface.
#
# age
#
#     The difference between the newest and oldest interfaces that this library
#     implements. In other words, the library implements all the interface
#     numbers in the range from number current - age to current.

current  = 4
revision = 1
age      = 0

BUILT_SOURCES = dimensions.hpp units.cpp units.hpp
CLEANFILES = $(BUILT_SOURCES)
EXTRA_DIST = *.pl *.txt *.md

pkginclude_HEADERS =\
 dense-table.hpp\
 dim-exps.hpp\
 dimval.hpp\
 ilist.hpp\
 integral.hpp\
 integral-stats.hpp\
 interpolant.hpp\
 interval.hpp\
 rk.hpp\
 sparse-table.hpp\
 util.hpp

nodist_pkginclude_HEADERS = dimensions.hpp units.hpp

lib_LTLIBRARIES = libnumeric.la
nodist_libnumeric_la_SOURCES = dimval.cpp sparse-table.cpp units.cpp
libnumeric_la_CPPFLAGS = $(CODE_COVERAGE_CPPFLAGS)
libnumeric_la_CXXFLAGS = $(CODE_COVERAGE_CXXFLAGS)
libnumeric_la_CFLAGS = $(CODE_COVERAGE_CFLAGS)
libnumeric_la_LDFLAGS = -version-info $(current):$(revision):$(age)
libnumeric_la_LIBADD = $(CODE_COVERAGE_LIBS)

dimensions.hpp : dimensions.txt dimensions.pl
	./dimensions.pl

units.hpp units.cpp : units.txt units.pl
	./units.pl

clean-local :
	@rm -fv *.gcno

