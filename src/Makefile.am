
include_HEADERS =\
		 dimensions.hpp\
		 dimval.hpp\
		 integral.hpp\
		 interpolant.hpp\
		 units.hpp

BUILT_SOURCES =\
	       dimensions.hpp\
	       units.cpp\
	       units.hpp

CLEANFILES = $(BUILT_SOURCES)

lib_LTLIBRARIES = libnumeric.la

libnumeric_la_SOURCES = units.cpp

libnumeric_la_CPPFLAGS = $(CODE_COVERAGE_CPPFLAGS)

libnumeric_la_CXXFLAGS = $(CODE_COVERAGE_CXXFLAGS)

libnumeric_la_CFLAGS = $(CODE_COVERAGE_CFLAGS)

libnumeric_la_LDFLAGS = -version-info 0:0:0

libnumeric_la_LIBADD = $(CODE_COVERAGE_LIBS)

dimensions.hpp : dimensions.txt dimensions
	./dimensions

units.hpp units.cpp : units.txt units
	./units

clean-local :
	@rm -fv *.gcno

