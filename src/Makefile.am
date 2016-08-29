
include_HEADERS =\
		 dimensions.hpp\
		 dimval.hpp\
		 integral.hpp\
		 interpolant.hpp

BUILT_SOURCES =\
	       dimensions.hpp

CLEANFILES =\
	    dimensions.hpp

lib_LTLIBRARIES = libnumeric.la
libnumeric_la_SOURCES = dimval.cpp
libnumeric_la_LDFLAGS = -version-info 0:0:0

dimensions.hpp : dimensions.txt dimensions
	./dimensions

