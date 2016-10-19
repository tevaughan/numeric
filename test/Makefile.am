
# Copyright 2016  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.

@CODE_COVERAGE_RULES@

EXTRA_DIST =\
 interpolant_test-badinput2.txt\
 interpolant_test-badinput.txt\
 interpolant_test.txt

TESTS = tests

check_PROGRAMS = tests

tests_SOURCES =\
 integral_test.cpp\
 interpolant_test.cpp\
 tests.cpp\
 units_test.cpp

# Need to move -lcln and -lginac to autoconf variable of some kind and add test
# for them to configure script.
tests_LDADD = ../src/libnumeric.la $(CODE_COVERAGE_LIBS)
tests_CPPFLAGS = $(CODE_COVERAGE_CPPFLAGS)
tests_CXXFLAGS = $(CODE_COVERAGE_CXXFLAGS)
tests_CFLAGS = $(CODE_COVERAGE_CFLAGS)

all-local : tests

clean-local :
	@rm -frv numeric-*-coverage
	@rm -fv *.gcda
	@rm -fv *.gcno

run-tests : tests
	./tests

