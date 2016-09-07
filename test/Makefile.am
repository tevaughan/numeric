
@CODE_COVERAGE_RULES@

EXTRA_DIST = interpolant_test.txt

TESTS = tests

check_PROGRAMS = tests

tests_SOURCES =\
	       integral_test.cpp\
	       interpolant_test.cpp\
	       units_test.cpp\
	       tests.cpp

tests_LDADD = ../src/libnumeric.la $(CODE_COVERAGE_LIBS)

tests_CPPFLAGS = $(CODE_COVERAGE_CPPFLAGS)

tests_CXXFLAGS = $(CODE_COVERAGE_CXXFLAGS)

tests_CFLAGS = $(CODE_COVERAGE_CFLAGS)

clean-local :
	@rm -frv numeric-*-coverage
	@rm -fv *.gcda
	@rm -fv *.gcno

