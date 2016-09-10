
@CODE_COVERAGE_RULES@

SUBDIRS = src test .
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = .clang-format

dox :
	doxygen

all-local : dox

clean-local :
	@rm -frv numeric-*-coverage*

