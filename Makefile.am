
@CODE_COVERAGE_RULES@

SUBDIRS = src test .
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = .clang-format

dox :
	doxygen

clean-local :
	@rm -frv numeric-*-coverage*

