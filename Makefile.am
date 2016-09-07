
@CODE_COVERAGE_RULES@

SUBDIRS = src test
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = .clang-format

doc :
	doxygen

clean-local :
	@rm -frv doc
	@rm -frv numeric-*-coverage*

