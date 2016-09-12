
@CODE_COVERAGE_RULES@

SUBDIRS = src test docs .
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = .clang-format

dox :
	@(cd src; make dox)
	@(cd docs; make dox)

clean-local :
	@rm -frv numeric-*-coverage*

