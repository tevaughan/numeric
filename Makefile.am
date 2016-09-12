
# Copyright 2016  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.

@CODE_COVERAGE_RULES@

SUBDIRS = src test docs .
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = .clang-format

dox :
	@(cd docs; make dox)

clean-local :
	@rm -frv numeric-*-coverage*

