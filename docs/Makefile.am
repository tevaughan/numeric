
# Copyright 2016-2017  Thomas E. Vaughan
#
# As part of the numeric library, this file may be distributed under the terms
# of the GNU LPGL, Version 3 or later.

SUBDIRS = examples .

EXTRA_DIST = Doxyfile

dox :
	(cd examples; make dox)
	doxygen

clean-local :
	@rm -frv doxygen-html/*

