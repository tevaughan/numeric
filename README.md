
# numeric

## C++-11 library for numerical computation

The C++-11 numeric library provides physical units for compile-time checking of
physical modeling, an adaptive quadrature integrator, and the ability to
instantiate a linear interpolant in various, convenient ways.

The following plot shows, as an example of the result of constructing a linear
interpolant for a Gaussian. One specifies a desired accuracy for the area, and
the constructor chooses the minimum number of points necessary to produce that
accuracy. Notice that points are dense where the curvature is large.

![Interpolant of Gaussian for Tolerance=1.0E-03 on Value of Integral](src/examples/interp_1.png)

## Documentation

***[Documentation can be viewed on the Github-Pages site for numeric.](https://tevaughan.github.io/numeric/doxygen-html)***

## To do:

 - Incorporate and adapt code from my linear regression project on github.

## License

Copyright 2016
Thomas E. Vaughan

Distributable under the terms of the GNU LGPL, Version 3 or later.

