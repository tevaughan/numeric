
# numeric

## C++-11 library for numerical computation

The "numeric" project is a C++-11 library. It provides:
- [Types that allow an expression's dimension to be determined at compile
  time.](https://tevaughan.github.io/numeric/doxygen-html/dimvals.html)
- An adaptive-quadrature integrator.
```cpp
  #include <cmath>
  #include <iostream>
  #include <integral.hpp>
  #include <units.hpp>
  using namespace num;
  using prob_density = decltype(1.0 / length());
  length sigma = 3 * cm;
  prob_density gaussian(length x) {
     double const z = x / sigma;
     return exp(-0.5 * z * z) / sigma / sqrt(2.0 * M_PI);
  }
  int main() {
     double const prob = integral(gaussian, -3 * cm, 3 * cm);
     std::cout << prob << std::endl;  // 0.682689, about 68% for one sigma
     return 0;
  }
```
- An interpolant class with various, convenient constructors (from first two
  columns of ASCII file, from list of points, or from a continuous function
  that the constructed interpolant should approximate to a given tolerance).

  The following plot shows, as an example, the result of constructing a linear
  interpolant for a Gaussian.  One specifies a desired accuracy for the area of
  the interpolant.  Class num::interpolant's constructor, whose signature is
  the same as that for the global function integral() shown above, chooses the
  minimum number of points necessary to produce that accuracy.  Notice that
  points are dense where the curvature is large.

![Interpolant of Gaussian for Tolerance=1.0E-03 on Value of Integral](docs/examples/interp_1.png)

## Documentation

***[Documentation on the Github-Pages site for numeric.](https://tevaughan.github.io/numeric/doxygen-html)***

## To do:

 - Implement interpolant::invert() to return inverse function.
   - This will produce an error if the function be not invertible.
   - invert() will be useful if one want a quick, approximate inverse of a
     function that is invertible but not analytically invertible.
   - After an interpolant of sufficient accuracy is constructed, an interpolant
     for the inverse function is easy to construct!
 - Incorporate and adapt code from my linear regression project on github.

## License

Copyright 2016
Thomas E. Vaughan

Distributable under the terms of the GNU LGPL, Version 3 or later.

