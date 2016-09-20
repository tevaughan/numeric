
# Introduction  {#mainpage}

The "numeric" project is a C++-11 library that provides physical units, an
adaptive quadrature integrator, and various ways to represent a functional
relation numerically.

## Physically Dimensioned Quantities

The implementation of physical units with compile-time checking for dimensions
allows one to do things like this:

```cpp
#include <iostream>
#include <units.hpp>

using namespace num;
using namespace std;

int main()
{
   // 'ft' stands for one foot.
   length ell = 2.1 * ft;                    // length of pendulum
   // 'gn' stands for one standard gravity.
   frequency f = sqrt(gn / ell);             // frequency of small oscillations
   cout << f << endl;                        // [3.9142 s^-1]
   return 0;
}
```

The compiler would produce an error if, say, a distance were assigned to a
speed.

[Details on how to use dimensions and units.](\ref dimvals)

## Integration

An adaptive-quadrature, trapezoid-rule, global integration function,
num::integral allows the definite integral of any continuous function of one
variable to be computed to a specified tolerance.
- The template function num::integral is designed to be able to interoperate
  with units.
- If dimensioned quantities be not used, then one need include only
  [integral.hpp](\ref integral.hpp) and need not link against `libnumeric.so`.

[Details on integration.](\ref integrals)

## Functional Relations

A function written in a computer language like C or C++ can be used to model
many a mathematical function.  Sometimes, though, a C++ function, when it is
called, can take more time to return than is practical. There are various ways
to trade memory in order to gain faster access to an approximation of the value
that the function would return for a given input.

### Table

If memory be abundant and fast, an approximant may be built from a table of
pairs of values, one from the function's domain and the other from the
function's range.  Then, for a given argument to the function, the approximant
returns the range's element corresponding to the nearest tabulated element in
the domain or to an interpolant over nearby tabulated elements.

#### Dense Table

The fastest version of this approach can be applied when fast memory is most
abundant and when there is no constraint on the time required to set up the
table before it is used.  This approach requires that every entry in the table
be separated from the next by the same distance in the domain.  Such an
approach allows the approximate element of the range to be returned in constant
time, independent of the number of entries in the table.  The num::dense\_table
class is an implementation of this approach.

#### Sparse Table

A slower version of this approach requires less memory.  In this approach,
there is no requirement of uniformity imposed on the distance in the domain
between one entry in the table and the next.  Such an approach allows the
approximate element of the range to be returned in a time proportional to the
logarithm of the number of entries in the table. The num::sparse\_table class
is an implementation of this approach.

### Interpolation

Regardless of whether a dense table or a sparse table be used to store values
of a function, the question of how to make use of the stored values remains.
The simplest

TBS

When the interpolant is constructed from a continuous function, the same
algorithm employed in num::integral is used, but the intermediate results of
the adaptive subdivision are stored for use as points for interpolation.  So
the interpolant is intelligently constructed to have a high density of points
only where necessary to ensure an accurate representation of the function and
its integral. Notice how the points are distributed in the following
interpolant for a Gaussian:

![Interpolant of Gaussian for Tolerance=1.0E-03 on Value of Integral](interp_1.png)

The num::interpolant class is designed to interoperate with units.

However, if dimensioned quantities be not used, then one need only include
[interpolant.hpp](\ref interpolant.hpp) and need not link against
`libnumeric.so`.

[Details of linear interpolation.](\ref linterp)

## Installation

Either

 - Grab the latest 'tar.gz' file from [the releases
   page](https://github.com/tevaughan/numeric/releases),

 - or clone from [github](https://github.com/tevaughan/numeric).

Build and install:
```
./configure
make
sudo make install
```

All classes and functions inter-operate!

## License

Copyright 2016
Thomas E. Vaughan

Distributable under the terms of the GNU LGPL, Version 3 or later.

