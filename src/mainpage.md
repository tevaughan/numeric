
# Introduction  {#mainpage}

The C++-11 numeric library provides physical units, an adaptive quadrature
integrator, and the ability to instantiate a linear interpolant in various,
convenient ways.

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
variable to be computed a specified tolerance.

The global num::integral is designed to be able to interoperate with units.

However, if dimensioned quantities be not used, then one need only include
[integral.hpp](\ref integral.hpp) and need not link against `libnumeric.so`.

[Details of integration.](\ref integrals)

## Linear Interpolation

There are two basic ways to construct a linear interpolant:

1. from a list of points and
2. from a continuous function and an interval on its domain.

Once a linear interpolant is constructed, it can be integrated efficiently and
exactly over any interval. An interpolant can also be multiplied by another
interpolant or by a scalar without any loss of information.

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

