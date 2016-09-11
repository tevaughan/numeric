
# Introduction  {#mainpage}

The numeric library provides

 - a compiled library (linked with `-lnumeric`) for dimensioned quantities and

 - a header-only library for
   - numerical integration and
   - interpolation.

## Physically dimensioned quantities.

Physically dimensioned quantities allow one to do things like this:
```cpp
#include <units.hpp>
using namespace num;
acceleration const g = 9.81 * m / pow<2>(s); // acceleration of gravity
length ell = 1 * m;                          // length of pendulum
frequency f = sqrt(g / ell);                 // frequency of small oscillations
```
The compiler will tell you if you make a mistake dimensionally.

[Details of how to use dimensions and units.](\ref dimvals)

## Linear interpolation

A linear interpolant can be initialized directly from a file with two
space-delimited columns.  Then it can be scaled, multiplied by another
interpolant with no loss of information, integrated, etc.

```cpp
#include <iostream>
#include <interpolant.hpp>
#include <units.hpp>

using namespace num;
using namespace std;

// Read in first two columns of space-delimited file.
// First column is interpreted as a double-precision length in nanometers.
// Second column is interpreted as a double-precision time in seconds.
// Explicit namespace scoping for num::time might be required to avoid
// conflict with C library's time() function.
interpolant<length, num::time> i("file.txt", nm, s);

// Read in first two columns of space-delimited file.
// Each column is interpreted as a double.
interpolantd j("file.txt");

int main()
{
   cout << i(500 * nm) << endl;
   cout << j(500) << endl;
   // Class interpolant knows how to integrate itself exactly.
   cout << i.integral(300 * nm, 900 * nm) << endl;

   // Product of interpolants is interpolant with no loss of information.
   interpolant<length, double> k("file2.txt", nm, 1.0);
   interpolant<length, num::time> m = i * k;

   return 0;
}
```

A linear interpolant can also be initialized from an existing, continuous
function.

```.cpp
#include <cmath> // for erf()
#include <interpolant.hpp>

using namespace num;

int main()
{
   // 'interpolantd' is an alias for 'interpolant<double, double>'.
   interpolantd const i(erf, -1.0, 2.0);

   // Now you can use, for example i(1.3), to do a quick table-lookup of erf(),
   // etc.  This is useful especially for a continuous function that takes a
   // long time to calculate, when it has to be used over and over again in a
   // loop.

   return 0;
}
```

## Integration

```cpp
#include <iostream>
#include <integral.hpp>
#include <units.hpp>

using namespace num;
using namespace std;

area square(length x) { return x * x; }

int main()
{
   volume const i = integral(square, 0 * cm, 1 * cm);
   cout << i << endl; // "[3.33333e-07 m^3]"
   return 0;
}
```

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

