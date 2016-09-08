
# numeric

## C++-11 library for numerical computation

Included at present are physically dimensioned quantities, linear
interpolation, and integration by adaptive quadrature.

 - Linear interpolation and integration are header-only template facilities.

 - Physically dimensioned quantities require linking with '-lnumeric'.

See [the releases page](https://github.com/tevaughan/numeric/releases).

Grab the 'tar.gz' file.

Build and install:
```
tar fvxz numeric-x.y.tar.gz
cd numeric-x.y
./configure
make
sudo make install
```

All classes and functions inter-operate!

## Physically dimensioned quantities.

Base class 'dimval' has public default constructor for zero-valued quantity.
```c++
#include <units.hpp>
dimval<-1,1,0,0,0> foo; // A speed with initial value of zero.
```
Each of dimval's template value parameters is an integer representing the
exponent of one of five fundamental units.

1. second
2. meter
3. kilogram
4. coulomb
5. kelvin

Aliases are available for convenience.
```c++
#include <units.hpp>
num::speed foo; // A speed with initial value of zero.
```

A wide array of units can be used for initialization and computation.
```c++
#include <units.hpp>
num::speed v = 2.3 * mi / hr;  // Speed with initial value 2.3 miles/hour.
```

sqrt() and fabs() are supported.

Member functions allow for integer power or root.
```c++
#include <units.hpp>
using namespace num;
length x = 2 * cm;
area a = x.pow<2>();
length y = 0.5 * sqrt(a);
volume v = x.pow<3>();
length z = 0.1 * v.root<3>();
```

## Linear interpolation

A linear interpolant can be initialized directly from a file with two
space-delimited columns.  Then it can be scaled, multiplied by another
interpolant with no loss of information, integrated, etc.

```c++
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

```c++
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

```c++
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

## To do:

 - Incorporate and adapt code from my linear regression project on github.

# License

Copyright 2016
Thomas E. Vaughan

Distributable under the terms of the GNU LGPL, Version 3 or later.

