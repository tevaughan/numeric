
# Physically Dimensioned Quantities {#dimvals}

Class num::dimval is a model of a dimensioned quantity.  In the physical
sciences, every product of a number and a physical unit is a dimensioned
quantity.

A physical unit is not a number but can be multiplied by a number or by another
unit.  For example, if a ruler be twelve inches long, then one might write the
length of the ruler as \f$12~\text{in}\f$.  In this case, \f$\text{in}\f$
refers to a physical unit, an inch.  In the example, \f$12~\text{in}\f$
indicates that the number 12 is multiplied by an inch.  In ordinary language,
one has that twelve times an inch is twelve inches.  In this example, the
quantity under consideration is a length, which is one kind of dimensioned
quantity.  The same 12-inch length could be expressed as \f$30.48~\text{cm}\f$.
So different units can correspond to the same dimension, which in the present
example is the dimension of length.

Every instance of num::dimval corresponds to one or another dimension, such as
length, speed, mass, force, etc.  The dimension is specified by the template
value parameters of num::dimval, which is a template class.  Each of dimval's
template value parameters is an integer representing the exponent of one of
five fundamental units.

1. second
2. meter
3. kilogram
4. coulomb
5. kelvin

Class num::dimval has public default constructor for a zero-valued quantity.
```cpp
#include <units.hpp>
dimval<-1,1,0,0,0> foo; // A speed with initial value of zero.
```

For many a dimension, an alias is available for convenience.
```cpp
#include <units.hpp>
num::speed foo; // A speed with initial value of zero.
```

The definition of the aliased dimensions is stored in a text file
[dimensions.txt](https://github.com/tevaughan/numeric/blob/master/src/dimensions.txt)
that is used at library build time to generate the corresponding C++ file
dimensions.hpp.  This file can easily be edited before the library is built to
increase the number of known dimensions.

For initialization of a dimensioned quantity, a wide array of units can be used
for initialization and computation.  The definition for every unit known by the
library is stored in a text file
[units.txt](https://github.com/tevaughan/numeric/blob/master/src/dimensions.txt)
that is used at library build time to generate the corresponding C++ file
units.hpp.
```cpp
#include <units.hpp>
using namespace num;
speed v = 2.3 * u::mi / u::hr;  // Speed with initial value 2.3 miles/hour.
```

sqrt() and fabs() are supported.

A member function allows for an integer power, and another member function
allows for an integer root.
```cpp
#include <units.hpp>
using namespace num;
length x = 2 * u::cm;
area a = pow<2>(x);
length y = 0.5 * sqrt(a);
volume v = pow<3>(x);
length z = 0.1 * root<3>(v);
```

Standard-library output streams are supported.
```cpp
  #include <iostream>
  #include <units.hpp>
  using namespace num;
  int main() {
     using time = num::time;  // in case of conflict with C library's time()
     length x0 = 2 * u::m;         // meters
     speed  v  = 5 * u::m / u::s;  // meters per second
     time   t  = 3 * u::min;       // minutes
     // Output distance in meters (the default) and in yards.
     std::cout << x0 + v * t << std::endl;                  // "[902 m]"
     std::cout << (x0 + v * t) / yd << " yd" << std::endl;  // "986.439 yd"
     return 0;
  }
```

