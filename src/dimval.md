
# Physically Dimensioned Quantities {#dimvals}

Base class [dimval](\ref num::dimval) has public default constructor for
zero-valued quantity.

```cpp
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
```cpp
#include <units.hpp>
num::speed foo; // A speed with initial value of zero.
```

A wide array of units can be used for initialization and computation.
```cpp
#include <units.hpp>
num::speed v = 2.3 * mi / hr;  // Speed with initial value 2.3 miles/hour.
```

sqrt() and fabs() are supported.

Member functions allow for integer power or root.
```cpp
#include <units.hpp>
using namespace num;
length x = 2 * cm;
area a = pow<2>(x);
length y = 0.5 * sqrt(a);
volume v = pow<3>(x);
length z = 0.1 * root<3>(v);
```

