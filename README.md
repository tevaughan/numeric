
# numeric

## C++-11 library for numerical computation

Included at present are physically dimensioned quantities, linear
interpolation, and integration by adaptive quadrature.

These inter-operate!

## Physically dimensioned quantities.

Base class 'dimval' has public default constructor for zero-valued quantity.
```c++
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
speed foo; // A speed with initial value of zero.
```

A wide array of units can be used for initialization and computation.
```c++
speed v = 2.3 * mi / hr;  // Speed with initial value 2.3 miles/hour.
```

sqrt() and fabs() are supported.

Member functions allow for integer power or root.
```c++
length x = 2 * cm;
area a = x.pow<2>();
length y = 0.5 * a.root<2>().
```
## Linear interpolation

TBS

## Integration

TBS

## To do:

 - Make new class whose constructor takes a smooth function and a range of its
   argument.  After construction, the instance will contain a linear
   interpolant for each of the function's first derivative, the function
   itself, and the function's integral.

 - Incorporate and adapt code from my linear regression project on github.

