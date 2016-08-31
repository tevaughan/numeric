
# numeric

## C++-11 library for numerical computation

Included at present are physically dimensioned quantities, linear
interpolation, and integration by adaptive quadrature.

These inter-operate!

To do:

 - Make new class whose constructor takes a smooth function and a range of its
   argument.  After construction, the instance will contain a linear
   interpolant for each of the function's first derivative, the function
   itself, and the function's integral.

 - Incorporate and adapt code from my linear regression project on github.

