
# numeric

C++-based library for numerical computation

Included at present are physically dimensioned quantities, linear
interpolation, and integration by adaptive quadrature.  These interoperate.

At the moment, a recent clang++ seems to be required because g++ seems to limit
an identifier to consist of ASCII characters.  However, numeric uses the
character μ in more than one identifier.

