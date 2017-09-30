
// Copyright 2016-2017  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   integral.hpp
/// \brief  Definition of num::integral().

#ifndef NUMERIC_INTEGRAL_HPP
#define NUMERIC_INTEGRAL_HPP

#include <rk.hpp> // for rk_quad

namespace num
{
   /// Numerically integrate a function, and return the result.
   ///
   /// Use fifth-order Runge-Kutta with adaptive stepsize. See
   /// rk_quad::rk_quad().
   ///
   /// The initial-guess parameter is the inverse fraction of the domain
   /// (between limits of integration) to use as the guess for the initial step
   /// size.
   ///
   /// \tparam X   Type of argument to function.
   /// \tparam X1  Type of lower limit of integration (convertible to X).
   /// \tparam X2  Type of upper limit of integration (convertible to X).
   /// \tparam Y   Type returned by function that is to be integrated.
   /// \return     Numeric integral of function.
   template <typename X, typename Y, typename X1, typename X2>
   PRD<X, Y> integral(
         /** Function to be integrated.  */ std::function<Y(X)> f,
         /** Lower limit of integration. */ X1                  aa,
         /** Upper limit of integration. */ X2                  bb,
         /** Error tolerance.            */ double              t = 1.0E-06,
         /** Initial guess parameter.    */ unsigned            n = 16)
   {
      return rk_quad<X, PRD<X,Y>>(f, aa, bb, t, n).def_int();
   }

   /// Numerically integrate a function, and return the result.
   ///
   /// Use fifth-order Runge-Kutta with adaptive stepsize. See
   /// rk_quad::rk_quad().
   ///
   /// The initial-guess parameter is the inverse fraction of the domain
   /// (between limits of integration) to use as the guess for the initial step
   /// size.
   ///
   /// \tparam X   Type of argument to function.
   /// \tparam X1  Type of lower limit of integration (convertible to X).
   /// \tparam X2  Type of upper limit of integration (convertible to X).
   /// \tparam Y   Type returned by function that is to be integrated.
   /// \return     Numeric integral of function.
   template <typename X, typename Y, typename X1, typename X2>
   PRD<X, Y> integral(
         /** Function to be integrated.               */ Y (*f)(X),
         /** Lower limit of integration.              */ X1       a,
         /** Upper limit of integration.              */ X2       b,
         /** Error tolerance.                         */ double   t = 1.0E-06,
         /** Initial number of evenly spaced samples. */ unsigned n = 16)
   {
      return num::integral(std::function<Y(X)>(f), a, b, t, n);
   }
}

#endif // ndef NUMERIC_INTEGRAL_HPP

