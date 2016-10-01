
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   util.hpp
/// \brief  Basic utilities.

#ifndef NUMERIC_UTIL_HPP
#define NUMERIC_UTIL_HPP

namespace num
{
   /// Type of ratio of two tings.
   /// \tparam Y  Type of numerator.
   /// \tparam X  Type of denominator.
   template <typename Y, typename X>
   using RAT = decltype(Y() / X());

   /// Type of product of two things.
   /// \tparam X  Type of left factor.
   /// \tparam Y  Type of right factor.
   template <typename X, typename Y>
   using PRD = decltype(X() * Y());

   /// Integer-template power of a double.
   /// \tparam P  Exponent.
   template <int P>
   double itpow(/** Base. */ double x)
   {
      double r = 1.0;
      if (P != 0) {
         int constexpr UP = (P > 0 ? P : -P);
         if (UP % 2 == 1) {
            r = x * itpow<UP - 1>(x);
         } else {
            r = itpow<UP / 2>(x * x);
         }
      }
      return P >= 0 ? r : 1.0 / r;
   }

   /// Integer-template power of a double.  This is an alias for \ref itpow().
   /// \tparam P  Exponent.
   template <int P>
   double pow(/** Base. */ double x)
   {
      return itpow<P>(x);
   }

   /// Integer power of a double.
   inline double ipow(/** Base. */ double x, /** Exponent. */ int e)
   {
      double r = 1.0;
      if (e != 0) {
         int const ue = (e > 0 ? e : -e);
         if (ue % 2 == 1) {
            r = x * ipow(x, ue - 1);
         } else {
            r = ipow(x * x, ue / 2);
         }
      }
      return e >= 0 ? r : 1.0 / r;
   }

   /// Integer power of a double.  This is an alias for \ref ipow().
   inline double pow(/** Base. */ double x, /** Exponent. */ int e)
   {
      return ipow(x, e);
   }
}

#endif // ndef NUMERIC_UTIL_HPP

