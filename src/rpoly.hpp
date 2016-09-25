
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

// FIXME: This should probably be cpoly, for continuous-variable polynomial.
// However, the complexity of template expressions is such that to do
// polynomials with units seems to require an implementation via dynamic
// dimval, which isn't done yet.

/// \file   rpoly.hpp
/// \brief  Definition of num::rpoly.

#ifndef NUMERIC_RPOLY_HPP
#define NUMERIC_RPOLY_HPP

#include <array> // for array
#include <cmath> // for pow()
#include <tuple> // for tuple, tuple_cat()

#include <rfunc.hpp>

namespace num
{
   /// Extend dimval-style pow() to regular double.
   template <int P>
   double pow(double x)
   {
      std::pow(x, P);
   }

   /// Model of a polynomial of a continuous variable.
   ///
   /// rpoly inherits from rfunc because rpoly is a model of a continuous
   /// function.
   ///
   /// \tparam A  Type of variable.
   /// \tparam T  Type of each term.
   /// \tparam D  Degree of polynomial.
   template <typename A, typename T, unsigned D>
   class rpoly : public rfunc<A, rpoly<A, T, D>>
   {
#if 0
      /// Type of coefficient in highest-degree term.
      using HC = decltype(T() / pow<D>(A()));

      HC hc_; ///< Coefficient in highest-degree term.

      /// Sub-polynomial with all terms but highest degree.
      rpoly<A, T, D - 1> lower_;

   public:
      /// Type of coefficient set.
      using COEFS =
            decltype(std::tuple_cat(std::make_tuple(hc_), lower_.coefs()));

      rpoly(HC const &h, Types... lower)
         : hc_(h), lower_(lower...)
      {
      }

      auto coefs() const
      {

      }

      /// Type of derivative.
      using DERIV = rpoly<A, decltype(T() / A()), D - 1>;

      /// Type of integral.
      using INTEG = rpoly<A, decltype(T() * A()), D + 1>;

      /// Return function representing derivative.
      DERIV derivative() const
      {
         return DERIV(D * hc_, 
      }

      /// Return function representing integral from specified lower bound.
      INTEG integral(A const &lb) const { return d()->integral(lb); }
#endif
   };

   template <typename A, typename T>
   class rpoly<A, T, 0> : public rfunc<A, rpoly<A, T, 0>>
   {
      T c_; ///< Coefficient of unity.

   public:
      /// Initialize coefficient.
      rpoly(T const &cc) : c_(cc) {}

      /// Type of integral.
      using INTEG = rpoly<A, decltype(T() * A()), 1>;

      /// Return function representing derivative.
      double derivative() const { return rpoly<A, double, 0>(0.0); }

      /// Return function representing integral from specified lower bound.
      INTEG integral(A const &lb) const { return INTEG(c_, c_ * lb); }
   };
}

#endif // ndef NUMERIC_RPOLY_HPP

