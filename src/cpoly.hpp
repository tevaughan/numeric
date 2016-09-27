
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   cpoly.hpp
/// \brief  Definition of num::cpoly.

#ifndef NUMERIC_CPOLY_HPP
#define NUMERIC_CPOLY_HPP

#include <array>   // for array
#include <cstring> // for memset()
#include <vector>  // for vector

#include <cfunc.hpp>

namespace num
{
   /// Model of a polynomial of a continuous variable.
   ///
   /// This class is not fully compatible with statdim.  The type \a C of each
   /// term may be a statdim, in which case type \a V must be dyndim.
   /// Similarly, if type \a C be dyndim, then type \a V must be dyndim.
   ///
   /// Class cpoly inherits from cfunc because cpoly is a model of a continuous
   /// function.
   ///
   /// \tparam D  Degree of polynomial.
   /// \tparam V  Type of variable.
   /// \tparam C  Type of each term.
   template <unsigned D, typename V=double, typename C=double>
   class cpoly : public cfunc<V, C, cpoly<D, V, C>>
   {
      enum {
         N = D + 1 ///< Number of coefficients.
      };

      /// Coefficients of polynomial.  Each of these is of the same type as the
      /// variable.  So if the variable be of type dyndim, then so will each
      /// coefficient be.
      std::array<V, N> c_;

   public:
      /// By default, initialize coefficients to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Copy coefficients from an array.
      cpoly(std::array<V, N> const& cc) : c_(cc) {}

      /// Copy coefficients from a vector.
      cpoly(std::vector<V> const& cc)
      {
         if (cc.size() != c_.size()) {
            throw "Wrong number of coefficients in vector initializer.";
         }
         for (unsigned i = 0; i < c_.size(); ++i) {
            c_[i] = cc[i];
         }
      }

      /// Return reference to mutable coefficients.
      std::array<V, N>& c() { return c_; }

      /// Return reference to immutable coefficients.
      std::array<V, N> const& c() const { return c_; }

      /// Type of derivative.
      using DERIV = cpoly<D - 1, V, decltype(C() / V())>;

      /// Type of integral.
      using INTEG = cpoly<D + 1, V, decltype(C() * V())>;

      /// Evaluate polynomial.
      C operator()(V const& v) const
      {
         C val = c_[0];
         V vp = v;
         for (unsigned i = 1; i < c_.size(); ++i) {
            val += c_[i] * vp;
            vp = v * vp;
         }
      }

      /// Return function representing derivative.
      DERIV derivative() const
      {
         DERIV d;
         for (unsigned i = 1; i < c_.size(); ++i) {
            d.c()[i - 1] = i * c_[i];
         }
         return d;
      }

      /// Return function representing integral from specified lower bound.
      INTEG integral(V const &lb) const
      {
         INTEG i;
         i.c()[0] = lb;
         for (unsigned j = 0; j < c_.size(); ++j) {
            i.c()[j + 1] = c_[j] / (j + 1);
         }
         return i;
      }
   };

   /// Specialization of cpoly for degree zero, which has a single coefficient,
   /// the constant term.
   template <typename V, typename C>
   class cpoly<0, V, C> : public cfunc<V, C, cpoly<0, V, C>>
   {
      C c_; ///< Coefficient of unity.

   public:
      /// Initialize coefficient.
      cpoly(C const &cc) : c_(cc) {}

      /// Type of integral.
      using INTEG = cpoly<1, V, decltype(C() * V())>;

      /// Return function representing integral from specified lower bound.
      INTEG integral(V const &lb) const
      {
         INTEG i;
         i.c()[0] = lb;
         i.c()[1] = c_;
         return i;
      }
   };
}

#endif // ndef NUMERIC_CPOLY_HPP

