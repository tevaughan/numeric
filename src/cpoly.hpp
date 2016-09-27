
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
   /// However, the specialization for the zeroth-degree polynomial (the
   /// constant) is fully compatible with statdim.
   ///
   /// Class cpoly inherits from cfunc because cpoly is a model of a continuous
   /// function.
   ///
   /// \tparam D  Degree of polynomial.
   /// \tparam V  Type of variable. By default, primitive type double.
   /// \tparam C  Type of each term. By default, primitive type double.
   template <unsigned D, typename V = double, typename C = double>
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
      cpoly(std::array<V, N> const &cc) : c_(cc) {}

      /// Copy coefficients from a vector.
      cpoly(std::vector<V> const &cc)
      {
         if (cc.size() != N) {
            throw "Wrong number of coefficients in vector initializer.";
         }
         for (unsigned i = 0; i < N; ++i) {
            c_[i] = cc[i];
         }
      }

      /// Return reference to mutable coefficients.
      std::array<V, N> &c() { return c_; }

      /// Return reference to immutable coefficients.
      std::array<V, N> const &c() const { return c_; }

      /// Type of derivative.
      using DERIV = cpoly<D - 1, V, decltype(C() / V())>;

      /// Type of integral.
      using INTEG = cpoly<D + 1, V, decltype(C() * V())>;

      /// Evaluate polynomial.
      C operator()(V const &v /**< Independent variable. */) const
      {
         C val = c_[0];
         V vp = v; // Initialize pth power of independent variable.
         for (unsigned i = 1; i < c_.size(); ++i) {
            val += c_[i] * vp;
            vp *= v;
         }
         return val;
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
      INTEG integral(V const &lb /**< Lower bound of integration. */) const
      {
         INTEG i;
         // Initialize constant with zero in right units (in case type C be
         // dimval).
         i.c()[0] = 0.0 * (*this)(lb)*lb;
         V lbn = lb; // Initialize nth power of lower bound.
         for (unsigned j = 0; j < c_.size(); ++j) {
            V const cn = c_[j] / (j + 1);
            i.c()[j + 1] = cn;
            i.c()[0] -= cn * lbn;
            lbn *= lb;
         }
         return i;
      }
   };

   /// Specialization of cpoly for degree zero, which has a single coefficient,
   /// the constant term.
   template <typename V, typename C>
   class cpoly<0, V, C> : public cfunc<V, C, cpoly<0, V, C>>
   {
      enum {
         N = 0 + 1 ///< Number of coefficients.
      };

      std::array<C, N> c_; ///< Coefficient of unity.

   public:
      /// By default, initialize coefficient to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Copy coefficient from an array.
      cpoly(std::array<C, N> const &cc) : c_(cc) {}

      /// Copy coefficient from a vector.
      cpoly(std::vector<V> const &cc)
      {
         if (cc.size() != N) {
            throw "Wrong number of coefficients in vector initializer.";
         }
         c_[0] = cc[0];
      }

      /// Copy coefficient from instance of its type.
      cpoly(C const &cc) { c_[0] = cc; }

      /// Return reference to mutable coefficient.
      std::array<C, N> &c() { return c_; }

      /// Return reference to mutable coefficient.
      std::array<C, N> const &c() const { return c_; }

      /// Allow assignment from instance of coefficient's type.
      cpoly &operator=(C const &cc)
      {
         c_[0] = cc;
         return *this;
      }

      /// Enable default assignment.
      cpoly &operator=(cpoly const &cp) = default;

      /// Convert to reference to mutable coefficient.
      operator C &() { return c_[0]; }

      /// Convert to reference to immutable coefficient.
      operator C const &() const { return c_[0]; }

      /// Type of integral.
      using INTEG = cpoly<1, V, decltype(C() * V())>;

      /// Evaluate polynomial. Ignore independent variable because present
      /// specialization models a constant function.
      C operator()(V const &) const
      {
         return c_[0];
      }

      /// Return function representing integral from specified lower bound.
      INTEG integral(V const &lb) const
      {
         INTEG i;
         i.c()[0] = -c_[0] * lb;
         i.c()[1] = +c_[0];
         return i;
      }
   };
}

#endif // ndef NUMERIC_CPOLY_HPP

