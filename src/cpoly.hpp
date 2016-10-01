
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
#include <util.hpp>

namespace num
{
   /// Model of a polynomial of a continuous variable.
   ///
   /// Class cpoly inherits from \ref cfunc because cpoly is a model of a
   /// continuous function.
   ///
   /// \tparam D  Degree of polynomial.
   /// \tparam V  Type of variable. By default, primitive type double.
   /// \tparam C  Type of each term. By default, primitive type double.
   template <unsigned D, typename V = double, typename C = double>
   class cpoly : public cfunc<V, C, cpoly<D, V, C>>
   {
      static unsigned constexpr N = D + 1; ///< Number of coefficients.
      std::array<C, N> c_;                 ///< Normalized coefficients.

   public:
      /// By default, initialize coefficients to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Initialize coefficients.
      template <typename T>
      cpoly(std::array<T, N> const &a)
      {
         for (unsigned i = 0; i < N; ++i) {
            c_[i] = a[i] * pow(V(1.0), i);
         }
      }

      /// Initialize coefficients.
      template <typename T>
      cpoly(std::vector<T> const &v)
      {
         if (v.size() != N) {
            throw "Wrong number of coefficients.";
         }
         for (unsigned i = 0; i < N; ++i) {
            c_[i] = v[i] * pow(V(1.0), i);
         }
      }

      /// Return number of coefficients in polynomial.
      static unsigned constexpr num_coefs() { return N; }

      /// Type of coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      using ctype = decltype(C(1.0) / pow<I>(V(1.0)));

      /// Fetch coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      ctype<I> coef() const
      {
         static_assert(I <= D, "Array access must be in bounds.");
         return c_[I] / pow<I>(V(1.0));
      }

      /// Set coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      void set_coef(ctype<I> const &c)
      {
         static_assert(I <= D, "Array access must be in bounds.");
         c_[I] = c * pow<I>(V(1.0));
      }

      /// Evaluate polynomial.
      C operator()(/** Value of variable in polynomial. */ V const &v) const
      {
         double const vn = v / V(1.0); // Normalized value of variable.
         double vp = vn;               // Initial power of normalized value.
         C val = c_[0];
         for (unsigned i = 1; i < c_.size(); ++i) {
            val += c_[i] * vp;
            vp *= vn;
         }
         return val;
      }

      using dterm = decltype(C(1.0) / V(1.0)); ///< Type of derivative's term.
      using iterm = decltype(C(1.0) * V(1.0)); ///< Type of integral's term.
      using deriv = cpoly<D - 1, V, dterm>;    ///< Type of derivative.
      using integ = cpoly<D + 1, V, iterm>;    ///< Type of integral.
      friend class cpoly<D - 1, V, dterm>;     ///< Integral needs access.
      friend class cpoly<D + 1, V, iterm>;     ///< Derivative needs access.

      /// Return function representing derivative.
      deriv derivative() const
      {
         deriv d;
         for (unsigned i = 1; i < c_.size(); ++i) {
            d.c_[i - 1] = i * c_[i] / V(1.0);
         }
         return d;
      }

      /// Return function representing integral from specified lower bound.
      integ integral(/** Lower bound of integration. */ V const &lb) const
      {
         integ i;
         // Initialize constant with zero in right units.
         i.c_[0] = 0.0 * (*this)(lb)*lb;
         double lbn = lb / V(1.0); // Initialize power of normalized lower bound.
         for (unsigned j = 0; j < c_.size(); ++j) {
            iterm const cn = c_[j] / (j + 1) * V(1.0);
            i.c_[j + 1] = cn;
            i.c_[0] -= cn * lbn;
            lbn *= lb;
         }
         return i;
      }
   };

   /// Degree-zero specialization of \ref cpoly.  %cpoly<0,V,C> has but a
   /// single coefficient, the constant term.
   template <typename V, typename C>
   class cpoly<0, V, C> : public cfunc<V, C, cpoly<0, V, C>>
   {
      static unsigned constexpr N = 0 + 1; ///< Number of coefficients.
      std::array<C, N> c_;                 ///< Coefficient of unity.

   public:
      /// By default, initialize coefficient to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Copy coefficient from instance of its type.
      cpoly(C const &cc) { c_[0] = cc; }

      /// Return number of coefficients in polynomial.
      static unsigned constexpr num_coefs() { return N; }

      /// Fetch coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      C coef() const
      {
         static_assert(I == 0, "Array access must be in bounds.");
         return c_[I];
      }

      /// Set coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      void set_coef(C const &c)
      {
         static_assert(I == 0, "Array access must be in bounds.");
         c_[I] = c;
      }

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

      /// Evaluate polynomial. Ignore independent variable because present
      /// specialization models a constant function.
      C operator()(V const &) const { return c_[0]; }

      using iterm = decltype(C(1.0) * V(1.0)); ///< Type of integral's term.
      using integ = cpoly<1, V, iterm>;        ///< Type of integral.
      friend class cpoly<1, V, iterm>;         ///< Integral needs access here.

      /// Return function representing integral from specified lower bound.
      integ integral(V const &lb) const
      {
         integ i;
         i.c_[0] = -c_[0] * lb;
         i.c_[1] = +c_[0] / V(1.0);
         return i;
      }
   };
}

#include <dimval.hpp>

namespace num
{
   /// Disallowed specialization of \ref cpoly.
   template <unsigned D, typename V>
   class cpoly<D, V, dyndim>
   {
      cpoly(); ///< Only constructor is unimplemented and private.
   };

   /// Disallowed specialization of \ref cpoly.
   template <unsigned D, typename C>
   class cpoly<D, dyndim, C>
   {
      cpoly(); ///< Only constructor is unimplemented and private.
   };

   /// Specialization of \ref cpoly for \a V = \ref dyndim and \a C = \ref
   /// dyndim.
   ///
   /// Note carefully that, in order to use dyndim, both the variable type and
   /// the term time must be \ref dyndim.  One must initialize the
   /// coefficients, in order from 0 to \a D, inclusive, before using any
   /// facility.
   ///
   /// \tparam D  Degree of polynomial.
   template <unsigned D>
   class cpoly<D, dyndim, dyndim>
         : public cfunc<dyndim, dyndim, cpoly<D, dyndim, dyndim>>
   {
      static unsigned constexpr N = D + 1; ///< Number of coefficients.
      std::array<dyndim, N> c_;            ///< Coefficients.

   public:
      /// By default, initialize coefficients to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Initialize from array of coefficients.
      cpoly(std::array<dyndim, N> const& cc) : c_(cc) {}

      /// Initialize from vector of coefficients.
      cpoly(std::vector<dyndim> const& cc)
      {
         if (cc.size() != c_.size()) {
            throw "Initializer is wrong length.";
         }
         for (unsigned i = 0; i < c_.size(); ++i) {
            c_[i] = cc[i];
         }
      }

      /// Return number of coefficients in polynomial.
      static unsigned constexpr num_coefs() { return N; }

      /// Fetch coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      dyndim coef() const
      {
         static_assert(I <= D, "Array access must be in bounds.");
         return c_[I];
      }

      /// Fetch coefficient for term of degree \a i.
      dyndim const &coef(/** Degree of term. */ unsigned i) const
      {
         return c_[i];
      }

      /// Set coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      void set_coef(dyndim const &c)
      {
         static_assert(I <= D, "Array access must be in bounds.");
         c_[I] = c;
      }

      /// Set coefficient for term if degree \a i.
      dyndim &coef(/** Degree of term. */ unsigned i) { return c_[i]; }

      /// Evaluate polynomial.
      dyndim operator()(/** Value of variable. */ dyndim const &v) const
      {
         dyndim vp = v; // Initial power of variable's value.
         dyndim val = c_[0];
         for (unsigned i = 1; i < c_.size(); ++i) {
            val += c_[i] * vp;
            vp *= v;
         }
         return val;
      }

      using deriv = cpoly<D - 1, dyndim, dyndim>; ///< Type of derivative.
      using integ = cpoly<D + 1, dyndim, dyndim>; ///< Type of integral.
      friend class cpoly<D - 1, dyndim, dyndim>;  ///< Integral needs access.
      friend class cpoly<D + 1, dyndim, dyndim>;  ///< Derivative needs access.

      /// Return function representing derivative.
      deriv derivative() const
      {
         deriv d;
         for (unsigned i = 1; i < c_.size(); ++i) {
            d.c_[i - 1] = i * c_[i];
         }
         return d;
      }

      /// Return function representing integral from specified lower bound.
      integ integral(/** Lower bound. */ dyndim const &lb) const
      {
         integ i;
         // Initialize constant with zero in right units.
         i.c_[0] = 0.0 * (*this)(lb)*lb;
         dyndim lbn = lb; // Initialize power of lower bound.
         for (unsigned j = 0; j < c_.size(); ++j) {
            dyndim const cn = c_[j] / (j + 1);
            i.c_[j + 1] = cn;
            i.c_[0] -= cn * lbn;
            lbn *= lb;
         }
         return i;
      }
   };

   /// Degree-zero specialization of \ref cpoly with \a V = \ref dyndim and \a
   /// C = \ref dyndim.  \c cpoly<0,dyndim,dyndim> has but a single
   /// coefficient, the constant term.
   ///
   /// Note that if either of \a V and \a C be \ref dyndim, then the other must
   /// be, too.
   template <>
   class cpoly<0, dyndim, dyndim>
         : public cfunc<dyndim, dyndim, cpoly<0, dyndim, dyndim>>
   {
      static unsigned constexpr N = 0 + 1; ///< Number of coefficients.
      std::array<dyndim, N> c_;            ///< Coefficient of unity.

   public:
      /// By default, initialize coefficient to zero.
      cpoly() { memset(c_.data(), 0, sizeof(c_)); }

      /// Copy coefficient from instance of its type.
      cpoly(dyndim const &cc) { c_[0] = cc; }

      /// Return number of coefficients in polynomial.
      static unsigned constexpr num_coefs() { return N; }

      /// Fetch coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      dyndim coef() const
      {
         static_assert(I == 0, "Array access must be in bounds.");
         return c_[I];
      }

      /// Set coefficient for term of degree \a I.
      /// \tparam I  Degree of term.
      template <unsigned I>
      void set_coef(dyndim const &c)
      {
         static_assert(I == 0, "Array access must be in bounds.");
         c_[I] = c;
      }

      /// Allow assignment from instance of coefficient's type.
      cpoly &operator=(dyndim const &cc)
      {
         c_[0] = cc;
         return *this;
      }

      /// Enable default assignment.
      cpoly &operator=(cpoly const &cp) = default;

      /// Convert to reference to mutable coefficient.
      operator dyndim &() { return c_[0]; }

      /// Convert to reference to immutable coefficient.
      operator dyndim const &() const { return c_[0]; }

      /// Evaluate polynomial. Ignore independent variable because present
      /// specialization models a constant function.
      dyndim operator()(dyndim const &) const { return c_[0]; }

      using integ = cpoly<1, dyndim, dyndim>; ///< Type of integral.
      friend class cpoly<1, dyndim, dyndim>;  ///< Integral needs access here.

      /// Return function representing integral from specified lower bound.
      integ integral(dyndim const &lb) const
      {
         integ i;
         i.c_[0] = -c_[0] * lb;
         i.c_[1] = +c_[0];
         return i;
      }
   };
}

#endif // ndef NUMERIC_CPOLY_HPP

