
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   dimval.hpp
/// \brief  Source code for num::dimval and for conversions dealing with angle.

#ifndef NUMERIC_DIMVAL_HPP
#define NUMERIC_DIMVAL_HPP

#include <cmath>      // for M_PI
#include <functional> // for function
#include <iostream>   // for ostream
#include <utility>    // for forward()

#include <template.hpp> // for RAT and PRD

namespace num
{
   double constexpr deg = M_PI / 180.0;     ///< Radians per degree.
   double constexpr arcmin = deg / 60.0;    ///< Radians per arcminute.
   double constexpr arcsec = arcmin / 60.0; ///< Radians per arcsecond.

   /// \return   Number of radians, converted from number of degrees.
   inline double degs(double n /**< Number of degrees. */) { return n * deg; }

   /// \return   Number of radians, converted from number of arcminutes.
   inline double arcmins(double n /**< Number of arcminutes. */)
   {
      return n * arcmin;
   }

   /// \return   Number of radians, converted from number of arcseconds.
   inline double arcsecs(double n /**< Number of arcseconds. */)
   {
      return n * arcsec;
   }

   // Forward declaration needed to allow interpolant to call A(0.0).  The
   // dimval's friend declaration for interpolant seems to work only if this
   // forward declaration be present.
   template <typename A, typename R>
   class interpolant;

   /// Model of a dimensioned value.
   /// \tparam TI  Exponent of time
   /// \tparam D   Exponent of distance.
   /// \tparam M   Exponent of mass
   /// \tparam C   Exponent of charge.
   /// \tparam TE  Exponent of temperature.
   template <int TI, int D, int M, int C, int TE>
   class dimval
   {
   protected:
      double v_; ///< Value in MKS.

      /// Construct from double that is known to contain value in MKS.
      explicit dimval(double vv /** Numeric coefficient of MKS unit. */)
            : v_(vv)
      {
      }

   public:
      /// By default, construct a zero-valued quantity.
      dimval() : v_(0.0) {}

      /// Make every kind of dimval be a friend to every other.
      template <int OTI, int OD, int OM, int OC, int OTE>
      friend class dimval;

      template <typename I>
      friend class integral_stats;

      /// Type of std::function that can be integrated.  A single-argument
      /// function is required.
      template <typename R, typename A>
      using func = std::function<R(A)>;

      /// Allow interpolant to construct from known MKS quantity.
      template <typename A, typename R>
      friend class interpolant;

      /// Allow integral() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend PRD<R, A> integral(func<R, A> f, A1 a, A2 b, double t,
                                unsigned n);

      /// Allow integral_rk() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend PRD<R, A> integral_rk(func<R, A> f, A1 x1, A2 x2, double t, int n,
                                   interpolant<A, R> *fi,
                                   interpolant<A, PRD<R, A>> *ii);

      /// Allow rkqs() to construct from known MKS quantity.
      template <typename X, typename Y>
      friend void rkqs(Y &y, decltype(Y() / X()) const &dydx, X &x,
                       X const &htry, double eps, Y const &yscal, X &hdid,
                       X &hnext, func<RAT<Y, X>, X> deriv);

      /// Add dimensioned values.
      dimval operator+(dimval dv) const { return dimval(v_ + dv.v_); }

      /// Unary position.
      friend dimval operator+(dimval dv) { return dv; }

      /// Unary negation.
      friend dimval operator-(dimval dv) { return dimval(-dv.v_); }

      /// Subtract dimensioned values.
      dimval operator-(dimval dv) const { return dimval(v_ - dv.v_); }

      /// Additive assignment.
      dimval &operator+=(dimval dv)
      {
         v_ += dv.v_;
         return *this;
      }

      /// Subtractive assignment.
      dimval &operator-=(dimval dv)
      {
         v_ -= dv.v_;
         return *this;
      }

      /// Multiply oppositely dimensioned values.
      double operator*(dimval<-TI, -D, -M, -C, -TE> dv) const
      {
         return v_ * dv.v_;
      }

      /// Multiply dimensioned values.
      template <int OTI, int OD, int OM, int OC, int OTE>
      dimval<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>
      operator*(dimval<OTI, OD, OC, OM, OTE> dv) const
      {
         return dimval<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>(v_ * dv.v_);
      }

      /// Multiplication by number on right side.
      dimval operator*(double s) const { return dimval(v_ * s); }

      /// Multiplication of dimval by number on left side.
      friend dimval operator*(double s, dimval dv)
      {
         return dimval(s * dv.v_);
      }

      /// Division resulting in double.
      double operator/(dimval dv) const { return v_ / dv.v_; }

      /// Divide dimensioned values.
      template <int OTI, int OD, int OM, int OC, int OTE>
      dimval<TI - OTI, D - OD, M - OM, C - OC, TE - OTE>
      operator/(dimval<OTI, OD, OM, OC, OTE> dv) const
      {
         return dimval<TI - OTI, D - OD, M - OM, C - OC, TE - OTE>(v_ / dv.v_);
      }

      /// Divide number by dimensioned value.
      friend dimval<-TI, -D, -M, -C, -TE> operator/(double s, dimval dv)
      {
         return dimval(dv.v_ / s).pow<-1>();
      }

      /// Division by number.
      dimval operator/(double s) const { return dimval(v_ / s); }

      /// Multiplicative assignment against double.
      dimval &operator*=(double s)
      {
         v_ *= s;
         return *this;
      }

      /// Divisive assignment against double.
      dimval &operator/=(double s)
      {
         v_ /= s;
         return *this;
      }

      /// Less-than comparison.
      bool operator<(dimval dv) const { return v_ < dv.v_; }

      /// Greater-than comparison.
      bool operator>(dimval dv) const { return v_ > dv.v_; }

      /// Less-than-or-equal-to comparison.
      bool operator<=(dimval dv) const { return v_ <= dv.v_; }

      /// Greater-than-or-equal-to comparison.
      bool operator>=(dimval dv) const { return v_ >= dv.v_; }

      /// Equality comparison.
      bool operator==(dimval dv) const { return v_ == dv.v_; }

      /// Inequality comparison.
      bool operator!=(dimval dv) const { return v_ != dv.v_; }

      /// Type of integer power of present instance.
      /// \tparam E  Integer exponent indicating power.
      template <int E>
      using dimval_power = dimval<TI * E, D * E, M * E, C * E, TE * E>;

      /// \return Integer power.
      template <int E>
      dimval_power<E> pow() const
      {
         return dimval_power<E>(std::pow(v_, E));
      }

      /// \return Integer power.
      template <int E>
      friend dimval_power<E> pow(dimval const &dv)
      {
         return dv.pow<E>();
      }

      /// Type of integer root of present instance.
      /// \tparam E  Integer indicating degree of root.
      template <int E>
      using dimval_root = dimval<TI / E, D / E, M / E, C / E, TE / E>;

      /// \return Integer root.
      template <int E>
      dimval_root<E> root() const
      {
         static_assert(E > 0, "zero or negative root");
         static_assert(D / E * E == D, "illegal root along distance");
         static_assert(M / E * E == M, "illegal root along mass");
         static_assert(C / E * E == C, "illegal root along charge");
         static_assert(TI / E * E == TI, "illegal root along time");
         static_assert(TE / E * E == TE, "illegal root along temperature");
         double constexpr e = 1.0 / E;
         return dimval_root<E>(std::pow(v_, e));
      }

      /// \return Integer root.
      template <int E>
      friend dimval_root<E> root(dimval const &dv)
      {
         return dv.root<E>();
      }

      /// \return Square root.
      friend dimval_root<2> sqrt(dimval dv) { return dv.root<2>(); }

      /// \return Square root.
      dimval_root<2> sqrt() const { return root<2>(); }

      /// \return Absolute value.
      friend dimval fabs(dimval dv) { return dimval(fabs(dv.v_)); }

      /// Write dimensioned value to output stream.
      friend std::ostream &operator<<(std::ostream &os, dimval dv)
      {
         os << "[" << dv.v_;
         if (M == 1) {
            os << " kg";
         } else if (M != 0) {
            os << " kg^" << M;
         }
         if (D == 1) {
            os << " m";
         } else if (D != 0) {
            os << " m^" << D;
         }
         if (TI == 1) {
            os << " s";
         } else if (TI != 0) {
            os << " s^" << TI;
         }
         if (C == 1) {
            os << " C";
         } else if (C != 0) {
            os << " C^" << C;
         }
         if (TE == 1) {
            os << " K";
         } else if (TE != 0) {
            os << " K^" << TE;
         }
         return os << "]";
      }
   };
}

#endif // ndef NUMERIC_DIMVAL_HPP

