
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#ifndef NUMERIC_DIMVAL_HPP
#define NUMERIC_DIMVAL_HPP

#include <cmath>      // for M_PI
#include <functional> // for function
#include <iostream>   // for ostream
#include <utility>    // for forward()

namespace num
{
   double constexpr deg = M_PI / 180.0;     ///< Radians per degree.
   double constexpr arcmin = deg / 60.0;    ///< Radians per arcminute.
   double constexpr arcsec = arcmin / 60.0; ///< Radians per arcsecond.

   /// \param n  Number of degrees.
   /// \return   Corresponding number of radians.
   inline double degs(double n) { return n * deg; }

   /// \param n  Number of arcminutes.
   /// \return   Corresponding number of radians.
   inline double arcmins(double n) { return n * arcmin; }

   /// \param n  Number of arcseconds.
   /// \return   Corresponding number of radians.
   inline double arcsecs(double n) { return n * arcsec; }

   /// Model of a dimensioned value.
   /// \tparam TI  Exponent of time.
   /// \tparam D   Exponent of distance.
   /// \tparam M   Exponent of mass.
   /// \tparam C   Exponent of charge.
   /// \tparam TE  Exponent of temperature.
   template <int TI, int D, int M, int C, int TE>
   class dimval
   {
   protected:
      double v_; ///< Value in MKS.

      /// Construct from double that is known to contain value in MKS.
      explicit dimval(double vv) : v_(vv) {}

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

      /// Allow integral() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend auto integral(func<R, A> f, A1 a, A2 b, double t, unsigned n)
            -> decltype(std::forward<func<R, A>>(f)(A()) * A());

      /// Allow interpolant to construct from known MKS quantity.
      template <typename A, typename R>
      friend class interpolant;

      /// Add dimensioned values.
      dimval operator+(dimval dv) const { return dimval(v_ + dv.v_); }

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

      /// \return Integer power.
      template <int E>
      dimval<TI * E, D * E, M * E, C * E, TE * E> pow() const
      {
         return dimval<TI * E, D * E, M * E, C * E, TE * E>(std::pow(v_, E));
      }

      /// \return Integer root.
      template <unsigned E>
      dimval<TI / E, D / E, M / E, C / E, TE / E> root() const
      {
         static_assert(E, "zeroth root");
         static_assert(D / E * E == D, "illegal root along distance");
         static_assert(M / E * E == M, "illegal root along mass");
         static_assert(C / E * E == C, "illegal root along charge");
         static_assert(TI / E * E == TI, "illegal root along time");
         static_assert(TE / E * E == TE, "illegal root along temperature");
         double constexpr e = 1.0 / E;
         return dimval<TI / E, D / E, M / E, C / E, TE / E>(std::pow(v_, e));
      }

      /// \return Square root.
      friend dimval<TI / 2, D / 2, M / 2, C / 2, TE / 2> sqrt(dimval dv)
      {
         return dv.root<2>();
      }

      /// \return Square root.
      dimval<TI / 2, D / 2, M / 2, C / 2, TE / 2> sqrt() const
      {
         return root<2>();
      }

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

