
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

   template <typename A, typename R>
   class ilist;

   template <typename X, typename Y>
   class rk_quad;

   /// Base class for dimval and dyndim.
   ///
   /// \tparam DER  Type of descendant, either dimval or dyndim.  This is the
   ///              curiously recurring template pattern (CRTP).
   template <typename DER>
   class dimval_base
   {
   protected:
      double v_; ///< Value in MKS.

      /// Construct from double that is known to contain value in MKS.
      explicit dimval_base(double vv /**< Numeric coefficient of MKS unit. */)
         : v_(vv)
      {
      }

      /// Return reference to present instance as instance of DER.
      DER &d() { return *static_cast<DER *>(this); }

      /// Return reference to present instance as instance of DER const.
      DER const &d() const { return *static_cast<DER const *>(this); }

      /// Return reference to derived instance as base type.
      static dimval_base &base(DER &der)
      {
         return static_cast<dimval_base &>(der);
      }

      /// Return reference to derived instance as base type.
      static dimval_base const &base(DER const &der)
      {
         return static_cast<dimval_base const &>(der);
      }

   public:
      /// By default, construct a zero-valued quantity.
      dimval_base() : v_(0.0) {}

      /// Unary position.
      friend DER const &operator+(dimval_base const &dv) { return dv.d(); }

      /// Unary negation.
      friend DER operator-(dimval_base const &dv)
      {
         DER r = dv.d(); // Copy derived type in case of dyndim.
         base(r).v_ = -dv.v_;
         return r;
      }

      /// Multiplication by number on right side.
      DER operator*(double s) const
      {
         DER r = d(); // Copy derived type in case of dyndim.
         base(r).v_ *= s;
         return r;
      }

      /// Multiplication of dimval by number on left side.
      friend DER operator*(double s, dimval_base const &dv)
      {
         DER r = dv.d(); // Copy derived type in case of dyndim.
         base(r).v_ *= s;
         return r;
      }

      /// Division by number.
      DER operator/(double s) const
      {
         DER r = d(); // Copy derived type in case of dyndim.
         base(r).v_ /= s;
         return r;
      }

      /// Multiplicative assignment against double.
      DER &operator*=(double s)
      {
         v_ *= s;
         return d();
      }

      /// Divisive assignment against double.
      DER &operator/=(double s)
      {
         v_ /= s;
         return d();
      }

      /// Add dimensioned values.
      template <typename ODER>
      DER operator+(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Addition requires same dimension.";
         }
         DER r = d();
         base(r).v_ += dv.v_;
         return r;
      }

      /// Subtract dimensioned values.
      template <typename ODER>
      DER operator-(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Subtraction requires same dimension.";
         }
         DER r = d();
         base(r).v_ -= dv.v_;
         return r;
      }

      /// Additive assignment.
      template <typename ODER>
      DER &operator+=(dimval_base<ODER> const &dv)
      {
         if (!d().same_dim(dv.d())) {
            throw "Addition requires same dimension.";
         }
         v_ += dv.v_;
         return d();
      }

      /// Subtractive assignment.
      template <typename ODER>
      DER &operator-=(dimval_base<ODER> const &dv)
      {
         if (!d().same_dim(dv.d())) {
            throw "Subtraction requires same dimension.";
         }
         v_ -= dv.v_;
         return d();
      }

      /// Less-than comparison.
      template <typename ODER>
      bool operator<(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ < dv.v_;
      }

      /// Greater-than comparison.
      template <typename ODER>
      bool operator>(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ > dv.v_;
      }

      /// Less-than-or-equal-to comparison.
      template <typename ODER>
      bool operator<=(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ <= dv.v_;
      }

      /// Greater-than-or-equal-to comparison.
      template <typename ODER>
      bool operator>=(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ >= dv.v_;
      }

      /// Equality comparison.
      template <typename ODER>
      bool operator==(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ == dv.v_;
      }

      /// Inequality comparison.
      template <typename ODER>
      bool operator!=(dimval_base<ODER> const &dv) const
      {
         if (!d().same_dim(dv.d())) {
            throw "Comparison requires same dimension.";
         }
         return v_ != dv.v_;
      }

      /// Write dimensioned value to output stream.
      friend std::ostream &operator<<(std::ostream &os, dimval_base const &dvb)
      {
         DER const& dv = dvb.d();
         os << "[" << dvb.v_;
         if (dv.eM() == 1) {
            os << " kg";
         } else if (dv.eM() != 0) {
            os << " kg^" << dv.eM();
         }
         if (dv.eD() == 1) {
            os << " m";
         } else if (dv.eD() != 0) {
            os << " m^" << dv.eD();
         }
         if (dv.eTI() == 1) {
            os << " s";
         } else if (dv.eTI() != 0) {
            os << " s^" << dv.eTI();
         }
         if (dv.eC() == 1) {
            os << " C";
         } else if (dv.eC() != 0) {
            os << " C^" << dv.eC();
         }
         if (dv.eTE() == 1) {
            os << " K";
         } else if (dv.eTE() != 0) {
            os << " K^" << dv.eTE();
         }
         return os << "]";
      }
   };

   /// Model of a dimensioned value.
   /// \tparam TI  Exponent of time
   /// \tparam D   Exponent of distance.
   /// \tparam M   Exponent of mass
   /// \tparam C   Exponent of charge.
   /// \tparam TE  Exponent of temperature.
   template <char TI, char D, char M, char C, char TE>
   class dimval : public dimval_base<dimval<TI, D, M, C, TE>>
   {
   protected:
      using PT = dimval_base<dimval>;
      using PT::PT;
      using PT::v_;

   public:
      using PT::operator*;
      using PT::operator/;

      static int constexpr eM() { return M; }   ///< Exponent of mass.
      static int constexpr eD() { return D; }   ///< Exponent of distance.
      static int constexpr eTI() { return TI; } ///< Exponent of time.
      static int constexpr eC() { return C; }   ///< Exponent of charge.
      static int constexpr eTE() { return TE; } ///< Exponent of temperature.

      /// Return true only if other dimval has same dimension.
      template <char OTI, char OD, char OM, char OC, char OTE>
      static bool constexpr same_dim(dimval<OTI, OD, OM, OC, OTE>)
      {
         return TI == OTI && D == OD && M == OM && C == OC && TE == OTE;
      }

      /// Make every kind of dimval be a friend to every other.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend class dimval;

      /// Allow integral_stats to construct from known MKS quantity.
      template <typename I>
      friend class integral_stats;

      /// Type of std::function that can be integrated.  A single-argument
      /// function is required.
      template <typename R, typename A>
      using func = std::function<R(A)>;

      /// Allow interpolant to construct from known MKS quantity.
      template <typename A, typename R>
      friend class interpolant;

      /// Allow rk_quad to construct from known MKS quantity.
      template <typename X, typename Y>
      friend class rk_quad;

      /// Allow integral() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend PRD<R, A> integral(func<R, A> f, A1 a, A2 b, double t,
                                unsigned n);

      /// Multiply oppositely dimensioned values.
      double operator*(dimval<-TI, -D, -M, -C, -TE> dv) const
      {
         return v_ * dv.v_;
      }

      /// Multiply dimensioned values.
      template <char OTI, char OD, char OM, char OC, char OTE>
      dimval<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>
      operator*(dimval<OTI, OD, OC, OM, OTE> dv) const
      {
         return dimval<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>(v_ * dv.v_);
      }

      /// Division resulting in double.
      double operator/(dimval dv) const { return v_ / dv.v_; }

      /// Divide dimensioned values.
      template <char OTI, char OD, char OM, char OC, char OTE>
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
   };
}

#endif // ndef NUMERIC_DIMVAL_HPP

