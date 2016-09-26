
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   dimval.hpp
///
/// \brief  Source code for conversions dealing with angle, for num::dim_exps,
///         for num::dimval_base, for num::statdim, and for num::dyndim.

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

   // Forward declaration needed to allow a class to declare interpolant as a
   // friend.
   template <typename A, typename R>
   class interpolant;

   // Forward declaration needed to allow a class to declare ilist as a friend.
   template <typename A, typename R>
   class ilist;

   // Forward declaration needed to allow a class to declare rk_quad as a
   // friend.
   template <typename X, typename Y>
   class rk_quad;

   /// Dimensional exponents.
   class dim_exps
   {
      union {
         uint64_t n_; ///< Dimensional exponents as single number.
         char e_[8];  ///< Dimensional exponents.
      };

   public:
      /// Initialize exponent for each dimension.
      dim_exps(char ti, ///< Exponent of time.
               char d,  ///< Exponent of distance.
               char m,  ///< Exponent of mass.
               char c,  ///< Exponent of electric charge.
               char te  ///< Exponent of temperature
               )
         : n_(0)
      {
         TI() = ti;
         D() = d;
         M() = m;
         C() = c;
         TE() = te;
      }

      /// Reference to mutable exponent of time.
      char &TI() { return e_[0]; }

      /// Reference to mutable exponent of distance.
      char &D() { return e_[1]; }

      /// Reference to mutable exponent of mass.
      char &M() { return e_[2]; }

      /// Reference to mutable exponent of electric charge.
      char &C() { return e_[3]; }

      /// Reference to mutable exponent of temperature.
      char &TE() { return e_[4]; }

      /// Reference to immutable exponent of time.
      char const &TI() const { return e_[0]; }

      /// Reference to immutable exponent of distance.
      char const &D() const { return e_[1]; }

      /// Reference to immutable exponent of mass.
      char const &M() const { return e_[2]; }

      /// Reference to immutable exponent of electric charge.
      char const &C() const { return e_[3]; }

      /// Reference to immutable exponent of temperature.
      char const &TE() const { return e_[4]; }

      /// Reference to immutable number representing all exponents.
      uint64_t const &n() const { return n_; }

      /// Add exponents for multiplication of two dimensioned quantities.
      dim_exps operator+(dim_exps ode /**< Other exponents. */) const
      {
         return dim_exps(TI() + ode.TI(), D() + ode.D(), M() + ode.M(),
                         C() + ode.C(), TE() + ode.TE());
      }

      /// Subtract exponents for division of two dimensioned quantities.
      dim_exps operator-(dim_exps ode /**< Other exponents. */) const
      {
         return dim_exps(TI() - ode.TI(), D() - ode.D(), M() - ode.M(),
                         C() - ode.C(), TE() - ode.TE());
      }

      /// Multiply exponents by integer as for integer power of dimensioned
      /// value.
      dim_exps operator*(int p /**< Integer power. */) const
      {
         return dim_exps(TI() * p, D() * p, M() * p, C() * p, TE() * p);
      }

      /// Divide exponents by integer as for integer root of dimensioned value.
      dim_exps operator/(int r /**< Integer root. */) const
      {
         if (r <= 0) {
            throw "zero or negative root";
         }
         if (D() / r * r != D()) {
            throw "illegal root along distance";
         }
         if (M() / r * r != M()) {
            throw "illegal root along mass";
         }
         if (C() / r * r != C()) {
            throw "illegal root along charge";
         }
         if (TI() / r * r != TI()) {
            throw "illegal root along time";
         }
         if (TE() / r * r != TE()) {
            throw "illegal root along temperature";
         }
         return dim_exps(TI() / r, D() / r, M() / r, C() / r, TE() / r);
      }

      /// Return true only if every exponent be same as correspondent in other
      /// set.
      bool operator==(dim_exps ode /**< Other exponents. */) const
      {
         return n_ == ode.n_;
      }

      /// Return true only if any exponent be different from correspondent in
      /// other set.
      bool operator!=(dim_exps ode /**< Other exponents. */) const
      {
         return n_ != ode.n_;
      }

      /// Return negative exponents as for reciprocal of dimensioned value.
      dim_exps neg() const { return dim_exps(-TI(), -D(), -M(), -C(), -TE()); }
   };

   /// Base class for statdim and dyndim.
   ///
   /// \tparam DER  Type of descendant, either statdim or dyndim.  This is the
   ///              curiously recurring template pattern (CRTP).
   template <typename DER>
   class dimval
   {
      /// Needed to allow statdim to access v_ in some circumstances.
      friend class dyndim;

      /// Needed to allow copying from one kind of dimval to another.
      template <typename ODER>
      friend class dimval;

   protected:
      double v_; ///< Value in MKS.

      /// Construct from double that is known to contain value in MKS.
      dimval(double vv /**< Numeric coefficient of MKS unit. */) : v_(vv) {}

      /// Copy from another base object. This needs to be done carefully in
      /// case the descendant has extra stuff of its own to copy. So copy
      /// constructor is protected.
      template <typename ODER>
      dimval(dimval<ODER> const &dvb)
         : v_(dvb.v_)
      {
      }

      /// Return reference to present instance as instance of DER.
      DER &d() { return *static_cast<DER *>(this); }

      /// Return reference to present instance as instance of DER const.
      DER const &d() const { return *static_cast<DER const *>(this); }

   public:
      /// Dimensional exponents.
      dim_exps exps() const { return d().exps(); }

      /// By default, construct a zero-valued quantity.
      dimval() : v_(0.0) {}

      /// Unary position.
      friend DER const &operator+(dimval const &dv) { return dv.d(); }

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

      /// Write dimensioned value to output stream.
      friend std::ostream &operator<<(std::ostream &os, dimval const &dvb)
      {
         dim_exps const e = dvb.exps();
         os << "[" << dvb.v_;
         if (e.M() == 1) {
            os << " kg";
         } else if (e.M() != 0) {
            os << " kg^" << int(e.M());
         }
         if (e.D() == 1) {
            os << " m";
         } else if (e.D() != 0) {
            os << " m^" << int(e.D());
         }
         if (e.TI() == 1) {
            os << " s";
         } else if (e.TI() != 0) {
            os << " s^" << int(e.TI());
         }
         if (e.C() == 1) {
            os << " C";
         } else if (e.C() != 0) {
            os << " C^" << int(e.C());
         }
         if (e.TE() == 1) {
            os << " K";
         } else if (e.TE() != 0) {
            os << " K^" << int(e.TE());
         }
         return os << "]";
      }
   };

   template <char TI, char D, char M, char C, char TE>
   class statdim;

   class dyndim;

   template <char TI, char D, char M, char C, char TE>
   statdim<TI, D, M, C, TE> &operator+=(statdim<TI, D, M, C, TE> &dv,
                                        dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   statdim<TI, D, M, C, TE> &operator-=(statdim<TI, D, M, C, TE> &dv,
                                        dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator==(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator!=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator<(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator<=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator>(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   bool operator>=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   /// Multiply statdim and dyndim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator*(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   /// Divide statdim by dyndim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator/(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

   /// Model of a statically dimensioned value.  For a dynamically dimensioned
   /// value, see dyndim.
   ///
   /// \tparam TI  Exponent of time
   /// \tparam D   Exponent of distance.
   /// \tparam M   Exponent of mass
   /// \tparam C   Exponent of charge.
   /// \tparam TE  Exponent of temperature.
   template <char TI, char D, char M, char C, char TE>
   class statdim : public dimval<statdim<TI, D, M, C, TE>>
   {
      /// Every kind of statdim must be a friend of every other.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend class statdim;

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

      friend statdim &operator+=
            <TI, D, M, C, TE>(statdim &dv, dyndim const &dd);

      friend statdim &operator-=
            <TI, D, M, C, TE>(statdim &dv, dyndim const &dd);

      friend bool operator==<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend bool operator!=<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend bool operator< //
            <TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend bool operator<=<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend bool operator><TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend bool operator>=<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend dyndim operator*<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      friend dyndim operator/<TI, D, M, C, TE>(statdim dv, dyndim const &dd);

      /// Dimensional exponents in a form useful for interaction with dyndim.
      static dim_exps const exps_;

   protected:
      using PT = dimval<statdim>; ///< Type of parent.
      using PT::PT;               ///< Inherit parent's constructor.
      using PT::v_;               ///< Inherit protected, numeric value.

   public:
      /// Dimensional exponents.
      static dim_exps exps() { return exps_; }

      /// Unary negation.
      friend statdim operator-(statdim dv) { return -dv.v_; }

      /// Multiply oppositely dimensioned values.
      double operator*(statdim<-TI, -D, -M, -C, -TE> dv) const
      {
         return v_ * dv.v_;
      }

      /// Multiplication by number on right side.
      statdim operator*(double s) const { return v_ * s; }

      /// Multiplication of statdim by number on left side.
      friend statdim operator*(double s, statdim dv) { return s * dv.v_; }

      /// Type of product of present statdim with other statdim.
      template <char OTI, char OD, char OM, char OC, char OTE>
      using prod = statdim<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>;

      /// Multiply dimensioned values.
      template <char OTI, char OD, char OM, char OC, char OTE>
      prod<OTI, OD, OM, OC, OTE>
      operator*(statdim<OTI, OD, OC, OM, OTE> dv) const
      {
         return v_ * dv.v_;
      }

      /// Division by number.
      statdim operator/(double s) const { return v_ / s; }

      /// Division resulting in double.
      double operator/(statdim dv) const { return v_ / dv.v_; }

      /// Type of quotient of present statdim and other statdim.
      template <char OTI, char OD, char OM, char OC, char OTE>
      using quot = statdim<TI - OTI, D - OD, M - OM, C - OC, TE - OTE>;

      /// Divide dimensioned values.
      template <char OTI, char OD, char OM, char OC, char OTE>
      quot<OTI, OD, OM, OC, OTE>
      operator/(statdim<OTI, OD, OM, OC, OTE> dv) const
      {
         return v_ / dv.v_;
      }

      /// Type of reciprocal of present statdim.
      using inv = statdim<-TI, -D, -M, -C, -TE>;

      /// Divide number by dimensioned value.
      friend inv operator/(double s, statdim dv)
      {
         return statdim(dv.v_ / s).pow<-1>();
      }

      /// Add dimensioned values.
      statdim operator+(statdim dv) const { return v_ + dv.v_; }

      /// Subtract dimensioned values.
      statdim operator-(statdim dv) const { return v_ - dv.v_; }

      /// Additive assignment.
      statdim &operator+=(statdim dv)
      {
         v_ += dv.v_;
         return *this;
      }

      /// Subtractive assignment.
      statdim &operator-=(statdim dv)
      {
         v_ -= dv.v_;
         return *this;
      }

      /// Less-than comparison.
      bool operator<(statdim dv) const { return v_ < dv.v_; }

      /// Greater-than comparison.
      bool operator>(statdim dv) const { return v_ > dv.v_; }

      /// Less-than-or-equal-to comparison.
      bool operator<=(statdim dv) const { return v_ <= dv.v_; }

      /// Greater-than-or-equal-to comparison.
      bool operator>=(statdim dv) const { return v_ >= dv.v_; }

      /// Equality comparison.
      bool operator==(statdim dv) const { return v_ == dv.v_; }

      /// Inequality comparison.
      bool operator!=(statdim dv) const { return v_ != dv.v_; }

      /// Type of integer power of present instance.
      /// \tparam E  Integer exponent indicating power.
      template <int P>
      using statdim_power = statdim<TI * P, D * P, M * P, C * P, TE * P>;

      /// Integer power.
      template <int P>
      statdim_power<P> pow() const
      {
         return statdim_power<P>(std::pow(v_, P));
      }

      /// Integer power.
      template <int P>
      friend statdim_power<P> pow(statdim dv)
      {
         return dv.pow<P>();
      }

      /// Type of integer root of present instance.
      /// \tparam E  Integer indicating degree of root.
      template <int R>
      using statdim_root = statdim<TI / R, D / R, M / R, C / R, TE / R>;

      /// Integer root.
      template <int R>
      statdim_root<R> root() const
      {
         static_assert(R > 0, "zero or negative root");
         static_assert(D / R * R == D, "illegal root along distance");
         static_assert(M / R * R == M, "illegal root along mass");
         static_assert(C / R * R == C, "illegal root along charge");
         static_assert(TI / R * R == TI, "illegal root along time");
         static_assert(TE / R * R == TE, "illegal root along temperature");
         double constexpr p = 1.0 / R;
         return statdim_root<R>(std::pow(v_, p));
      }

      /// Integer root.
      template <int R>
      friend statdim_root<R> root(statdim dv)
      {
         return dv.root<R>();
      }

      /// Square root.
      friend statdim_root<2> sqrt(statdim dv) { return dv.root<2>(); }

      /// Square root.
      statdim_root<2> sqrt() const { return root<2>(); }

      /// Absolute value.
      friend statdim fabs(statdim dv) { return fabs(dv.v_); }
   };

   // Definition of static member.
   template <char TI, char D, char M, char C, char TE>
   dim_exps const statdim<TI, D, M, C, TE>::exps_(TI, D, M, C, TE);

   /// Model of a dynamically dimensioned value.  For a statically dimensioned
   /// value, see statdim.
   class dyndim : public dimval<dyndim>
   {
      template <char TI, char D, char M, char C, char TE>
      friend statdim<TI, D, M, C, TE> &operator+=(statdim<TI, D, M, C, TE> &,
                                                  dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend statdim<TI, D, M, C, TE> &operator-=(statdim<TI, D, M, C, TE> &,
                                                  dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator==(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator!=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator<(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator<=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator>(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      template <char TI, char D, char M, char C, char TE>
      friend bool operator>=(statdim<TI, D, M, C, TE> dv, dyndim const &dd);

      using PT = dimval<dyndim>; ///< Type of parent.
      using PT::PT;              ///< Inherit constructor.
      dim_exps exps_;            ///< Storage for dimensional exponents.

      /// Initialize numeric value and dimensional exponents.
      dyndim(double v,  ///< Numeric value corresponding to MKS units.
             dim_exps e ///< Dimensional exponents.
             )
         : PT(v), exps_(e)
      {
      }

   public:
      /// Initialize from either statdim or another dyndim.
      template <typename DER>
      dyndim(dimval<DER> const &dvb)
         : PT(dvb), exps_(dvb.exps())
      {
      }

      dim_exps exps() const { return exps_; } ///< Dimensional exponents.

      /// Unary negation.
      friend dyndim operator-(dyndim const &dv)
      {
         return dyndim(-dv.v_, dv.exps_);
      }

      /// Multiplication by number on right side.
      dyndim operator*(double s) const { return dyndim(v_ * s, exps_); }

      /// Multiplication of dyndim by number on left side.
      friend dyndim operator*(double s, dyndim dv)
      {
         return dyndim(s * dv.v_, dv.exps_);
      }

      /// Multiply dimensioned values.
      template <typename DER>
      dyndim operator*(dimval<DER> const &dvb) const
      {
         return dyndim(v_ * dvb.v_, exps_ + dvb.exps());
      }

      /// Multiply statdim and dyndim.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator*(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
      {
         return dyndim(dv.v_ * dd.v_, dv.exps_ + dd.exps_);
      }

      /// Division by number.
      dyndim operator/(double s) const { return dyndim(v_ / s, exps_); }

      /// Divide dimensioned values.
      template <typename DER>
      dyndim operator/(dimval<DER> const &dvb) const
      {
         return dyndim(v_ / dvb.v_, exps_ - dvb.exps());
      }

      /// Divide number by dimensioned value.
      friend dyndim operator/(double s, dyndim const &dv)
      {
         return dyndim(s / dv.v_, dv.exps_.neg());
      }

      /// Divide statdim by dyndim.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator/(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
      {
         return dyndim(dv.v_ / dd.v_, dv.exps_ - dd.exps_);
      }

      /// Add dimensioned values.
      template <typename DER>
      dyndim operator+(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Addition requires same dimension.";
         }
         return dyndim(v_ + dv.v_, exps_);
      }

      /// Subtract dimensioned values.
      template <typename DER>
      dyndim operator-(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Subtraction requires same dimension.";
         }
         return dyndim(v_ - dv.v_, exps_);
      }

      /// Additive assignment.
      template <typename DER>
      dyndim &operator+=(dimval<DER> const &dv)
      {
         if (exps_ != dv.exps()) {
            throw "Addition requires same dimension.";
         }
         v_ += dv.v_;
         return *this;
      }

      /// Subtractive assignment.
      template <typename DER>
      dyndim &operator-=(dimval<DER> const &dv)
      {
         if (exps_ != dv.exps()) {
            throw "Subtraction requires same dimension.";
         }
         v_ -= dv.v_;
         return *this;
      }

      /// Less-than comparison.
      template <typename DER>
      bool operator<(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ < dv.v_;
      }

      /// Greater-than comparison.
      template <typename DER>
      bool operator>(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ > dv.v_;
      }

      /// Less-than-or-equal-to comparison.
      template <typename DER>
      bool operator<=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ <= dv.v_;
      }

      /// Greater-than-or-equal-to comparison.
      template <typename DER>
      bool operator>=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ >= dv.v_;
      }

      /// Equality comparison.
      template <typename DER>
      bool operator==(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ == dv.v_;
      }

      /// Inequality comparison.
      template <typename DER>
      bool operator!=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            throw "Comparison requires same dimension.";
         }
         return v_ != dv.v_;
      }

      /// Integer power.
      template <int P>
      dyndim pow() const
      {
         return dyndim(std::pow(v_, P), exps_ * P);
      }

      /// Integer power.
      dyndim pow(int p) const { return dyndim(std::pow(v_, p), exps_ * p); }

      /// Integer power.
      template <int P>
      friend dyndim pow(dyndim const &dv)
      {
         return dv.pow<P>();
      }

      /// Integer power.
      template <typename DER>
      friend dyndim pow(dimval<DER> const &dv, int p)
      {
         return dyndim(std::pow(dv.v_, p), dv.exps_ * p);
      }

      /// Integer root.
      template <int R>
      dyndim root() const
      {
         double constexpr p = 1.0 / R;
         return dyndim(std::pow(v_, p), exps_ / R);
      }

      /// Integer root.
      dyndim root(int r) const
      {
         return dyndim(std::pow(v_, 1.0 / r), exps_ / r);
      }

      /// Integer root.
      template <int R>
      friend dyndim root(dyndim const &dv)
      {
         return dv.root<R>();
      }

      /// Integer root.
      template <typename DER>
      friend dyndim root(dimval<DER> const &dv, int r)
      {
         return dyndim(std::pow(dv.v_, 1.0 / r), dv.exps_ / r);
      }

      /// Square root.
      friend dyndim sqrt(dyndim const &dv) { return dv.root<2>(); }

      /// Square root.
      dyndim sqrt() const { return root<2>(); }

      /// Absolute value.
      friend dyndim fabs(dyndim const &dv)
      {
         return dyndim(fabs(dv.v_), dv.exps_);
      }

      /// Convert to double, but throw exception if not dimensionless.
      double number() const
      {
         if (exps_.n() != 0) {
            throw "dyndim converts to number only if dimensionless.";
         }
         return v_;
      }
   };

   /// Additive assignment of dyndim to statdim.
   template <char TI, char D, char M, char C, char TE>
   statdim<TI, D, M, C, TE> &operator+=(statdim<TI, D, M, C, TE> &dv,
                                        dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for addition.";
      }
      dv.v_ += dd.v_;
      return dv;
   }

   /// Subtractive assignment of dyndim from statdim.
   template <char TI, char D, char M, char C, char TE>
   statdim<TI, D, M, C, TE> &operator-=(statdim<TI, D, M, C, TE> &dv,
                                        dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for subtraction.";
      }
      dv.v_ -= dd.v_;
      return dv;
   }

   /// Return true only if statdim be equal to dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator==(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ == dd.v_;
   }

   /// Return true only if statdim be unequal to dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator!=(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ != dd.v_;
   }

   /// Return true only if statdim be less than dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator<(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ < dd.v_;
   }

   /// Return true only if statdim be less than or equal to dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator<=(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ <= dd.v_;
   }

   /// Return true only if statdim be greater than dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator>(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ > dd.v_;
   }

   /// Return true only if statdim be greater than or equal to dyndim.
   template <char TI, char D, char M, char C, char TE>
   bool operator>=(statdim<TI, D, M, C, TE> dv, dyndim const &dd)
   {
      if (dv.exps_ != dd.exps_) {
         throw "Dimensions must be same for comparison.";
      }
      return dv.v_ >= dd.v_;
   }
}

#endif // ndef NUMERIC_DIMVAL_HPP

