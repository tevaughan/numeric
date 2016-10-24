
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   dimval.hpp
///
/// \brief  Definition of conversions dealing with angle, num::dimval_base,
///         num::statdim, and num::dyndim.

#ifndef NUMERIC_DIMVAL_HPP
#define NUMERIC_DIMVAL_HPP

#include <cassert>    // for assert()
#include <cmath>      // for M_PI
#include <functional> // for function
#include <iostream>   // for ostream

#include <ginac/ginac.h> // for ex, symbol

#include <dim-exps.hpp> // for dim_exps
#include <integral.hpp> // for integral()
#include <util.hpp>     // for RAT and PRD

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

   class dyndim;

   template <typename DER>
   class dimval;

   template <typename DER>
   dyndim pow(dimval<DER> const &dv, int p);

   /// Ultimate base class for dimensioned quantity.
   struct dimval_base
   {
      static GiNaC::symbol s;  ///< Symbol for second.
      static GiNaC::symbol m;  ///< Symbol for meter.
      static GiNaC::symbol kg; ///< Symbol for kilogram.
      static GiNaC::symbol C;  ///< Symbol for Coulomb.
      static GiNaC::symbol K;  ///< Symbol for Kelvin.

      /// Extract exponents from expression, which is treated as a polynomial
      /// over the MKS unit symbols.
      static dim_exps exps(/** Expression. */ GiNaC::ex const &e)
      {
         dim_exps r; // return value
         r.TI() = e.degree(s);
         r.D()  = e.degree(m);
         r.M()  = e.degree(kg);
         r.C()  = e.degree(C);
         r.TE() = e.degree(K);
         return r;
      }
   };

   /// Base class for statdim and dyndim.
   ///
   /// \tparam DER  Type of descendant, either statdim or dyndim.  This is the
   ///              curiously recurring template pattern (CRTP).
   template <typename DER>
   class dimval : public dimval_base
   {
      /// Needed to allow statdim to access v_ in some circumstances.
      friend class dyndim;

      /// Needed to allow copying from one kind of dimval to another.
      template <typename ODER>
      friend class dimval;

      /// Dynamic integer power of any dimval must be a friend because it needs
      /// to construct a dyndim.
      template <typename ODER>
      friend dyndim pow(dimval<ODER> const &dv, int p);

      /// Dynamic integer root of any dimval must be a friend because it needs
      /// to construct a dyndim.
      template <typename ODER>
      friend dyndim root(dimval<ODER> const &dv, int r);

   protected:
      double v_; ///< Value in MKS.

      /// Construct from double that is known to contain value in MKS.
      dimval(double vv /**< Numeric coefficient of MKS unit. */) : v_(vv) {}

   public:
      /// Return reference to present instance as instance of DER.
      DER &d() { return *static_cast<DER *>(this); }

      /// Return reference to present instance as instance of DER const.
      DER const &d() const { return *static_cast<DER const *>(this); }

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

      /// Convert dimensioned value to ginac expression.
      operator GiNaC::ex() const
      {
         dim_exps const e = exps(); // exponents
         GiNaC::ex      r = v_;     // return value
         r *= pow(kg, e.M());
         r *= pow(m, e.D());
         r *= pow(s, e.TI());
         r *= pow(C, e.C());
         r *= pow(K, e.TE());
         return r;
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

   template <char TI, char D, char M, char C, char TE>
   dyndim operator*(statdim<TI, D, M, C, TE> sd, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   dyndim operator*(dyndim const &dd, statdim<TI, D, M, C, TE> sd);

   template <char TI, char D, char M, char C, char TE>
   dyndim operator/(statdim<TI, D, M, C, TE> sd, dyndim const &dd);

   template <char TI, char D, char M, char C, char TE>
   dyndim operator/(dyndim const &dd, statdim<TI, D, M, C, TE> sd);

   /// Structure used to determine type of integer power of \ref statdim.
   /// \tparam P  Integer exponent indicating power.
   template <int P, char TI, char D, char M, char C, char TE>
   struct statdim_power {
      /// Type of integer power of statdim.
      using type = statdim<TI * P, D * P, M * P, C * P, TE * P>;
   };

   /// Specialization for integer power zero of statdim.
   template <char TI, char D, char M, char C, char TE>
   struct statdim_power<0, TI, D, M, C, TE> {
      /// Type of zeroth power of statdim.
      using type = double;
   };

   // Forward declaration needed to allow a class to declare rk_quad as a
   // friend.
   template <typename X, typename Y>
   class rk_quad;

   template <typename T>
   class tiny;

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

      /// Allow rk_quad to construct from known MKS quantity.
      template <typename X, typename Y>
      friend class rk_quad;

      /// Allow tiny to construct from known MKS quantity.
      template <typename T>
      friend class tiny;

      /// Type of std::function that can be integrated.  A single-argument
      /// function is required.
      template <typename R, typename A>
      using func = std::function<R(A)>;

      /// Allow \ref sparse_table to call constructor.
      template <typename X>
      friend class sparse_table;

      /// Allow \ref dense_table to call constructor.
      template <typename A, typename F>
      friend class dense_table;

      /// Allow integral() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend PRD<R, A> integral(func<R, A> f, A1 a, A2 b, double t,
                                unsigned n);

      /// Allow \ref dyndim to access protected, inherited \a v_ for copy
      /// construction.
      friend class dyndim;

      /// Global multiplication operator against dyndim needs access to v_ and
      /// eps_.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend dyndim operator*(statdim<OTI, OD, OM, OC, OTE> sd,
                              dyndim const &dd);

      /// Global multiplication operator against dyndim needs access to v_ and
      /// eps_.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend dyndim operator*(dyndim const &dd,
                              statdim<OTI, OD, OM, OC, OTE> sd);

      /// Global division operator against dyndim needs access to v_ and eps_.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend dyndim operator/(statdim<OTI, OD, OM, OC, OTE> sd,
                              dyndim const &dd);

      /// Global division operator against dyndim needs access to v_ and eps_.
      template <char OTI, char OD, char OM, char OC, char OTE>
      friend dyndim operator/(dyndim const &dd,
                              statdim<OTI, OD, OM, OC, OTE> sd);

      /// Dimensional exponents in a form useful for interaction with dyndim.
      static dim_exps const exps_;

   protected:
      using PT = dimval<statdim>; ///< Type of parent.
      using PT::PT;               ///< Inherit parent's constructor.
      using PT::v_;               ///< Inherit protected, numeric value.

   public:
      using PT::operator GiNaC::ex; ///< Inherit conversion to expression.

      /// By default, construct a zero-valued quantity (by way of parent's
      /// default constructor).
      statdim() = default;

      /// Construct statdim by copying value in dyndim if it have right
      /// dimension; otherwise throw error.
      statdim(dyndim const &dd);

      /// Use default copy construction for statdim of same dimension.
      statdim(statdim const &sd) = default;

      /// Convert from expression.
      explicit statdim(GiNaC::ex const &e)
      {
         using namespace GiNaC;
         static statdim const unit(1.0);
         ex const             n = e / unit;
         if (!is_a<numeric>(n)) {
            throw "expression not numeric";
         }
         v_ = ex_to<numeric>(n).to_double();
      }

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

      /// Multiplication by number on right side.
      statdim operator*(int s) const { return v_ * s; }

      /// Multiplication of statdim by number on left side.
      friend statdim operator*(double s, statdim dv) { return s * dv.v_; }

      /// Multiplication of statdim by number on left side.
      friend statdim operator*(int s, statdim dv) { return s * dv.v_; }

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

      /// Type of integer power of statdim.
      /// \tparam P  Integer exponent representing power.
      template <int P>
      using power_type = typename statdim_power<P, TI, D, M, C, TE>::type;

      /// Integer power.
      template <int P>
      power_type<P> pow() const
      {
         return power_type<P>(itpow<P>(v_));
      }

      /// Integer power.
      template <int P>
      friend power_type<P> pow(statdim dv)
      {
         return dv.pow<P>();
      }

      /// Type of integer root of present instance.
      /// \tparam E  Integer indicating degree of root.
      template <int R>
      using statdim_root = statdim<TI / R, D / R, M / R, C / R, TE / R>;

      /// Integer root.
      template <int R>
      friend statdim_root<R> root(statdim sd)
      {
         return sd.root<R>();
      }

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

   /// Integer root of dyndim.
   template <int R>
   dyndim root(dyndim const &dv);

   /// Model of a dynamically dimensioned value.  For a statically dimensioned
   /// value, see statdim.
   class dyndim : public dimval<dyndim>
   {
      template <char TI, char D, char M, char C, char TE>
      friend class statdim;

      /// Allow integral_stats to construct from known MKS quantity.
      template <typename I>
      friend class integral_stats;

      /// Allow tiny to construct from known MKS quantity.
      template <typename T>
      friend class tiny;

      /// Allow \ref sparse_table to call constructor.
      template <typename X>
      friend class sparse_table;

      /// Allow \ref dense_table to call constructor.
      template <typename A, typename F>
      friend class dense_table;

      /// Type of std::function that can be integrated.  A single-argument
      /// function is required.
      template <typename R, typename A>
      using func = std::function<R(A)>;

      /// Allow integral() to construct from known MKS quantity.
      template <typename R, typename A, typename A1, typename A2>
      friend PRD<R, A> integral(func<R, A> f, A1 a, A2 b, double t,
                                unsigned n);

      /// Global multiplication operator against statdim needs access to v_ and
      /// exps_.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator*(statdim<TI, D, M, C, TE> sd, dyndim const &dd);

      /// Global multiplication operator against statdim needs access to v_ and
      /// exps_.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator*(dyndim const &dd, statdim<TI, D, M, C, TE> sd);

      /// Global division operator against statdim needs access to v_ and
      /// exps_.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator/(statdim<TI, D, M, C, TE> sd, dyndim const &dd);

      /// Global division operator against statdim needs access to v_ and
      /// exps_.
      template <char TI, char D, char M, char C, char TE>
      friend dyndim operator/(dyndim const &dd, statdim<TI, D, M, C, TE> sd);

      /// Integer power.
      template <typename DER>
      friend dyndim pow(dimval<DER> const &dv, int p);

      /// Integer root.
      template <typename DER>
      friend dyndim root(dimval<DER> const &dv, int r);

      using PT = dimval<dyndim>; ///< Type of parent.
      using PT::PT;              ///< Inherit constructor.
      using PT::v_;              ///< Inherit field.
      dim_exps exps_;            ///< Storage for dimensional exponents.

      /// Initialize numeric value and dimensional exponents.
      dyndim(/** Number. */ double v, /** MKS unit. */ dim_exps e)
         : PT(v), exps_(e)
      {
      }

   public:
      using PT::operator*=;         ///< Inherit multiplicative assignment.
      using PT::operator/=;         ///< Inherit divisive assignment.
      using PT::operator GiNaC::ex; ///< Inherit conversion to expression.
      using dimval_base::exps;      ///< Inherit extraction of units.

      /// By default, construct a dimensionless quantity of magnitude zero.
      dyndim() = default;

      /// Copy from statdim.
      template <char TI, char D, char M, char C, char TE>
      dyndim(/** Object to be copied. */ statdim<TI, D, M, C, TE> const &sd)
         : PT(sd.v_), exps_(sd.exps())
      {
      }

      /// Use default copy constructor.
      dyndim(/** Object to be copied. */ dyndim const& dd) = default;

      /// Copy from expression.
      explicit dyndim(/** Object to be copied. */ GiNaC::ex const &e)
      {
         exps_ = exps(e);
         using namespace GiNaC;
         dyndim const unit(1.0, exps_);
         ex const     n = e / unit;
         if (!is_a<numeric>(n)) {
            throw "expression not numeric";
         }
         v_ = ex_to<numeric>(n).to_double();
      }

      dim_exps exps() const { return exps_; } ///< Dimensional exponents.

      /// Unary negation.
      friend dyndim operator-(dyndim const &dv)
      {
         return dyndim(-dv.v_, dv.exps_);
      }

      /// Multiplication by number on right side.
      dyndim operator*(double s) const { return dyndim(v_ * s, exps_); }

      /// Multiplication by number on right side.
      dyndim operator*(int s) const { return dyndim(v_ * s, exps_); }

      /// Multiplication of dyndim by number on left side.
      friend dyndim operator*(double s, dyndim dv)
      {
         return dyndim(s * dv.v_, dv.exps_);
      }

      /// Multiplication of dyndim by number on left side.
      friend dyndim operator*(int s, dyndim dv)
      {
         return dyndim(s * dv.v_, dv.exps_);
      }

      /// Multiply dyndim by dyndim.
      dyndim operator*(dyndim const &dd) const
      {
         return dyndim(v_ * dd.v_, exps_ + dd.exps_);
      }

      /// Multiplicative assignment by dyndim.
      template <typename DER>
      dyndim& operator*=(dimval<DER> const& dv)
      {
         v_ *= dv.v_;
         exps_ += dv.exps();
         return *this;
      }

      /// Divisive assignment by dyndim.
      template <typename DER>
      dyndim& operator/=(dimval<DER> const& dv)
      {
         v_ /= dv.v_;
         exps_ -= dv.exps();
         return *this;
      }

      /// Division by number.
      dyndim operator/(double s) const { return dyndim(v_ / s, exps_); }

      /// Division by number.
      dyndim operator/(int s) const { return dyndim(v_ / s, exps_); }

      /// Division by number.
      dyndim operator/(unsigned s) const { return dyndim(v_ / s, exps_); }

      /// Divide dyndim by dyndim.
      dyndim operator/(dyndim const &dd) const
      {
         return dyndim(v_ / dd.v_, exps_ - dd.exps_);
      }

      /// Divide number by dimensioned value.
      friend dyndim operator/(double s, dyndim const &dv)
      {
         return dyndim(s / dv.v_, dv.exps_.neg());
      }

      /// Add dimensioned values.
      template <typename DER>
      dyndim operator+(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            assert(0);
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
            assert(0);
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
            assert(0);
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
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ < dv.v_;
      }

      /// Greater-than comparison.
      template <typename DER>
      bool operator>(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ > dv.v_;
      }

      /// Less-than-or-equal-to comparison.
      template <typename DER>
      bool operator<=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            std::cerr << "dyndim::op<=: *this=" << *this << " dv=" << dv
                      << std::endl;
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ <= dv.v_;
      }

      /// Greater-than-or-equal-to comparison.
      template <typename DER>
      bool operator>=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ >= dv.v_;
      }

      /// Equality comparison.
      template <typename DER>
      bool operator==(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ == dv.v_;
      }

      /// Inequality comparison.
      template <typename DER>
      bool operator!=(dimval<DER> const &dv) const
      {
         if (exps_ != dv.exps()) {
            assert(0);
            throw "Comparison requires same dimension.";
         }
         return v_ != dv.v_;
      }

      /// Integer power.
      template <int P>
      dyndim pow() const
      {
         return dyndim(itpow<P>(v_), exps_ * P);
      }

      /// Integer power.
      dyndim pow(int p) const { return dyndim(ipow(v_, p), exps_ * p); }

      /// Integer power.
      template <int P>
      friend dyndim pow(dyndim const &dv)
      {
         return dv.pow<P>();
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

      /// Square root.
      friend dyndim sqrt(dyndim const &dv) { return dv.root<2>(); }

      /// Square root.
      dyndim sqrt() const { return root<2>(); }

      /// Absolute value.
      friend dyndim fabs(dyndim const &dv)
      {
         return dyndim(fabs(dv.v_), dv.exps_);
      }

      /// Convert to double, but throw exception if dimensioned.
      double number() const
      {
         if (exps_.n() != 0) {
            throw "dyndim converts to number only if dimensionless.";
         }
         return v_;
      }

      /// Convert to double, but throw exception if dimensioned.
      operator double() const { return number(); }
   };

   /// Integer power.
   template <typename DER>
   dyndim pow(dimval<DER> const &dv, int p)
   {
      return dyndim(ipow(dv.v_, p), dv.exps() * p);
   }

   /// Integer root.
   template <typename DER>
   dyndim root(dimval<DER> const &dv, int r)
   {
      return dyndim(std::pow(dv.v_, 1.0 / r), dv.exps() / r);
   }

   // Implementation of conversion from dyndim to statdim.
   template <char TI, char D, char M, char C, char TE>
   statdim<TI, D, M, C, TE>::statdim(dyndim const &dd)
   {
      if (exps_ != dd.exps()) {
         throw "Illegal dimension on construction of statdim.";
      }
      v_ = dd.v_;
   }

   /// Multiply a statdim by a dyndim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator*(statdim<TI, D, M, C, TE> sd, dyndim const &dd)
   {
      return dyndim(sd.v_ * dd.v_, sd.exps_ + dd.exps_);
   }

   /// Multiply a dyndim by a statdim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator*(dyndim const &dd, statdim<TI, D, M, C, TE> sd)
   {
      return dyndim(dd.v_ * sd.v_, dd.exps_ + sd.exps_);
   }

   /// Divide a statdim by a dyndim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator/(statdim<TI, D, M, C, TE> sd, dyndim const &dd)
   {
      return dyndim(sd.v_ / dd.v_, sd.exps_ - dd.exps_);
   }

   /// Divide a dyndim by a statdim.
   template <char TI, char D, char M, char C, char TE>
   dyndim operator/(dyndim const &dd, statdim<TI, D, M, C, TE> sd)
   {
      return dyndim(dd.v_ / sd.v_, dd.exps_ - sd.exps_);
   }

   /// Specialization of tiny for dyndim.
   template <>
   class tiny<dyndim>
   {
   public:
      /// Return dyndim with very small value and with units taken from
      /// argument. The number in the passed-in quantity \a u is ignore, but
      /// the units are extracted and included in the return value.
      ///
      /// \tparam D  Intended to by dyndim or statdim.
      /// \return    Very small number attached to units in \a u.
      template <typename D>
      static dyndim
      val(/** Value whose units are extracted. */ dimval<D> const &u)
      {
         return dyndim(1.0E-300, u.exps());
      }
   };
}

#endif // ndef NUMERIC_DIMVAL_HPP

