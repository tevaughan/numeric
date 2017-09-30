
// Copyright 2016-2017  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   dim-exps.hpp
/// \brief  Definition of class num::dim_exps.

#ifndef NUMERIC_DIM_EXPS_HPP
#define NUMERIC_DIM_EXPS_HPP

#include <cstdint> // for uint64_t

namespace num
{
   /// Dimensional exponents.
   class dim_exps
   {
      /// Type of storage for exponents.
      union s {
         char e[8];  ///< Dimensional exponents.
         uint64_t n; ///< Dimensional exponents as single number.
      };

      s s_; ///< Storage for exponents.

   public:
      /// By default, set all exponents to zero.
      dim_exps() : s_{{0, 0, 0, 0, 0, 0, 0, 0}} {}

      /// Initialize exponent for each dimension.
      dim_exps(/** time */ char ti, /** dist */ char ds, /** mass */ char ms,
               /** chrg */ char ch, /** temp */ char te)
         : s_{{ti, ds, ms, ch, te, 0, 0, 0}}
      {
      }

      /// Reference to mutable exponent of time.
      char &TI() { return s_.e[0]; }

      /// Reference to mutable exponent of distance.
      char &D() { return s_.e[1]; }

      /// Reference to mutable exponent of mass.
      char &M() { return s_.e[2]; }

      /// Reference to mutable exponent of charge.
      char &C() { return s_.e[3]; }

      /// Reference to mutable exponent of temperature.
      char &TE() { return s_.e[4]; }

      /// Reference to immutable exponent of time.
      char const &TI() const { return s_.e[0]; }

      /// Reference to immutable exponent of distance.
      char const &D() const { return s_.e[1]; }

      /// Reference to immutable exponent of mass.
      char const &M() const { return s_.e[2]; }

      /// Reference to immutable exponent of electric charge.
      char const &C() const { return s_.e[3]; }

      /// Reference to immutable exponent of temperature.
      char const &TE() const { return s_.e[4]; }

      /// Reference to immutable number representing all exponents.
      uint64_t const &n() const { return s_.n; }

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

      /// Additive assignment of exponents for multiplicative assignment of
      /// dimensioned quantity.
      dim_exps &operator+=(dim_exps ode /**< Other exponents. */)
      {
         TI() += ode.TI();
         D() += ode.D();
         M() += ode.M();
         C() += ode.C();
         TE() += ode.TE();
         return *this;
      }

      /// Subtractive assignment of exponents for divisive assignment of
      /// dimensioned quantity.
      dim_exps &operator-=(dim_exps ode /**< Other exponents. */)
      {
         TI() -= ode.TI();
         D() -= ode.D();
         M() -= ode.M();
         C() -= ode.C();
         TE() -= ode.TE();
         return *this;
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
         return s_.n == ode.s_.n;
      }

      /// Return true only if any exponent be different from correspondent in
      /// other set.
      bool operator!=(dim_exps ode /**< Other exponents. */) const
      {
         return s_.n != ode.s_.n;
      }

      /// Return negative exponents as for reciprocal of dimensioned value.
      dim_exps neg() const { return dim_exps(-TI(), -D(), -M(), -C(), -TE()); }
   };
}

#endif // ndef NUMERIC_DIM_EXPS_HPP

