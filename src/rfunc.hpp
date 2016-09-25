
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// FIXME: This should probably be cfunc, for continuous-variable function.

/// \file   rfunc.hpp
/// \brief  Definition of num::rfunc.

#ifndef NUMERIC_RFUNC_HPP
#define NUMERIC_RFUNC_HPP

namespace num
{
   /// Base class for model of a function of a continuous variable.  The
   /// immediate descendant's type is a template type parameter.  See
   /// [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
   ///
   /// The descendant may use virtual functions but is not required to, and so
   /// we may pay the cost of dynamic polymorphism on a per-descendant basis.
   ///
   /// \tparam A  Type of function's argument.
   /// \tparam D  Type of descendant.
   template <typename A, typename D>
   struct rfunc {
      /// Type of derivative.
      using DERIV = decltype(D().derivative());

      /// Type of indefinite integral.
      using INTEG = decltype(D().integral());

      /// Return pointer to descendant instance.
      D const& d() const { return static_cast<D const*>(this); }

      /// Return function representing derivative.
      DERIV derivative() const { return d()->derivative(); }

      /// Return function representing integral from specified lower bound.
      INTEG integral(A const& lb) const { return d()->integral(lb); }
   };
}

#endif // ndef NUMERIC_RFUNC_HPP

