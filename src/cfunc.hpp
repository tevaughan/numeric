
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   cfunc.hpp
/// \brief  Definition of num::cfunc.

#ifndef NUMERIC_CFUNC_HPP
#define NUMERIC_CFUNC_HPP

namespace num
{
   /// Base class for model of a function of a continuous variable.  The
   /// immediate descendant's type \a D is a template type parameter.  See
   /// [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
   ///
   /// The descendant may use virtual functions for its descendants but is not
   /// required to, and so we may pay the cost of dynamic polymorphism on a
   /// per-descendant basis.
   ///
   /// \tparam A  Type of function's argument.
   /// \tparam R  Type returned by function.
   /// \tparam D  Type of descendant.
   template <typename A, typename R, typename D>
   struct cfunc {
      /// Type of derivative.
      using DERIV = decltype(D().derivative());

      /// Type of indefinite integral.
      using INTEG = decltype(D().integral());

      /// Return pointer to descendant instance.
      D const &d() const { return static_cast<D const *>(this); }

      /// Evaluate function.
      R operator()(A const& a) const { return d()(a); }

      /// Return function representing derivative.
      DERIV derivative() const { return d()->derivative(); }

      /// Return function representing integral from specified lower bound.
      INTEG integral(A const &lb) const { return d()->integral(lb); }
   };
}

#endif // ndef NUMERIC_CFUNC_HPP

