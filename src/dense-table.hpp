
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   dense-table.hpp
/// \brief  Definition of num::dense-table.

#ifndef NUMERIC_DENSE_TABLE_HPP
#define NUMERIC_DENSE_TABLE_HPP

#include <vector> // for vector

namespace num
{
   /// A [piecewise function](https://en.wikipedia.org/wiki/Piecewise) that, in
   /// constant time, looks up the sub-function appropriate to the argument.
   ///
   /// For any argument \f$ a \f$, the corresponding sub-function and its
   /// sub-domain are found in constant time.  The sub-domains are contiguous
   /// and of common length \f$ \Delta a \f$.  If \f$ a \f$ be less than the
   /// least element of the first subdomain or greater than the greatest
   /// element of the last subdomain, then there is no corresponding
   /// sub-function, and zero is returned.
   ///
   /// Stored in a table of \f$ n \f$ records, the sub-function corresponding
   /// to \f$ a \f$ is \f$ f_i \f$, and the center of the corresponding
   /// subdomain is \f$ a_i \f$, where \f$ i \in \{0, 1, \ldots, n-1\} \f$.
   /// The sub-argument \f$ a - a_i \f$ is passed to \f$ f_i \f$, and the value
   /// returned by \f$ f_i \f$ is returned from the lookup.
   ///
   /// In the simplest, fastest table, each \f$ f_i \f$ ignores its argument
   /// and returns a constant value.  However, the infrastructure provided by
   /// dense_table allows for constant-time lookup of any kind of sub-function
   /// (for example, a cubic interpolant).
   ///
   /// To implement constant-time lookup, an instance of dense_table stores
   ///
   /// - \f$ a_0 \f$, the center of the first sub-domain,
   ///
   /// - the common difference \f$ \Delta a = a_{i+1} - a_i \f$ (for \f$ i \in
   ///   \{0, 1, \ldots, n-2\} \f$) between any center of a sub-domain and the
   ///   next center, and
   ///
   /// - a list of the \f$ n \f$ sub-functions \f$ f_0 \f$, \f$ f_1 \f$,
   ///   \f$ \ldots \f$, \f$ f_{n-1} \f$.
   ///
   /// Thus a regular grid is established across the value of the argument.
   ///
   /// The lookup occurs by way of operator()(). When \f$ a > a_0 -
   /// \frac{\Delta a}{2} \f$, and \f$ a < a_{n-1} + \frac{\Delta a}{2} \f$,
   /// the offset \f$ i \f$ is found as the integer truncation of \f$ \frac{a -
   /// a_0}{\Delta a} + 0.5 \f$.
   ///
   /// operator()() returns \f$ f_i(a - a_i) \f$.
   ///
   /// See sparse_table for a piecewise function that can approximate a
   /// function so precisely as dense_table but with fewer sub-functions of the
   /// same type \a F.  The cost of a smaller table in sparse_table is
   /// logarithmic-time lookup.
   ///
   /// \tparam A  Type of function's argument.
   /// \tparam F  Type of each sub-function \f$ f_i \f$ in the table.
   template <typename A, typename F>
   class dense_table
   {
      using I = decltype(1.0 / A());
      A a_frst_;         ///< Center \f$ a_0 \f$ of first sub-domain.
      A da_;             ///< Common difference between subsequent centers.
      I ida_;            ///< Inverse of common difference.
      std::vector<F> f_; ///< Table of sub-functions.

   public:
      /// Initialize members from a list of arguments.
      dense_table(
            /// Center \f$ a_0 \f$ of first sub-domain.
            A const &first,
            /// Common difference \f$ \Delta a \f$.
            A const &delta,
            /// Sub-function objects.  The first sub-function in \a vf is
            /// interpreted as \f$ f_0 \f$, the second as \f$ f_1 \f$, etc.
            std::vector<F> vf)
         : a_frst_(first), da_(delta), ida_(1.0 / da_), f_(vf)
      {
         if (vf.size() < 1) {
            throw "dense_table must have at least one record.";
         }
         if (delta <= A(0)) {
            throw "Length of sub-domain must be positive.";
         }
      }

      /// Center \f$ a_0 \f$ of first sub-domain.
      A const &a_frst() const { return a_frst_; }

      /// Center \f$ a_{n-1} \f$ of last sub-domain.
      A a_last() const { return a_frst_ + (f_.size() - 1) * da_; }

      /// Common difference \f$ \Delta a \f$ between subsequent centers.
      A const &da() const { return da_; }

      /// List of sub-functions \f$ \{f_0, f_1, \ldots, f_{n-1}\} \f$.
      std::vector<F> const &f() const { return f_; }

      /// Type of value returned by every sub-function.
      using R = decltype(F()(A()));

      /// Find the offset \f$ i \f$ of the sub-domain containing the argument
      /// \f$ a \f$, and return \f$ f_i(a - a_i) \f$.  If \f$ a < a_0 -
      /// \frac{\Delta a}{2} \f$, or \f$ a > a_{n-1} + \frac{\Delta a}{2} \f$,
      /// then return 0.
      ///
      /// \return \f$ f_i(a - a_i) \f$.
      R operator()(A const &a /**< Argument to function. */) const
      {
         if (a < a_frst() - 0.5 * da_ || a > a_last() + 0.5 * da_) {
            return R(0);
         }
         int const i = (a - a_frst()) * ida_ + 0.5;
         A const ai = a_frst() + i * da_;
         return f_[i](a - ai);
      }
   };
}

#endif // ndef NUMERIC_DENSE_TABLE_HPP

