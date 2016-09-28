
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   sparse-table.hpp
/// \brief  Definition of num::sparse-table.

#ifndef NUMERIC_SPARSE_TABLE_HPP
#define NUMERIC_SPARSE_TABLE_HPP

#include <algorithm> // for upper_bound()
#include <vector>    // for vector

namespace num
{
   /// A [piecewise function](https://en.wikipedia.org/wiki/Piecewise) that, in
   /// logarithmic time, looks up the sub-function appropriate to the argument.
   ///
   /// For any argument \f$ a \f$, the corresponding sub-function and its
   /// sub-domain are found in logarithmic time.  The sub-domains are
   /// contiguous and of variable length.  If \f$ a \f$ be less than the least
   /// element of the first subdomain or greater than the greatest element of
   /// the last subdomain, then there is no corresponding sub-function, and
   /// zero is returned.
   ///
   /// Stored in a table of \f$ n \f$ records, the sub-function corresponding
   /// to \f$ a \f$ is \f$ f_i \f$, and the center of the corresponding
   /// subdomain is \f$ a_i \f$, where \f$ i \in \{0, 1, \ldots, n-1\} \f$.
   /// The sub-argument \f$ a - a_i \f$ is passed to \f$ f_i \f$, and the value
   /// returned by \f$ f_i \f$ is returned from the lookup.
   ///
   /// In the simplest, fastest table, each \f$ f_i \f$ ignores its argument
   /// and returns a constant value.  However, the infrastructure provided by
   /// sparse_table allows for logarithmic-time lookup of any kind of
   /// sub-function (for example, a cubic interpolant).
   ///
   /// To implement logarithmic-time lookup, an instance of sparse_table stores
   /// \f$ n \f$ triplets \f$ (a_0, \Delta a_0, f_0) \f$, \f$ (a_1, \Delta a_1,
   /// f_1) \f$, \f$ \ldots \f$, \f$ (a_{n-1}, \Delta a_{n-1}, f_{n-1}) \f$,
   /// sorted so that \f$ a_0 < a_1 < \cdots < a_{n-1} \f$, where \f$ \Delta
   /// a_i \f$ is the length of the sub-domain whose center is \f$ a_i \f$.
   ///
   /// The lookup occurs by way of operator()().  When \f$ a > a_0 -
   /// \frac{\Delta a_0}{2} \f$, and \f$ a < a_{n-1} + \frac{\Delta a_{n-1}}{2}
   /// \f$, record \f$ i \f$ is found so that \f$ |a - a_i| \leq \frac{{\Delta
   /// a}_i}{2} \f$.
   ///
   /// operator()() returns \f$ f_i(a - a_i) \f$.
   ///
   /// See dense_table for a type that usually can approximate a given function
   /// so precisely as sparse_table but with faster lookup. The cost is that
   /// dense_table requires more storage for the same precision with the same
   /// sub-function type \a F.
   ///
   /// \tparam A  Type of sub-function's argument.
   /// \tparam F  Type of sub-function.
   template <typename A, typename F>
   class sparse_table
   {
   public:
      /// Type of record in table.
      struct rec {
         A a;  ///< Center of sub-domain.
         A da; ///< Length of sub-domain.
         F f;  ///< Sub-function.
      };

      using data = std::vector<rec>; ///< Type of data structure for table.

   private:
      data dat_; ///< Tabular data.

      /// Compare argument with record so that right record can be found.
      static bool comp(A const &a, rec const &r) { return a < r.a; }

   public:
      /// Initialize table of sub-domain centers, sub-domain lengths, and
      /// sub-functions.
      sparse_table(
            /// Center \f$ a_0 \f$ of first subdomain.
            A const &a0,
            /// List of pairs, each containing the length of a sub-domain and a
            /// sub-function for that sub-domain.  The first pair in \a vf is
            /// interpreted as \f$ (\Delta a_0, f_0) \f$, the second as \f$
            /// (\Delta a_1, f_1) \f$, etc.
            std::vector<std::pair<A, F>> vf)
         : dat_(vf.size())
      {
         if (vf.size() == 0) {
            throw "sparse_table must have at least one entry.";
         }
         if (vf[0].first <= A(0)) {
            throw "Length of sub-domain must be positive.";
         }
         dat_[0].a = a0;
         dat_[0].da = vf[0].first;
         dat_[0].f = vf[0].second;
         for (unsigned i = 1; i < vf.size(); ++i) {
            if (vf[i].first <= A(0)) {
               throw "Length of sub-domain must be positive.";
            }
            dat_[i].da = vf[i].first;
            dat_[i].f = vf[i].second;
            dat_[i].a = dat_[i - 1].a + 0.5 * (dat_[i - 1].da + dat_[i].da);
         }
      }

      /// Table of sub-domain centers, sub-domain lengths, and sub-functions:
      /// \f$ (a_0, \Delta a_0, f_0) \f$, \f$ (a_1, \Delta a_1, f_1) \f$, \f$
      /// \ldots \f$, \f$ (a_{n-1}, \Delta a_{n-1}, f_{n-1}) \f$.
      data const &dat() const { return dat_; }

      /// Type of value returned by every sub-function.
      using R = decltype(F()(A()));

      /// Find \f$ a_i \f$ whose sub-domain contains \f$ a \f$, and return \f$
      /// f_i(a - a_i) \f$.  If \f$ a < a_0 - \frac{\Delta a_0}{2} \f$, or \f$
      /// a > a_{n-1} + \frac{\Delta a_{n-1}}{2} \f$, then return 0.
      ///
      /// \return \f$ f_i(a - a_i) \f$.
      R operator()(A const &a /**< Argument to function. */) const
      {
         if (a < dat_.begin()->a - 0.5 * dat_.begin()->da ||
             a > dat_.rbegin()->a + 0.5 * dat_.rbegin()->da) {
            return R(0);
         }
         // In log time, find pointer to first record after argument a.
         auto p = std::upper_bound(dat_.begin(), dat_.end(), a, comp);
         if (p == dat_.end()) {
            // Argument a is after last center.
            auto q = dat_.rbegin();
            return q->f(a - q->a);
         } else if (p->a - a > 0.5 * p->da) {
            --p; // Argument a is too far from subsequent center.
         }
         return p->f(a - p->a);
      }
   };
}

#endif // ndef NUMERIC_SPARSE_TABLE_HPP

