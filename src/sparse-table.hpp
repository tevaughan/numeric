
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   sparse-table.hpp
/// \brief  Definition of num::sparse-table.

#ifndef NUMERIC_SPARSE_TABLE_HPP
#define NUMERIC_SPARSE_TABLE_HPP

#include <algorithm> // for upper_bound()
#include <iostream>  // for cerr, endl
#include <vector>    // for vector

namespace num
{
   template <typename A, typename F>
   class sparse_table;

   template <typename D>
   class dimval;

   template <typename A, typename F, typename D>
   sparse_table<A, decltype(F() * D())>
   operator*(sparse_table<A, F> const &tab, dimval<D> const &fac);

   template <typename A, typename F, typename D>
   sparse_table<A, decltype(D() * F())>
   operator*(dimval<D> const &fac, sparse_table<A, F> const &tab);

   template <typename A, typename F, typename D>
   sparse_table<A, decltype(F() / D())>
   operator/(sparse_table<A, F> const &tab, dimval<D> const &fac);

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
      /// Allow other type of sparse_table to access private members.
      /// \tparam OA  Type of other sparse_table's argument.
      /// \tparam OF  Type of other sparse_table's subfunction.
      template <typename OA, typename OF>
      friend class sparse_table;

      /// Allow for external operator to use private constructor.
      template <typename OA, typename OF, typename D>
      friend sparse_table<OA, decltype(OF() * D())>
      operator*(sparse_table<OA, OF> const &tab, dimval<D> const &fac);

      /// Allow for external operator to use private constructor.
      template <typename OA, typename OF, typename D>
      friend sparse_table<OA, decltype(D() * OF())>
      operator*(dimval<D> const &fac, sparse_table<OA, OF> const &tab);

      /// Allow for external operator to use private constructor.
      template <typename OA, typename OF, typename D>
      friend sparse_table<OA, decltype(OF() / D())>
      operator/(sparse_table<OA, OF> const &tab, dimval<D> const &fac);

   public:
      /// Type of record in table.
      struct rec
      {
         A a;  ///< Center of sub-domain.
         A da; ///< Length of sub-domain.
         F f;  ///< Sub-function.
      };

      using data = std::vector<rec>; ///< Type of data structure for table.

   private:
      data dat_; ///< Tabular data.

      /// Compare argument with record so that right record can be found.
      static bool acomp(A const &a, rec const &r) { return a < r.a; }

      /// Compare two records so that they can be sorted.
      static bool rcomp(rec const &r1, rec const &r2) { return r1.a < r2.a; }

      /// Construct from data member, already consistently initialized and
      /// sorted.
      explicit sparse_table(data &&d) : dat_(std::move(d)) {}

   public:
      /// Construct null table.
      sparse_table() {}

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
         std::sort(dat_.begin(), dat_.end(), rcomp);
         if (vf[0].first <= 0.0 * vf[0].first) {
            throw "Length of sub-domain must be positive.";
         }
         dat_[0].a  = a0;
         dat_[0].da = vf[0].first;
         dat_[0].f  = vf[0].second;
         for (unsigned i = 1; i < vf.size(); ++i) {
            if (vf[i].first <= 0.0 * vf[0].first) {
               throw "Length of sub-domain must be positive.";
            }
            dat_[i].da = vf[i].first;
            dat_[i].f  = vf[i].second;
            dat_[i].a  = dat_[i - 1].a + 0.5 * (dat_[i - 1].da + dat_[i].da);
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
      R operator()(/** Argument to function. */ A const &a) const
      {
         if (a < dat_.begin()->a - 0.5 * dat_.begin()->da ||
             a > dat_.rbegin()->a + 0.5 * dat_.rbegin()->da) {
            auto q = dat_.begin();
            return 0.0 * q->f(a - q->a);
         }
         // In log time, find pointer to first record after argument a.
         auto p = std::upper_bound(dat_.begin(), dat_.end(), a, acomp);
         if (p == dat_.end()) {
            // Argument a is after last center.
            auto q = dat_.rbegin();
            return q->f(a - q->a);
         } else if (p->a - a > 0.5 * p->da) {
            --p; // Argument a is too far from subsequent center.
         }
         return p->f(a - p->a);
      }

      /// Type of definite integral.
      using integral_type = decltype(A() * R());

      /// Integral of piece-wise function over all pieces.
      integral_type integral() const
      {
         integral_type rv = 0; // Return value.
         for (auto i : dat_) {
            rv += i.f.integral(-0.5 * i.da)(+0.5 * i.da);
         }
         return rv;
      }

      /// Integral of piece-wise function over range.
      integral_type
      integral(/** Beginning of range. */ A a, /** End of range. */ A b) const
      {
         double sign = 1.0;
         if (a > b) {
            sign = -1.0;
            std::swap(a, b);
         }
         integral_type rv = 0; // Return value.
         if (dat_.size() == 0) {
            return rv; // There are no pieces over which to integrate.
         }
         A const beg = dat_.begin()->a - 0.5 * dat_.begin()->da;
         A const end = dat_.rbegin()->a + 0.5 * dat_.rbegin()->da;
         if (a > end) {
            return rv; // Interval is entirely past the last piece.
         }
         if (b < beg) {
            return rv; // Interval is entirely before the first piece.
         }
         if (a < beg) {
            a = beg; // Interval is partially before the first piece.
         }
         if (b > end) {
            b = end; // Interval is partially after the last piece.
         }
         auto pa = std::upper_bound(dat_.begin(), dat_.end(), a, acomp);
         auto pb = std::upper_bound(dat_.begin(), dat_.end(), b, acomp);
         if (pa == dat_.end()) {
            // Beginning of interval is after last center.
            auto q = dat_.rbegin();
            return sign * q->f.integral(a - q->a)(b - q->a);
         } else if (pa->a - a > 0.5 * pa->da) {
            --pa; // Beginning of interval is too far from subsequent center.
         }
         if (pb == dat_.end() || pb->a - b > 0.5 * pb->da) {
            // End of interval is after last center or too far from subsequent
            // center.
            --pb;
         }
         auto const pend = pb + 1;
         for (auto i = pa; i != pend; ++i) {
            A aa; // Beginning of integration relative to center of piece.
            A bb; // End of integration relative to center of piece.
            if (i == pa) {
               aa = a - i->a;     // Beginning of interval relative to center.
            } else {              //
               aa = -0.5 * i->da; // Beginning of piece relative to center.
            }
            if (i == pb) {
               bb = b - i->a;     // End of interval relative to center.
            } else {              //
               bb = +0.5 * i->da; // End of piece relative to center.
            }
            rv += i->f.integral(aa)(bb); // Integral over current piece.
         }
         return sign * rv;
      }

      /// Multiply table by scale factor on right.
      sparse_table operator*(/** Factor. */ double fac) const
      {
         data d(dat_.size()); // Initializer for return value.
         for (unsigned i = 0; i < dat_.size(); ++i) {
            d[i].a  = dat_[i].a;
            d[i].da = dat_[i].da;
            d[i].f  = dat_[i].f * fac;
         }
         return sparse_table(std::move(d));
      }

      /// Multiply table by scale factor on left.
      friend sparse_table operator*(
            /** Factor. */ double fac, /** Table. */ sparse_table const &tab)
      {
         data d(tab.dat_.size()); // Initializer for return value.
         for (unsigned i = 0; i < tab.dat_.size(); ++i) {
            d[i].a  = tab.dat_[i].a;
            d[i].da = tab.dat_[i].da;
            d[i].f  = fac * tab.dat_[i].f;
         }
         return sparse_table(std::move(d));
      }

      /// Multiplicative assignment.
      sparse_table &operator*=(/** Factor. */ double rf)
      {
         for (unsigned i = 0; i < dat_.size(); ++i) {
            dat_[i].f *= rf;
         }
         return *this;
      }

      /// Divisive assignment.
      sparse_table &operator/=(/** Factor. */ double rf)
      {
         for (unsigned i = 0; i < dat_.size(); ++i) {
            dat_[i].f /= rf;
         }
         return *this;
      }
   };
}

#include <dimval.hpp>

namespace num
{
   /// Multiply table by dimval on right.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam F  Type of each functional piece of table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename F, typename D>
   sparse_table<A, decltype(F() * D())> operator*(
         /** Table.  */ sparse_table<A, F> const &tab,
         /** Factor. */ dimval<D> const &fac)
   {
      using PF = decltype(F() * D());
      // Initializer for return value.
      typename sparse_table<A, PF>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = tab.dat()[i].f * fac.d();
      }
      return sparse_table<A, PF>(std::move(d));
   }

   /// Multiply table by dimval on left.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam F  Type of each functional piece of table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename F, typename D>
   sparse_table<A, decltype(D() * F())> operator*(
         /** Factor. */ dimval<D> const &fac,
         /** Table.  */ sparse_table<A, F> const &tab)
   {
      using PF = decltype(D() * F());
      // Initializer for return value.
      typename sparse_table<A, PF>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = fac.d() * tab.dat()[i].f;
      }
      return sparse_table<A, PF>(std::move(d));
   }

   /// Divide table by dimval.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam F  Type of each functional piece of table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename F, typename D>
   sparse_table<A, decltype(F() / D())> operator/(
         /** Table.  */ sparse_table<A, F> const &tab,
         /** Factor. */ dimval<D> const &fac)
   {
      using PF = decltype(F() / D());
      // Initializer for return value.
      typename sparse_table<A, PF>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = tab.dat()[i].f / fac.d();
      }
      return sparse_table<A, PF>(std::move(d));
   }
}

#endif // ndef NUMERIC_SPARSE_TABLE_HPP

