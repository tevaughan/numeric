
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   sparse-table.hpp
/// \brief  Definition of num::sparse-table.

#ifndef NUMERIC_SPARSE_TABLE_HPP
#define NUMERIC_SPARSE_TABLE_HPP

#include <ginac/ginac.h> // for ex

#include <algorithm> // for upper_bound()
#include <vector>    // for vector

namespace num
{
   template <typename A>
   class sparse_table;

   template <typename D>
   class dimval;

   template <typename A, typename D>
   sparse_table<A> operator*(sparse_table<A> const &tab, dimval<D> const &fac);

   template <typename A, typename D>
   sparse_table<A> operator*(dimval<D> const &fac, sparse_table<A> const &tab);

   template <typename A, typename D>
   sparse_table<A> operator/(sparse_table<A> const &tab, dimval<D> const &fac);

   /// Base class for sparse_table.
   struct sparse_table_base
   {
      static GiNaC::symbol x; ///< Symbol for independent variable.
   };

   /// A [piecewise function](https://en.wikipedia.org/wiki/Piecewise) that, in
   /// logarithmic time, looks up the sub-function appropriate to the argument.
   ///
   /// For any argument \f$a\f$, the corresponding sub-function and its
   /// sub-domain are found in logarithmic time.  The sub-domains are
   /// contiguous and of variable length.  If \f$a\f$ be less than the least
   /// element of the first subdomain or greater than the greatest element of
   /// the last subdomain, then there is no corresponding sub-function, and
   /// zero is returned.
   ///
   /// Stored in a table of \f$n\f$ records, the sub-function corresponding to
   /// \f$a\f$ is \f$f_i\f$, and the center of the corresponding subdomain is
   /// \f$a_i\f$, where \f$i \in \{0, 1, \ldots, n-1\}\f$.  The argument
   /// \f$a\f$ is passed to \f$f_i\f$, and the value returned by \f$f_i\f$ is
   /// returned from the lookup.
   ///
   /// In the simplest, fastest table, each \f$f_i\f$ ignores its argument and
   /// returns a constant value.
   ///
   /// The infrastructure provided by sparse_table allows for logarithmic-time
   /// lookup of any kind of sub-function (for example, a cubic interpolant).
   /// To implement logarithmic-time lookup, an instance of sparse_table stores
   /// \f$n\f$ triplets
   /// \f$(a_0, \Delta a_0, f_0)\f$,
   /// \f$(a_1, \Delta a_1, f_1)\f$,
   /// \f$\ldots\f$,
   /// \f$(a_{n-1}, \Delta a_{n-1}, f_{n-1})\f$,
   /// sorted so that \f$a_0 < a_1 < \cdots < a_{n-1}\f$, where \f$\Delta
   /// a_i\f$ is the length of the sub-domain whose center is \f$a_i\f$.
   ///
   /// The lookup occurs by way of operator()().  When
   /// \f$a > a_0 - \frac{\Delta a_0}{2}\f$, and
   /// \f$a < a_{n-1} + \frac{\Delta a_{n-1}}{2}\f$,
   /// record \f$i\f$ is found so that
   /// \f$|a - a_i| \leq \frac{{\Delta a}_i}{2}\f$.
   ///
   /// operator()() returns \f$f_i(a)\f$.
   ///
   /// See dense_table for a type that usually can approximate a given function
   /// so precisely as sparse_table but with faster lookup. The cost is that
   /// dense_table requires more storage for the same precision with the same
   /// sub-function type \a F.
   ///
   /// \tparam A  Type of sub-function's argument.
   template <typename A>
   class sparse_table : public sparse_table_base
   {
      /// Allow other type of sparse_table to access private members.
      /// \tparam OA  Type of other sparse_table's argument.
      template <typename OA>
      friend class sparse_table;

      /// Allow for external operator to use private constructor.
      template <typename OA, typename D>
      friend sparse_table<OA>
      operator*(sparse_table<OA> const &tab, dimval<D> const &fac);

      /// Allow for external operator to use private constructor.
      template <typename OA, typename D>
      friend sparse_table<OA>
      operator*(dimval<D> const &fac, sparse_table<OA> const &tab);

      /// Allow for external operator to use private constructor.
      template <typename OA, typename D>
      friend sparse_table<OA>
      operator/(sparse_table<OA> const &tab, dimval<D> const &fac);

   public:
      /// Type of record in table.
      struct rec
      {
         A         a;  ///< Center of sub-domain.
         A         da; ///< Length of sub-domain.
         GiNaC::ex f;  ///< Sub-function.
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
            std::vector<std::pair<A, GiNaC::ex>> vf)
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
      /// \f$ (a_0, \Delta a_0, f_0) \f$,
      /// \f$ (a_1, \Delta a_1, f_1) \f$,
      /// \f$ \ldots \f$,
      /// \f$ (a_{n-1}, \Delta a_{n-1}, f_{n-1}) \f$.
      data const &dat() const { return dat_; }

      /// Find \f$a_i\f$ whose sub-domain contains \f$a\f$, and return
      /// \f$f_i(a)\f$.  If
      /// \f$a < a_0     - \frac{\Delta a_{0}}{2}\f$, or
      /// \f$a > a_{n-1} + \frac{\Delta a_{n-1}}{2}\f$,
      /// then return 0.
      ///
      /// \return \f$ f_i(a) \f$.
      GiNaC::ex operator()(/** Argument to function. */ A const &a) const
      {
         auto const &frst = *dat_.begin();
         auto const &last = *dat_.rbegin();
         if (a < frst.a - 0.5 * frst.da || a > last.a + 0.5 * last.da) {
            return 0;
         }
         // In log time, find pointer to first record after argument a.
         auto p = std::upper_bound(dat_.begin(), dat_.end(), a, acomp);
         if (p == dat_.end()) {
            return last.f.subs(x == a); // Argument a is after last center.
         } else if (p->a - a > 0.5 * p->da) {
            --p; // Argument a is too far from subsequent center.
         }
         return p->f.subs(x == a);
      }

      /// Integral of piece-wise function over all pieces.
      GiNaC::ex integral() const
      {
         GiNaC::ex rv = 0; // Return value.
         for (auto i : dat_) {
            rv += GiNaC::integral(x, i.a - 0.5 * i.da, i.a + 0.5 * i.da, i.f);
         }
         return rv.evalf();
      }

      /// Integral of piece-wise function over range.
      GiNaC::ex
      integral(/** Beginning of range. */ A a, /** End of range. */ A b) const
      {
         double sign = 1.0;
         if (a > b) {
            sign = -1.0;
            std::swap(a, b);
         }
         GiNaC::ex rv = 0; // Return value.
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
            return (sign * GiNaC::integral(x, a, b, dat_.rbegin()->f)).evalf();
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
            A aa;                        // Current beg of integration.
            A bb;                        // Current end of integration.
            if (i == pa) {               //
               aa = a;                   // Beg of interval.
            } else {                     //
               aa = i->a - 0.5 * i->da;  // Beg of piece.
            }                            //
            if (i == pb) {               //
               bb = b;                   // End of interval.
            } else {                     //
               bb = i->a + 0.5 * i->da;  // End of piece.
            }                            //
            // Integral over current piece.
            rv += GiNaC::integral(x, aa, bb, i->f);
         }
         return (sign * rv).evalf();
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

      /// Multiply table by other table.
      sparse_table
      operator*(sparse_table const &st)
      {
         // Initializer for return value.  The initial size might be too large
         // and is reduced if necessary at the end.
         data d(dat_.size() + st.dat().size());
         unsigned i  = 0; // offset into returned array
         unsigned i1 = 0; // offset into dat_
         unsigned i2 = 0; // offset into st.dat_
         while (i1 < dat_.size() && i2 < st.dat().size()) {
            auto const &d1   = dat_[i1];
            auto const &d2   = st.dat()[i2];
            auto const  hda1 = 0.5 * d1.da;
            auto const  hda2 = 0.5 * d2.da;
            auto const  l1   = d1.a - hda1;
            auto const  l2   = d2.a - hda2;
            auto const  r1   = d1.a + hda1;
            auto const  r2   = d2.a + hda2;
            if (l1 < l2) {
               // d1 starts before d2.
               if (r1 < l2) {
                  // d1 ends before d2 starts.
                  ++i1;
               } else if (r1 < r2) {
                  // d1 ends before d2 ends.
                  d[i].a  = 0.5 * (l2 + r1);
                  d[i].da = r1 - l2;
                  d[i].f  = d1.f * d2.f;
                  ++i1;
                  ++i;
               } else {
                  // d1 ends after d2 ends.
                  d[i].a  = d2.a;
                  d[i].da = d2.da;
                  d[i].f  = d1.f * d2.f;
                  ++i2;
                  ++i;
               }
            } else {
               // d2 starts before d1.
               if (r2 < l1) {
                  // d2 ends before d1 starts.
                  ++i2;
               } else if (r2 < r1) {
                  // d2 ends before d1 ends.
                  d[i].a  = 0.5 * (l1 + r2);
                  d[i].da = r2 - l1;
                  d[i].f  = d1.f * d2.f;
                  ++i2;
                  ++i;
               } else {
                  // d2 ends after d1 ends.
                  d[i].a  = d1.a;
                  d[i].da = d1.da;
                  d[i].f  = d1.f * d2.f;
                  ++i1;
                  ++i;
               }
            }
         }
         d.resize(i); // Shrink down to size if necessary.
         return sparse_table(std::move(d));
      }
   };
}

#include <dimval.hpp>

namespace num
{
   /// Multiply table by dimval on right.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename D>
   sparse_table<A> operator*(
         /** Table.  */ sparse_table<A> const &tab,
         /** Factor. */ dimval<D> const &fac)
   {
      // Initializer for return value.
      typename sparse_table<A>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = tab.dat()[i].f * fac.d();
      }
      return sparse_table<A>(std::move(d));
   }

   /// Multiply table by dimval on left.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename D>
   sparse_table<A> operator*(
         /** Factor. */ dimval<D> const &fac,
         /** Table.  */ sparse_table<A> const &tab)
   {
      // Initializer for return value.
      typename sparse_table<A>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = fac.d() * tab.dat()[i].f;
      }
      return sparse_table<A>(std::move(d));
   }

   /// Divide table by dimval.
   /// \tparam A  Type of argument to function modeled by table.
   /// \tparam D  Derived type of dimval.
   template <typename A, typename D>
   sparse_table<A> operator/(
         /** Table.  */ sparse_table<A> const &tab,
         /** Factor. */ dimval<D> const &fac)
   {
      // Initializer for return value.
      typename sparse_table<A>::data d(tab.dat().size());
      for (unsigned i = 0; i < d.size(); ++i) {
         d[i].a  = tab.dat()[i].a;
         d[i].da = tab.dat()[i].da;
         d[i].f  = tab.dat()[i].f / fac.d();
      }
      return sparse_table<A>(std::move(d));
   }
}

#endif // ndef NUMERIC_SPARSE_TABLE_HPP

