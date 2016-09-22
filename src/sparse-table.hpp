
/// \file   sparse-table.hpp
/// \brief  Definition of num::sparse-table.

#ifndef NUMERIC_SPARSE_TABLE_HPP
#define NUMERIC_SPARSE_TABLE_HPP

#include <vector> // for vector

namespace num
{
   /// A [piecewise function](https://en.wikipedia.org/wiki/Piecewise) that, in
   /// logarithmic time, looks up the appropriate sub-function.
   ///
   /// The argument \f$ a \f$ to the function is of type \a A.
   ///
   /// For any \f$ a \f$, what is found in logarithmic time are the
   /// corresponding sub-function and its sub-domain.  The corresponding
   /// sub-function, stored in a table, is \f$ f_i \f$ of type \a F, and the
   /// center of the corresponding subdomain is \f$ a_i \f$, where \f$ i \in
   /// \{0, 1, \ldots, n-1\} \f$.  The sub-argument \f$ a - a_i \f$ is passed
   /// to \f$ f_i \f$, and the value returned by \f$ f_i \f$ is returned from
   /// the lookup.
   ///
   /// In the simplest, fastest table, each \f$ f_i \f$ ignores its argument
   /// and returns a constant value.  However, the infrastructure provided by
   /// sparse_table allows for logarithmic-time lookup of any kind of
   /// sub-function (for example, a cubic interpolant).
   ///
   /// To enable logarithmic-time lookup, an instance of sparse_table stores a
   /// table of \f$ n \f$ triplets \f$ (a_0, {\Delta a}_0, f_0), (a_1, {\Delta
   /// a}_1, f_1), \ldots, (a_{n-1}, {\Delta a}_{n-1}, f_{n-1}) \f$, sorted so
   /// that \f$ a_0 < a_1 < \cdots < a_{n-1} \f$.
   ///
   /// The lookup occurs by way of operator()().  When \f$ a > a_0 \f$ and \f$
   /// a < a_{n-1} \f$, record \f$ i \f$ is found so that
   /// \f[
   ///    -\frac{{\Delta a}_i}{2} < a - a_i \leq +\frac{{\Delta a}_i}{2}.
   /// \f]
   /// What operator()() returns is \f$ f_i(a - a_i) \f$.
   ///
   /// See dense_table for a type that usually can approximate a given function
   /// so precisely as sparse_table but with faster lookup. The cost is that
   /// dense_table requires more storage for the same precision with the same
   /// sub-function type \a F.
   ///
   /// \tparam A  Type of function's argument.
   /// \tparam F  Type of each sub-function \f$ f_i \f$ in the table.
   template <typename A, typename F>
   class sparse_table
   {
      using I = decltype(1.0 / A());
      A a_frst_; ///< First tabulated argument.
      A da_;     ///< Common difference between subsequent tabulated arguments.
      I ida_;    ///< Inverse of difference between subsequent tabulated args.
      std::vector<F> f_; ///< Table of sub-functions.

   public:
      /// Initialize members from a list of arguments.
      sparse_table(A const &first, ///< First tabulated argument \f$ a_0 \f$.
                  A const &delta, ///< Common difference \f$ \Delta a \f$.
                  /// Sub-function objects \f$ f_0, \ldots, f_{n-1} \f$.
                  std::vector<F> vf)
            : a_frst_(first), da_(delta), ida_(1.0 / da_), f_(vf)
      {
      }

      /// First tabulated argument \f$ a_0 \f$.
      A const &a_frst() const { return a_frst_; }

      /// Last tabulated argument \f$ a_{n-1} \f$.
      A a_last() const { return a_frst_ + (f_.size() - 1) * da_; }

      /// Common difference \f$ \Delta a \f$ between subsequent tabulated
      /// arguments.
      A const &da() const { return da_; }

      /// List of sub-functions \f$ \{f_0, f_1, \ldots, f_{n-1}\} \f$.
      std::vector<F> const &f() const { return f_; }

      /// For \f$ f_i \f$ and \f$ a_i \f$ associated with an argument \f$ a \f$
      /// to the function, calculate \f$ f_i(a - a_i) \f$.
      ///
      /// - If \f$ a < a_0 \f$, then choose \f$ i = 0 \f$.
      ///
      /// - If \f$ a > a_{n-1} \f$, then choose \f$ i = n - 1 \f$.
      ///
      /// - Otherwise, choose \f$ i \f$ to be the integer nearest to \f$
      ///   \frac{a - a_0}{\Delta a} + 0.5 \f$.
      ///
      /// \return \f$ f_i(a - a_i) \f$.
      auto operator()(A const &a /**< Argument to function. */) const
            -> decltype(F()(a))
      {
         int i;
         A ai;
         if (a < a_frst_) {
            i = 0;
            ai = a_frst_;
         } else {
            A const last = a_last();
            if (a > last) {
               i = f_.size() - 1;
               ai = last;
            } else {
               i = int((a - a_frst_) * ida_ + 0.5);
               ai = a_frst_ + i * da_;
            }
         }
         return f_[i](a - ai);
      }
   };
}

#endif // ndef NUMERIC_SPARSE_TABLE_HPP

