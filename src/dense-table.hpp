
/// \file   dense-table.hpp
/// \brief  Definition of num::dense-table.

#ifndef NUMERIC_DENSE_TABLE_HPP
#define NUMERIC_DENSE_TABLE_HPP

#include <vector> // for vector

namespace num
{
   /// Simple, constant-time-lookup model of a function of a continuous
   /// variable.
   ///
   /// To enable constant-time lookup, dense_table contains a regular grid
   /// across the value of the function's argument.  So each pair of subsequent
   /// table entries is separated by the same distance \f$ \Delta a \f$ in the
   /// space of the function's argument.
   ///
   /// At the \f$ i^{\text{th}} \f$ entry, corresponding to the value \f$ a_i
   /// \f$ of the function's argument, an instance of dense_table stores a
   /// function object \f$ f_i \f$ of type \a F.
   ///
   /// When the lookup occurs, by way of dense_table::operator(A const&), let
   /// the value of the passed-in argument be \f$ a \f$.  When \f$ a \f$ is
   /// neither before the first table entry nor after the last, the index \f$ i
   /// \f$ is chosen so that \f[ -\frac{\Delta a}{2} < a - a_i \leq
   /// +\frac{\Delta a}{2}. \f] What is returned by dense_table::operator() is
   /// the value \f$ f_i(a - a_i) \f$.
   ///
   /// For log-time-lookup, see sparse_table.
   ///
   /// \tparam A  Type of function's argument.
   ///
   /// \tparam F  Type of each function \f$ f_i \f$, of which an instance of
   ///            dense_table is a piece-wise construction.
   template<typename A, typename F>
   class dense_table
   {
      using I = decltype(1.0 / A());
      A a_frst_; ///< First tabulated argument.
      A da_;     ///< Difference between subsequent tabulated args.
      I ida_;    ///< Inverse of difference between subsequent tabulated args.
      std::vector<F> f_; ///< Tabulated return values.

   public:
      /// Initialize table from the first tabulated argument, the difference
      /// between subsequent arguments, and the list of tabulated function
      /// objects.
      dense_table(A const &a,       ///< First tabulated argument.
                  A const &d,       ///< Difference between subsequent args.
                  std::vector<F> vf ///< Tabulated function objects.
                  )
            : a_frst_(a)
            , da_(d)
            , ida_(1.0 / da_)
            , f_(vf)
      {
      }

      /// First tabulated argument.
      A const& a_frst() const { return a_frst_; }

      /// Last tabulated argument.
      A a_last() const { return a_frst_ + (f_.size() - 1) * da_; }

      /// Difference between subsequent tabulated arguments.
      A const& da() const { return da_; }

      /// Tabulated values of function.
      std::vector<F> const& f() const { return f_; }

      /// Value returned by tabulated function object associated with supplied
      /// argument.
      auto operator()(A const& a /**< Desired argument. */) const -> decltype(F()(a))
      {
         if (a < a_frst_) {
            return f_[0](a - a_frst_);
         }
         A const last = a_last();
         if (a > last) {
            return (*f_.rbegin())(a - last);
         }
         int const i((a - a_frst_) * ida_ + 0.5);
         A const ai = a_frst_ + i * da_;
         return f_[i](a - ai);
      }
   };
}

#endif // ndef NUMERIC_DENSE_TABLE_HPP

