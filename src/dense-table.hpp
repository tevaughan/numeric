
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
   /// dense_table stores a function's value at each point on a regular grid
   /// across the value of the function's argument.  So each pair of subsequent
   /// table entries is separated by the same distance \f$\Delta a\f$ in the
   /// space of the function's argument.
   ///
   /// For log-time-lookup, see sparse_table.
   ///
   /// \tparam A  Type of function's argument.
   /// \tparam R  Type of function's return value.
   template<typename A, typename R>
   class dense_table
   {
      using I = decltype(1.0 / A());
      A a0_;    ///< First tabulated argument.
      A a_min_; ///< Minimum argument for tabulated value.
      A a_max_; ///< Maximum argument for tabulated value.
      A da_;    ///< Difference between subsequent tabulated args.
      I ida_;   ///< Inverse of difference between subsequent tabulated args.
      std::vector<R> r_; ///< Tabulated return values.

   public:
      /// Initialize table from the first tabulated argument, the difference
      /// between subsequent arguments, and the list of tabulated return
      /// values.
      dense_table(A const &a,       ///< First tabulated argument.
                  A const &d,       ///< Difference between subsequent args.
                  std::vector<R> vr ///< Tabulated return values.
                  )
            : a0_(a)
            , a_min_(a - 0.5 * d)
            , a_max_(a + (vr.size() - 0.5) * d)
            , da_(d)
            , ida_(1.0 / da_)
            , r_(vr)
      {
      }

      /// First tabulated argument.
      A const& a0() const { return a0_; }

      /// Minimum argument associated with a tabulated value.
      A const& a_min() const { return a_min_; }

      /// Maximum argument associated with a tabulated value.
      A const& a_max() const { return a_max_; }

      /// Difference between subsequent tabulated arguments.
      A const& da() const { return da_; }

      /// Tabulated values of function.
      std::vector<R> const& r() const { return r_; }

      /// Tabulated value associated with supplied argument; zero if supplied
      /// argument be associated with no tabulated value.
      R const& operator()(A const& a /**< Desired argument. */) const
      {
         if (a < a_min_ || a > a_max_) {
            return R(0);
         }
         return r_[int((a - a0_) * ida_ + 0.5)];
      }
   };
}

#endif // ndef NUMERIC_DENSE_TABLE_HPP

