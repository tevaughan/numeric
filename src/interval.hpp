
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file interval.hpp
///
/// \brief  Source for class interval and class subinterval_stack, which are
///         used by function integral() and by class interpolant.

#ifndef NUMERIC_INTERVAL_HPP
#define NUMERIC_INTERVAL_HPP

#include <algorithm>  // for sort
#include <functional> // for function
#include <vector>     // for vector

namespace num
{
   /// Type of each element on the stack used by integral() and by class
   /// interpolant. An instance of interval stores the value of a function's
   /// argument and of the function's return value at the end points of an
   /// interval.
   /// \tparam A  Type of argument to function.
   /// \tparam R  Type returned by function.
   template <typename A, typename R>
   struct interval {
      A a;  ///< Function's argument at left end of interval.
      A b;  ///< Function's argument at right end of interval.
      R fa; ///< Function's return value at left end of interval.
      R fb; ///< Function's return value at right end of interval.
   };

   /// Stack used by integral() and by class interpolant.
   ///
   /// On construction, an instance of interval_stack contains a partition of
   /// an interval into subintervals of equal length.  As integral() or one of
   /// interpolant's constructors executes, however, each of these subintervals
   /// is adaptively subdivided as necessary to achieve the desired accuracy.
   ///
   /// The instance of interval_stack is constructed so that the largest
   /// function values are at the top of the stack (the back of the vector).
   /// This allows the adaptive-quadrature integrator, which starts at the top
   /// of the stack, to minimize the time spent in subdividing. If, on
   /// construction, the top of the stack had an interval corresponding to a
   /// region where the function had smaller absolute values, then the
   /// integrator would spend time subidividing for a relative precision that
   /// would be thrown away when larger function absolute values are
   /// encountered farther down the stack.
   ///
   /// \tparam A  Type of argument to function.
   /// \tparam R  Type returned by function.
   template <typename A, typename R>
   struct subinterval_stack : public std::vector<interval<A, R>> {
      using I = interval<A, R>;

      // Allow sorting of intervals by function value.
      static bool compare(I const &iv1, I const &iv2)
      {
         R const a1 = fabs(iv1.fa);
         R const a2 = fabs(iv2.fa);
         R const b1 = fabs(iv1.fb);
         R const b2 = fabs(iv2.fb);
         if (a1 > b1) {
            if (a2 > b2) {
               return a1 < a2;
            } else {
               return a1 < b2;
            }
         } else {
            if (a2 > b2) {
               return b1 < a2;
            } else {
               return b1 < b2;
            }
         }
      };

      /// Initialize stack on construction.
      ///
      /// Evaluate function at each of n equally spaced points across the
      /// interval.  The function is evaluated at a minimum of two points (the
      /// end points).  The n-1 corresponding instances of interval are pushed
      /// onto the stack.
      ///
      /// \tparam F  Type of function.
      /// \param  n  Initial number of points at which function is evaluated.
      /// \param  a  Starting value of function's argument.
      /// \param  b  Ending value of function's argument.
      /// \param  f  Function.
      template <typename F>
      subinterval_stack(unsigned n, A a, A b, F &&f)
      {
         if (n < 2) {
            n = 2;
         }
         A const d = (b - a) / (n - 1);
         A ta = a;
         // Loop starts at 1 (not zero) to ensure n-1 iterations (not n).
         for (unsigned j = 1; j < n; ++j) {
            A const tb = ta + d;
            this->push_back({ta, tb, f(ta), f(tb)});
            ta = tb;
         }
         sort(this->begin(), this->end(), compare);
      }
   };
}

#endif // ndef NUMERIC_INTERVAL_HPP

