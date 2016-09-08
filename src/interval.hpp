
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

#include <stack> // for stack

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
   /// \tparam A  Type of argument to function.
   /// \tparam R  Type returned by function.
   template <typename A, typename R>
   struct subinterval_stack : public std::stack<interval<A, R>> {
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
            this->push(interval<A, R>{ta, tb, f(ta), f(tb)});
            ta = tb;
         }
      }
   };
}

#endif // ndef NUMERIC_INTERVAL_HPP

