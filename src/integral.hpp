
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   integral.hpp
/// \brief  Definition of num::integral().

#ifndef NUMERIC_INTEGRAL_HPP
#define NUMERIC_INTEGRAL_HPP

#include <cmath>      // for fabs()
#include <functional> // for function
#include <iostream>   // for cerr
#include <limits>     // for numeric_limits

#include <integral-stats.hpp> // for integral_stats
#include <interval.hpp>       // for interval and subinterval_stack
#include <util.hpp>           // for PRD

namespace num
{
   /// Numerically integrate a function, and return the result.
   ///
   /// Use the trapezoidal rule with adaptive quadrature.
   ///
   /// The specified tolerance is on the relative accuracy of the estimated
   /// mean value of the function over each subdivided interval of the variable
   /// of integration.  Every interval is recursively subdivided into equal
   /// halves.  At each stage in the subdivision, the estimated mean height of
   /// the function over the undivided interval is compared to the refined
   /// estimate after subdivision of the interval.  If the magnitude of the
   /// difference between these estimates be smaller than the product of the
   /// tolerance and the refined estimate, then recursion at the current
   /// interval ceases.
   ///
   /// Recursion is not done by way of the function-call stack but rather by
   /// way of a stack with dynamically allocated memory.
   ///
   /// The basic idea for the algorithm is found in the notes of Jim Lambers.
   /// http://www.math.usm.edu/lambers/mat460/fall09/lecture30.pdf
   ///
   /// The present, original C++-11 implementation is due entirely to Thomas E.
   /// Vaughan, who wrote it in 2016.  The meaning of the tolerance is
   /// different in Vaughan's implementation, which also adds a test for the
   /// numerical limitation of double-precision numbers.
   ///
   /// \tparam R   Type of instance returned by function.
   /// \tparam A   Type of argument to function.
   /// \tparam A1  Type of lower limit of integration; A1 must convert to A.
   /// \tparam A2  Type of upper limit of integration; A2 must convert to A.
   /// \return     Numeric integral of function.
   template <typename R, typename A, typename A1, typename A2>
   PRD<R, A>
   integral(std::function<R(A)> f, ///< Function to be integrated.
            A1 aa,                 ///< Lower limit of integration.
            A2 bb,                 ///< Upper limit of integration.
            double t = 1.0E-06,    ///< Error tolerance.
            unsigned n = 16 ///< Initial number of evenly spaced samples.
            )
   {
      double constexpr eps = std::numeric_limits<double>::epsilon();
      double constexpr min_tol = 1000.0 * eps;
      double tol = t;
      if (tol <= 0.0) {
         throw "tolerance not positive";
      } else if (tol < min_tol) {
         tol = min_tol;
      }
      double sign = 1.0;
      A a = aa;
      A b = bb;
      if (a > b) {
         std::swap(a, b);
         sign = -1.0;
      }
      subinterval_stack<A, R> s(n, a, b, f);
      using I = PRD<R, A>;
      integral_stats<I> stats(0.0 * aa * f(aa));
      while (s.size()) {
         using interval = interval<A, R>;
         interval const r = *s.rbegin();
         s.pop_back();
         A const midp = 0.5 * (r.a + r.b);    // midpoint of interval
         R const fmid = f(midp);              // function value at midpoint
         A const len = r.b - r.a;             // length of interval
         R const mean = 0.5 * (r.fa + r.fb);  // mean of function values
         R const rmean = 0.5 * (mean + fmid); // refined mean
         // Trapezoidal rule len*mean is estimate for current subinterval.
         // Composite trapezoidal rule 0.5*len*(mean + fmid) is refined
         // estimate obtained by subdividing the current subinterval.
         // Magnitude of difference between them is u1.  Stop refining when u1
         // become smaller than tolerance t times magnitude of refined
         // estimate.
         R const u0 = fabs(rmean);
         R const u1 = fabs(mean - rmean);
         R const u2 = fabs(mean + rmean);
         R const u3 = u0 * tol;
         I const ds = rmean * len;
         if (u1 <= u3) {
            // Estimated error sufficiently small.  Stop refining estimate.
            stats.add(ds, u1 * len);
         } else if (u1 <= u2 * tol) {
            // Calculated error too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
         } else if (len <= fabs(midp) * tol) {
            // Length of interval too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
         } else if (fabs(ds) <= fabs(stats.area()) * tol) {
            // Increment to integral too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
         } else {
            // Continue refining estimate.
            s.push_back(interval{r.a, midp, r.fa, fmid});
            s.push_back(interval{midp, r.b, fmid, r.fb});
         }
      }
      I const farea = fabs(stats.area());
      I const sigma = stats.stdev(); // statistical error
      I const derr = farea * t;      // desired error
      I const rerr = farea * eps;    // round-off error
      I eerr; // greater of statistical and round-off error
      if (sigma < rerr) {
         eerr = rerr;
      } else {
         eerr = sigma;
      }
      if (eerr > derr) {
         std::cerr << "integral: WARNING: Estimated error " << eerr / farea
                   << " is greater than tolerance " << t << "." << std::endl;
      }
      return sign * stats.area();
   }

   /// Integrate an ordinary function by way of its function pointer.
   /// Transform the passed-in function pointer into an instance of function<>,
   /// and call integrate(function<>,...).
   ///
   /// \tparam R   Type of instance returned by function.
   /// \tparam A   Type of argument to function.
   /// \tparam A1  Type of lower limit of integration; A1 must convert to A.
   /// \tparam A2  Type of upper limit of integration; A2 must convert to A.
   /// \return     Numeric integral of function.
   template <typename R, typename A, typename A1, typename A2>
   PRD<R, A>
   integral(R (*f)(A),          ///< Function to be integrated.
            A1 a,               ///< Lower limit of integration.
            A2 b,               ///< Upper limit of integration.
            double t = 1.0E-06, ///< Error tolerance.
            unsigned n = 16     ///< Initial number of evenly spaced samples.
            )
   {
      return integral(std::function<R(A)>(f), a, b, t, n);
   }
}

#endif // ndef NUMERIC_INTEGRAL_HPP

