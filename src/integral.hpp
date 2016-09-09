
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file integral.hpp
///
/// \brief  Source for integral(), which numerically integrates via adaptive
///         quadrature and the trapezoidal rule.

#ifndef NUMERIC_INTEGRAL_HPP
#define NUMERIC_INTEGRAL_HPP

#include <algorithm>  // for sort()
#include <cmath>      // for fabs()
#include <functional> // for function
#include <iostream>   // for cerr
#include <limits>     // for numeric_limits
#include <utility>    // for foward()
#include <vector>     // for vector

#include "interval.hpp" // for interval and subinterval_stack

/// Namespace for C++-11 library enabling numerical computation.
namespace num
{
   /// Numerically integrate a function, and return the result.
   ///
   /// Use the trapezoidal rule with adaptive quadrature.
   ///
   /// NOTE: If the function itself call integral(), as for a double, a triple,
   /// or a higher-level integral, then the tolerance specified for the
   /// innermost integral must be the tightest, and the ratio of the tolerance
   /// between any integral and the integral at the next level out must be at
   /// most about 0.05.  So, for example, a double integral with tolerance
   /// 1.0E-03 should have the outer tolerance set to 1.0E-03 and the inner
   /// tolerance set at most to 5.0E-05.  Otherwise, the outer integral will
   /// fruitlessly make too many subdivisions and consume too much time,
   /// possibly never to converge.
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
   /// \param  f   Function to be integrated.
   /// \param  aa  Lower limit of integration.
   /// \param  bb  Upper limit of integration.
   /// \param  t   Error tolerance.
   /// \param  n   Initial number of evenly spaced samples of the function.
   /// \param  ms  Maximum number of samples of function.
   /// \return     Numeric integral of function.
   template <typename R, typename A, typename A1, typename A2>
   auto integral(std::function<R(A)> f, A1 aa, A2 bb, double t = 1.0E-06,
                 unsigned n = 16)
         // See Page 28 of Effective Modern C++ by Scott Meyers.
         -> decltype(std::forward<std::function<R(A)>>(f)(A()) * A())
   {
      if (t <= 0.0) {
         throw "tolerance not positive";
      }
      using I = decltype(std::forward<std::function<R(A)>>(f)(A()) * A());
      double sign = 1.0;
      A a = aa;
      A b = bb;
      if (a > b) {
         std::swap(a, b);
         sign = -1.0;
      }
      subinterval_stack<A, R> s(n, a, b, f);
      std::vector<I> areas;          // area of each trapezoid
      std::vector<I> areas_roundoff; // estimate roundoff errors
      double constexpr eps = std::numeric_limits<double>::epsilon();
      double const c = 2.0 * eps;
      while (s.size()) {
         using interval = interval<A, R>;
         interval const r = s.top();
         s.pop();
         A const len = r.b - r.a;             // length of interval
         A const midp = 0.5 * (r.a + r.b);    // midpoint of interval
         R const mean = 0.5 * (r.fa + r.fb);  // mean of function values
         R const fmid = f(midp);              // function value at midpoint
         R const rmean = 0.5 * (mean + fmid); // refined mean
         // Trapezoidal rule len*mean is estimate for current subinterval.
         // Composite trapezoidal rule 0.5*len*(mean + fmid) is refined
         // estimate obtained by subdividing the current subinterval.
         // Magnitude of difference between them is u1.  Stop refining when u1
         // become smaller than tolerance t times magnitude of refined
         // estimate.
         R const u1 = fabs(mean - rmean);
         R const u2 = fabs(mean) + fabs(rmean);
         R const u3 = fabs(rmean) * t;
         I const ds = rmean * len;
         if (midp * c > len || u2 * c > u3) {
            // Stop refining estimate because of roundoff error.
            areas_roundoff.push_back((mean - rmean) * len);
            areas.push_back(ds);
         } else if (u1 <= u3) {
            // Stop refining estimate because desired accuracy has been
            // reached.
            areas_roundoff.push_back((mean - rmean) * len);
            areas.push_back(ds);
         } else {
            // Continue refining estimate.
            s.push(interval{r.a, midp, r.fa, fmid});
            s.push(interval{midp, r.b, fmid, r.fb});
         }
      }
      auto comp = [](I a1, I a2) { return fabs(a1) < fabs(a2); };
      std::sort(areas.begin(), areas.end(), comp);
      std::sort(areas_roundoff.begin(), areas_roundoff.end(), comp);
      I sum(0.0);
      I sum_roundoff(0.0);
      for (auto a : areas) {
         sum += a;
      }
      for (auto a : areas_roundoff) {
         sum_roundoff += a;
      }
      I const fsr = fabs(sum_roundoff);
      I const fs = fabs(sum);
      if (fsr > fs * t) {
         std::cerr << "integral: WARNING: Fractional round-off error ";
         if (fs > I(0.0)) {
            std::cerr << fsr / fs;
         }
         std::cerr << " > tolerance " << t << " for sum " << sum << std::endl;
      }
      return sign * sum;
   }

   /// Integrate an ordinary function by way of its function pointer.
   /// Transform the passed-in function pointer into an instance of function<>,
   /// and call integrate(function<>,...).
   ///
   /// \tparam R   Type of instance returned by function.
   /// \tparam A   Type of argument to function.
   /// \tparam A1  Type of lower limit of integration; A1 must convert to A.
   /// \tparam A2  Type of upper limit of integration; A2 must convert to A.
   /// \param  f   Function to be integrated.
   /// \param  a   Lower limit of integration.
   /// \param  b   Upper limit of integration.
   /// \param  t   Error tolerance.
   /// \param  n   Initial number of evenly spaced samples of the function.
   /// \return     Numeric integral of function.
   template <typename R, typename A, typename A1, typename A2>
   auto integral(R (*f)(A), A1 a, A2 b, double t = 1.0E-06, unsigned n = 16)
         // See Page 28 of Effective Modern C++ by Scott Meyers.
         -> decltype(std::forward<R (*)(A)>(f)(A()) * A())
   {
      return integral(std::function<R(A)>(f), a, b, t, n);
   }
}

#endif // ndef NUMERIC_INTEGRAL_HPP

