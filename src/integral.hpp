
/// \file integral.hpp
///
/// \brief  Source for integral(), which numerically integrates via adaptive
///         quadrature and the trapezoidal rule.

#ifndef NUMERIC_INTEGRAL_HPP
#define NUMERIC_INTEGRAL_HPP

#include <cmath>      // for fabs()
#include <functional> // for function
#include <iostream>   // for cerr
#include <limits>     // for numeric_limits
#include <stack>      // for stack
#include <utility>    // for foward()

/// Namespace for C++-11 library enabling numerical computation.
namespace num
{
   /// Type of each element on the stack used by integral().
   /// \tparam A  Type of argument to function integrated by integral().
   /// \tparam R  Type returned by function.
   template <typename A, typename R>
   struct integration_interval {
      A a;  ///< Beginning of interval.
      A b;  ///< End of interval.
      R fa; ///< Value of integrand at beginning of interval.
      R fb; ///< Value of integrand at end of interval.
   };

   /// Stack used by integral().
   ///
   /// On construction, an instance of integration_subinterval_stack contains a
   /// partition of the interval of integration into subintervals of equal
   /// length.  As integral() executes, however, each of these is adaptively
   /// subdivided as necessary to achieve the desired accuracy.
   ///
   /// \tparam A  Type of argument to function integrated by integral().
   /// \tparam R  Type returned by function.
   template <typename A, typename R>
   struct integration_subinterval_stack
         : public std::stack<integration_interval<A, R>> {
      /// Initialize stack on construction.
      ///
      /// Evaluate integrand at each of n equally spaced points across the
      /// interval of the variable of integration.  The integrand is evaluated
      /// at a minimum of two points (the end points).  The n-1 corresponding
      /// instances of integration_interval are pushed onto the stack.
      ///
      /// \tparam F  Type of function object representing integrand.
      /// \param  n  Initial number of points at which function is evaluated.
      /// \param  a  Starting value of variable of integration.
      /// \param  b  Ending value of variable of integration.
      /// \param  f  Function that represents integrand.
      template <typename F>
      integration_subinterval_stack(unsigned n, A a, A b, F &&f)
      {
         if (n < 2) {
            n = 2;
         }
         A const d = (b - a) / (n - 1);
         A ta = a;
         // Loop starts at 1 (not zero) to ensure n-1 iterations (not n).
         for (unsigned j = 1; j < n; ++j) {
            A const tb = ta + d;
            this->push(integration_interval<A, R>{ta, tb, f(ta), f(tb)});
            ta = tb;
         }
      }
   };

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
   /// \param  a   Lower limit of integration.
   /// \param  b   Upper limit of integration.
   /// \param  t   Error tolerance.
   /// \param  n   Initial number of evenly spaced samples of the function.
   /// \return     Numeric integral of function.
   template <typename R, typename A, typename A1, typename A2>
   auto integral(std::function<R(A)> f, A1 a, A2 b, double t = 1.0E-06,
                 unsigned n = 16)
         // See Page 28 of Effective Modern C++ by Scott Meyers.
         -> decltype(std::forward<std::function<R(A)>>(f)(a) * a)
   {
      using I = decltype(std::forward<std::function<R(A)>>(f)(a) * a);
      integration_subinterval_stack<A, R> s(n, a, b, f);
      // Allow return value to be either a double or a dimval, whose default
      // constructor produces a zero value.
      I sum(0.0); // return value
      bool tol_achieved = true;
      double tol_worst = t;
      I ds_worst(0.0);
      double constexpr eps = std::numeric_limits<double>::epsilon();
      double constexpr min = std::numeric_limits<double>::min();
      while (s.size()) {
         using interval = integration_interval<A, R>;
         interval const r = s.top();
         s.pop();
         A const len = r.b - r.a;             // length of interval
         A const midp = 0.5 * (r.a + r.b);    // midpoint of interval
         R const mean = 0.5 * (r.fa + r.fb);  // mean of function values
         R const fmid = f(midp);              // function value at midpoint
         R const rmean = 0.5 * (mean + fmid); // refined mean
         // The trapezoidal rule gives len*mean.  The composite trapezoidal
         // rule gives 0.5*len*(mean + fmid).  The first is the current
         // estimate for the current subinterval.  The second is the refined
         // estimate obtained by subdividing the current subinterval.
         //
         // The magnitude of the difference between the current estimate and
         // the refined estimate is u1.  We are done refining our estimate when
         // u1 become smaller than a certain fraction of the magnitude of the
         // refined estimate.  That fraction is just the specified tolerance t.
         R const u1 = fabs(mean - rmean);
         R const u2 = fabs(rmean) * t;
         if (u1 < u2) {
            // Stop refining estimate because desired accuracy has been
            // reached.
            sum += rmean * len;
         } else {
            // Check to see if we should stop refining estimate because of
            // numerical inability to make improvement.
            R const u3 = fabs(u1 - u2);
            double constexpr ulp = 5.0;
            // If the two sides of the inequality u1 < u2 be different by fewer
            // than about ulp last-place units (in terms of machine epsilon),
            // or if the difference between the sides be subnormal, then stop
            // subdividing.  See example at
            // <http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>.
            // integral() must be friend of dimval for R(min) to work when R be
            // a dimval.
            if (u3 <= R(min) || u3 < eps * (u1 + u2) * ulp) {
               // Stop refining estimate because of inability to reach desired
               // accuracy.
               I const ds = rmean * len;
               sum += ds;
               tol_achieved = false;
               R const a2 = fabs(rmean);
               I const ads = fabs(ds);
               if (ads > ds_worst) {
                  ds_worst = ads;
               }
               // integral() must be friend of dimval for R(min) to work when R
               // be a dimval.
               if (a2 > R(min)) {
                  double const tol_cur = u1 / a2;
                  if (tol_cur > tol_worst) {
                     tol_worst = tol_cur;
                  }
               }
            } else {
               // Continue refining estimate.
               s.push(interval{r.a, midp, r.fa, fmid});
               s.push(interval{midp, r.b, fmid, r.fb});
            }
         }
      }
      I const as = fabs(sum);
      if (!tol_achieved && as > I(min) && ds_worst / as > t) {
         std::cerr << "integral: WARNING: worst-case final refinement was "
                   << tol_worst << " > tolerance=" << t << ".\n"
                   << "                   total fractional contribution: "
                   << ds_worst / as << std::endl;
      }
      return sum;
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
         -> decltype(std::forward<R (*)(A)>(f)(a) * a)
   {
      return integral(std::function<R(A)>(f), a, b, t, n);
   }
}

#endif // ndef NUMERIC_INTEGRAL_HPP

