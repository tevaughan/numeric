
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file integral.hpp
///
/// \brief  Source for num::integral(), which numerically integrates via
///         adaptive quadrature and the trapezoidal rule.

#ifndef NUMERIC_INTEGRAL_HPP
#define NUMERIC_INTEGRAL_HPP

#include <cmath>      // for fabs()
#include <functional> // for function
#include <iostream>   // for cerr
#include <limits>     // for numeric_limits
#include <utility>    // for foward()

#include <integral-stats.hpp> // for integral_stats
#include <interval.hpp>       // for interval and subinterval_stack
#include <util.hpp>           // for RAT and PRD

namespace num
{
   template <typename X, typename Y>
   class interpolant;

   template <typename X, typename Y>
   struct ilist;

   /// Runge-Kutta solver optimized for quadrature.  rk_quad is intended to be
   /// a descendant of a class implementing a more general Runge-Kutta solver.
   ///
   /// \tparam X  Type of the independent variable.
   /// \tparam Y  Type of the variable that is accumulated during integration.
   template <typename X, typename Y>
   class rk_quad
   {
      using DYDX = RAT<Y, X>; ///< Type of the integrand.

      // These constant expressions ought to be moved to the base class.  There
      // are more such constants used by the general solver.  The ones listed
      // here are the ones used by the version optimized for quadrature.

      /// Fraction of step size to Substep 3.
      static double constexpr a3 = 0.3;

      /// Fraction of step size to Substep 4.
      static double constexpr a4 = 0.6;

      /// Fraction of step size to Substep 5.
      static double constexpr a5 = 1.0;

      /// Fraction of step size to Substep 6.
      static double constexpr a6 = 0.875;

      /// Weight for derivative at Substep 1.
      static double constexpr c1 = 37.0 / 378.0;

      /// Weight for derivative at Substep 3.
      static double constexpr c3 = 250.0 / 621.0;

      /// Weight for derivative at Substep 4.
      static double constexpr c4 = 125.0 / 594.0;

      /// Weight for derivative at Substep 6.
      static double constexpr c6 = 512.0 / 1771.0;

      /// Error coefficient at Substep 1.
      static double constexpr dc1 = c1 - 2825.0 / 27648.0;

      /// Error coefficient at Substep 3.
      static double constexpr dc3 = c3 - 18575.0 / 48384.0;

      /// Error coefficient at Substep 4.
      static double constexpr dc4 = c4 - 13525.0 / 55296.0;

      /// Error coefficient at Substep 5.
      static double constexpr dc5 = -277.0 / 14336.0;

      /// Error coefficient at Substep 6.
      static double constexpr dc6 = c6 - 0.25;

      /// Function to be integrated.  This is called \a deriv because
      /// Runge-Kutta integrates a derivative.
      std::function<DYDX(X)> deriv;

      X x;             ///< Independent variable.
      Y y;             ///< Variable accumulated during integration.
      DYDX dydx;       ///< Value of \a deriv at beginning of interval.
      double tol;      ///< Error tolerance.
      bool store_dydx; ///< True if values from \a deriv should be stored.
      bool store_y;    ///< True if values of \a y should be stored.
      ilist<X, DYDX> fi_list;  ///< Storage for values from \a deriv.
      interpolant<X, DYDX> fi; ///< Interpolant for integrand.
      interpolant<X, Y> ii;    ///< Interpolant for integral.
      int nok;                 ///< Number of propagations with planned h.
      int nbad;                ///< Number of propagations with unplanned h.

      /// Given the value for variable \a y and the value for its derivative \a
      /// dydx, use the fifth-order Cash-Karp Runge-Kutta method to advance the
      /// solution for \a y over an interval \a h, and return the incremented
      /// variable as \a yout.  Also, return, by way of the embedded
      /// fourth-order method, an estimate \a yerr of the local truncation
      /// error in \a yout.  The supplied function \a deriv returns the
      /// derivative (just like dydx) at any value of the independent variable.
      ///
      /// This function is based on a simplified version of `rkck()` found on
      /// Page 719 in Numerical Recipes in C, Second Edition.  The
      /// simplification is due mainly to the fact that the present function is
      /// used for quadrature, and so \a deriv does not require \a y as input.
      ///
      void rkck(X const &h, ///< Length of interval.
                Y &yout,    ///< Accumulated value at end of interval.
                Y &yerr     ///< Estimate of local truncation error.
                )
      {
         // 1st step is given as dydx on input.
         // 2nd step not needed because deriv() does not need y as input.
         X const x3 = x + a3 * h;
         X const x4 = x + a4 * h;
         X const x5 = x + a5 * h;
         X const x6 = x + a6 * h;
         RAT<Y, X> const ak3 = deriv(x3); // 3rd step.
         RAT<Y, X> const ak4 = deriv(x4); // 4th step.
         RAT<Y, X> const ak5 = deriv(x5); // 5th step.
         RAT<Y, X> const ak6 = deriv(x6); // 6th step.
         if (store_dydx) {
            fi_list.push_back({x3, ak3});
            fi_list.push_back({x4, ak4});
            fi_list.push_back({x5, ak5});
            fi_list.push_back({x6, ak6});
         }
         // Accumulate increments with proper weights.
         yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
         // Estimate error as difference between fourth- and fifth-order
         // methods.
         yerr = h *
                (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
      }

      /// Fifth-order Runge-Kutta step with monitoring of local truncation
      /// error to ensure accurcy and adjust stepsize. Input are the
      /// accumulated variable \a y and its derivative \a dydx at the starting
      /// value of the independent variable \a x.  Also input are the stepsize
      /// \a htry to be attempted and the required accuracy \a tol.  On output,
      /// \a y and \a x are replaced by their new values, \a hdid is the
      /// stepsize that was actually accomplished, and \a hnext is the
      /// estimated next stepsize.  The supplied function \a deriv returns the
      /// derivative (just like dydx) at any value of the independent variable.
      ///
      /// This function is based on `rkqs()` found on Page 719 in Numerical
      /// Recipes in C, Second Edition.
      ///
      void rkqs(X const &htry,  ///< Stepsize to be attempted.
                Y const &yscal, ///< Scaling used to monitor accuracy.
                X &hdid,        ///< Stepsize that was accomplished.
                X &hnext        ///< Estimated next stepsize.
                )
      {
         static double constexpr SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25;
         static double const ERRCON = pow(5.0 / SAFETY, 1.0 / PGROW);
         double err;
         Y yerr;
         Y ytemp;
         X h = htry; // Set stepsize to the initial trial value.
         while (true) {
            rkck(h, ytemp, yerr); // Take a step.
            err = fabs(yerr / yscal / tol);
            if (err <= 1.0) {
               break;
            }
            X const htemp = SAFETY * h * pow(err, PSHRNK);
            // Truncation error is too large.  Reduce stepsize no more than a
            // factor of 10.
            X const tenth = 0.1 * h;
            static X const zero(0.0);
            if (h >= zero) {
               if (htemp > tenth) {
                  h = htemp;
               } else {
                  h = tenth;
               }
            } else {
               if (htemp < tenth) {
                  h = htemp;
               } else {
                  h = tenth;
               }
            }
            X const xnew = x + h;
            if (xnew == x) {
               std::cerr << "rkqs: WARNING: stepsize underflow" << std::endl;
               h = std::numeric_limits<double>::epsilon() * 100.0 * x;
               break;
            }
         }
         // Increase stepsize no more than a factor of 5.
         if (err > ERRCON) {
            hnext = SAFETY * h * pow(err, PGROW);
         } else {
            hnext = 5.0 * h;
         }
         x += (hdid = h);
         y = ytemp;
      }

      /// This function is based on `odeint()` on Page 721 of Numerical Recipes
      /// in C, Second Edition.
      ///
      void init(X x1, ///< Lower limit of integration.
                X x2, ///< Upper limit of integration.
                int n ///< Inverse of initial step size.
                )
      {
         double constexpr eps = std::numeric_limits<double>::epsilon();
         double constexpr min_tol = 100.0 * eps;
         if (tol <= 0.0) {
            throw "tolerance not positive";
         } else if (tol < min_tol) {
            tol = min_tol;
         }
         X const hmin(0.0);
         if (n < 2) {
            n = 2;
         }
         X const h1 = (X(x2) - X(x1)) / (n - 1);
         X h;
         if (x2 - x1 > X(0.0)) {
            h = +fabs(h1);
         } else {
            h = -fabs(h1);
         }
         ilist<X, Y> ii_list;
         while (true) {
            dydx = deriv(x);
            static Y const TINY(1.0E-300);
            // General-purpose scaling used to monitor accuracy.
            Y const yscal = fabs(y) + fabs(dydx * h) + TINY;
            if (store_dydx) {
               fi_list.push_back({x, dydx});
            }
            if (store_y) {
               ii_list.push_back({x, y});
            }
            X const xh = x + h;
            using XSQR = PRD<X, X>;
            if ((xh - x2) * (xh - x1) > XSQR(0.0)) {
               h = x2 - x; // Decrease stepsize to avoid overshoot.
            }
            X hdid, hnext;
            rkqs(h, yscal, hdid, hnext);
            if (hdid == h) {
               ++nok;
            } else {
               ++nbad;
            }
            if ((x - x2) * (x2 - x1) >= XSQR(0.0)) {
               if (store_dydx) {
                  fi = interpolant<X, DYDX>(fi_list);
               }
               if (store_y) {
                  ii_list.push_back({x, y});
                  ii = interpolant<X, Y>(ii_list);
               }
               return; // We are done; exit normally.
            }
            if (fabs(hnext) <= hmin) {
               std::cerr << "rk_quad: WARNING: step size too small"
                         << std::endl;
               return;
            }
            h = hnext;
         }
         std::cerr << "rk_quad: WARNING: too many steps" << std::endl;
      }

   public:
      /// Numerically integrate a function, and store the result in rk_quad::y.
      /// Use fifth-order Runge-Kutta with adaptive stepsize for quadrature.
      /// The initial guess for the step size is used at the lower limit of
      /// integration.
      ///
      /// Optionally store an approximant to the function integrated and an
      /// approximant to the indefinite integral, zeroed at the lower limit of
      /// integration.
      ///
      /// The optional parameter \a n indicates that the initial step should be
      /// 1/n of the interval of integration.
      ///
      /// \tparam X1  Type of lower limit of integration; X1 must convert to X.
      /// \tparam X2  Type of upper limit of integration; X2 must convert to X.
      template <typename X1, typename X2>
      rk_quad(std::function<DYDX(X)> f, ///< Function to be integrated.
              X1 x1,                    ///< Lower limit of integration.
              X2 x2,                    ///< Upper limit of integration.
              double t = 1.0E-06,       ///< Error tolerance.
              int n = 16,               ///< Inverse of initial step size.
              /// If true, store approximant to function \a f.
              bool s_dydx = false,
              /// If true, store approximant to indefinite integral of \a f.
              bool s_y = false)
         : deriv(f)
         , x(x1)
         , y(0.0)
         , tol(t)
         , store_dydx(s_dydx)
         , store_y(s_y)
         , nok(0)
         , nbad(0)
      {
         init(x1, x2, n);
      }

      /// Numerically integrate a function, and store the result in rk_quad::y.
      /// Use fifth-order Runge-Kutta with adaptive stepsize for quadrature.
      /// The initial guess for the step size is used at the lower limit of
      /// integration.
      ///
      /// Optionally store an approximant to the function integrated and an
      /// approximant to the indefinite integral, zeroed at the lower limit of
      /// integration.
      ///
      /// The optional parameter \a n indicates that the initial step should be
      /// 1/n of the interval of integration.
      ///
      /// \tparam X1  Type of lower limit of integration; X1 must convert to X.
      /// \tparam X2  Type of upper limit of integration; X2 must convert to X.
      template <typename X1, typename X2>
      rk_quad(DYDX (*f)(X),       ///< Function to be integrated.
              X1 x1,              ///< Lower limit of integration.
              X2 x2,              ///< Upper limit of integration.
              double t = 1.0E-06, ///< Error tolerance.
              int n = 16,         ///< Inverse size of initial step.
              /// If true, store approximant to function \a f.
              bool s_dydx = false,
              /// If true, store approximant to indefinite integral.
              bool s_y = false)
         : deriv(f)
         , x(x1)
         , y(0.0)
         , tol(t)
         , store_dydx(s_dydx)
         , store_y(s_y)
         , nok(0)
         , nbad(0)
      {
         init(x1, x2, n);
      }

      /// Value of definite integral.
      Y const &def_int() const { return y; }

      /// Tolerance used for computing definite integral.
      double tolerance() const { return tol; }

      /// Interpolant representing the function that was integrated.
      interpolant<X, DYDX> const &interp_func() const { return fi; }

      /// Interpolant representing indefinite integral.
      interpolant<X, Y> const &interp_indef_int() const { return ii; }
   };

   /// Short alias for Runge-Kutta solver for ordinary double-precision values.
   using rk_quadd = rk_quad<double, double>;

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

