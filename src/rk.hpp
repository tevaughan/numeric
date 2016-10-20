
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   rk.hpp
/// \brief  Definition of num::rk_base and num::rk_quad.

#ifndef NUMERIC_RK_HPP
#define NUMERIC_RK_HPP

#include <functional> // for function
#include <limits>     // for numeric_limits

#include <ilist.hpp>       // for ilist
#include <interpolant.hpp> // for interpolant
#include <util.hpp>        // for RAT

namespace num
{
   /// Base class for rk_quad and, eventually, for the class implementing the
   /// general Runge-Kutta integrator.
   class rk_base
   {
   protected:
      // These constant expressions are used by rk_quad and will also be used
      // by the general Runge-Kutta integrator, when it is implemented.  There
      // are more such constants used only by the general integrator.
      //
      // It might be interesting to see how much slow-down occurs if we
      // implement all calculations via GiNaC's (CLN's) numeric type instead of
      // double.

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
   };

   /// Runge-Kutta integrator optimized for quadrature.
   /// \tparam X  Type of the independent variable.
   /// \tparam Y  Type of the variable that is accumulated during integration.
   template <typename X, typename Y>
   class rk_quad : rk_base
   {
      using DYDX = RAT<Y, X>; ///< Type of the integrand.

      /// Function to be integrated.  This is called \a deriv because
      /// Runge-Kutta integrates a derivative.
      std::function<DYDX(X)> deriv;

      X      x;     ///< Independent variable.
      Y      y;     ///< Variable accumulated during integration.
      DYDX   dydx;  ///< Value of \a deriv at beginning of interval.
      double tol;   ///< Error tolerance.
      bool   store; ///< True if intermediate values should be stored.
      ilist<X, DYDX>       fi_list; ///< Storage for values from \a deriv.
      interpolant<X, DYDX> fi;      ///< Interpolant for integrand.
      interpolant<X, Y>    ii;      ///< Interpolant for integral.
      int nok;                      ///< Number of propagations with planned h.
      int nbad; ///< Number of propagations with unplanned h.

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
      void
      rkck(/** Length of interval.                   */ X const &h,
           /** Accumulated value at end of interval. */ Y &      out,
           /** Estimate of local truncation error.   */ Y &      e)
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
         // Accumulate increments with proper weights.
         out = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
         // Estimate eor as difference between fourth- and fifth-order
         // methods.
         e = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
      }

      /// See implementation of `rkqs()` on Paqe 719 in Numerical Recipes in C,
      /// Second Edition.
      static double constexpr SAFETY = 0.9;

      /// Reduce stepsize no more than a factor of 10.
      /// \return  New stepsize.
      static X reduce_step_size(
            /** Old stepsize. */ X h, /** Truncation error. */ double err)
      {
         static double constexpr PSHRNK = -0.25;
         X const        htemp           = SAFETY * h * pow(err, PSHRNK);
         X const        tenth           = 0.1 * h;
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
         return h;
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
      void
      rkqs(/** Stepsize to be attempted.         */ X const &htry,
           /** Scaling used to monitor accuracy. */ Y const &yscal,
           /** Stepsize that was accomplished.   */ X &      hdid,
           /** Estimated next stepsize.          */ X &      hnext)
      {
         double err;
         Y      yerr;
         Y      ytemp;
         X      h = htry; // initial trial value
         while (true) {
            rkck(h, ytemp, yerr); // Take a step.
            err = fabs(yerr / yscal / tol);
            if (err <= 1.0) {
               break;
            }
            // Truncation error is too large.
            h = reduce_step_size(h, err);
            if (x + h == x) {
               throw "stepsize underflow";
            }
         }
         // Increase stepsize no more than a factor of 5.
         static double constexpr PGROW = -0.2;
         static double const ERRCON    = pow(5.0 / SAFETY, 1.0 / PGROW);
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
      void
      init(/** Lower limit of integration.   */ X   x1,
           /** Upper limit of integration.   */ X   x2,
           /** Inverse of initial step size. */ int n)
      {
         double constexpr eps     = std::numeric_limits<double>::epsilon();
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
         X       h;
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
            if (store) {
               fi_list.push_back({x, dydx});
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
               if (store) {
                  fi_list.push_back({x, deriv(x)});
                  ii_list.push_back({x, y});
                  fi = interpolant<X, DYDX>(fi_list);
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
      /// Type of function to be integrated.
      using func = std::function<DYDX(X)>;

      /// Type of ordinary C function to be integrated.
      typedef DYDX (*cfunc)(X);

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
      rk_quad(
            /** Function to be integrated.            */ func   f,
            /** Lower limit of integration.           */ X1     x1,
            /** Upper limit of integration.           */ X2     x2,
            /** Error tolerance.                      */ double t = 1.0E-06,
            /** Inverse of initial step size.         */ int    n = 16,
            /** Whether to store intermediate values. */ bool   s = false)
         : deriv(f), x(x1), y(0.0), tol(t), store(s), nok(0), nbad(0)
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
      rk_quad(
            /** Function to be integrated.            */ cfunc  f,
            /** Lower limit of integration.           */ X1     x1,
            /** Upper limit of integration.           */ X2     x2,
            /** Error tolerance.                      */ double t = 1.0E-06,
            /** Inverse of initial step size.         */ int    n = 16,
            /** Whether to store intermediate values. */ bool   s = false)
         : deriv(f), x(x1), y(0.0), tol(t), store(s), nok(0), nbad(0)
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
}

#endif // ndef NUMERIC_RK_HPP
