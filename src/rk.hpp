
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

#include <ilist.hpp>        // for ilist
#include <sparse-table.hpp> // for sparse_table
#include <util.hpp>         // for RAT

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
      /// Type returned by function to be integrated.
      using DYDX = RAT<Y, X>;

   public:
      /// Type of function to be integrated.
      using func = std::function<DYDX(X)>;

      /// Type of list of values returned by function to be integrated.
      using dlist = ilist<X, DYDX>;

      /// Type of list of partial integrations of function.
      using ylist = ilist<X, Y>;

   private:
      /// Function to be integrated.  This is called \a deriv because
      /// Runge-Kutta integrates a derivative.
      func deriv;

      X      x;     ///< Independent variable.
      Y      y;     ///< Variable accumulated during integration.
      DYDX   dydx;  ///< Value of \a deriv at beginning of interval.
      double tol;   ///< Error tolerance.
      bool   store; ///< True if intermediate values should be stored.
      dlist  dl;    ///< Storage for values from \a deriv.
      ylist  yl;    ///< Storage for integrated values.
      int    nok;   ///< Number of propagations with planned h.
      int    nbad;  ///< Number of propagations with unplanned h.

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
            rkck(h, ytemp, yerr); // Take a trial step.
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

      /// Make sure that tolerance is neither negative nor too small.
      void check_tol()
      {
         double constexpr eps     = std::numeric_limits<double>::epsilon();
         double constexpr min_tol = 100.0 * eps;
         if (tol <= 0.0) {
            throw "tolerance not positive";
         } else if (tol < min_tol) {
            tol = min_tol;
         }
      }

      /// Make sure that n is reasonable, and pick initial guess at stepsize.
      X initial_h(
            /** Lower limit of integration. */ X const &x1,
            /** Upper limit of integration. */ X const &x2,
            /** Number of equal-size steps. */ int &    n)
      {
         if (n < 2) {
            n = 2;
         }
         X const h1 = (X(x2) - X(x1)) / (n - 1);
         if (x2 - x1 > X(0.0)) {
            return +fabs(h1);
         } else {
            return -fabs(h1);
         }
      }

      /// This function is based on `odeint()` on Page 721 of Numerical Recipes
      /// in C, Second Edition.
      void
      init(/** Lower limit of integration. */ X   x1,
           /** Upper limit of integration. */ X   x2,
           /** Number of equal-size steps. */ int n)
      {
         check_tol();
         X h = initial_h(x1, x2, n);
         using namespace std;
         static const PRD<X, X> XSQR_0 = 0.0;
         while (true) {
            dydx = deriv(x);
            static Y const TINY(1.0E-300);
            // General-purpose scaling used to monitor accuracy.
            Y const yscal = fabs(y) + fabs(dydx * h) + TINY;
            if (store) {
               dl.push_back({x, dydx});
               yl.push_back({x, y}); // y=0 on first time through loop.
            }
            X const xh = x + h;
            if ((xh - x2) * (xh - x1) > XSQR_0) {
               h = x2 - x; // Decrease stepsize to avoid overshoot.
            }
            X hdid, hnext;
            rkqs(h, yscal, hdid, hnext);
            if (hdid == h) {
               ++nok;
            } else {
               ++nbad;
            }
            if ((x - x2) * (x2 - x1) >= XSQR_0) {
               if (store) {
                  dl.push_back({x, deriv(x)});
                  yl.push_back({x, y});
               }
               return; // We are done; exit normally.
            }
            if (fabs(hnext) <= 0.0) {
               cerr << "rk_quad: WARNING: step size too small" << endl;
               return;
            }
            h = hnext;
         }
         cerr << "rk_quad: WARNING: too many steps" << endl;
      }

   public:
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

      /// List values returned by function to be integrated. Each of these has
      /// a corresponding element in the list returned by intermed_int().
      dlist const &intermed_fnc() const { return dl; }

      /// List of partial values of definite integral. Each of these has a
      /// corresponding element in the list returned by intermed_fnc().
      ylist const &intermed_int() const { return yl; }

      /// Construct the quadratic interpolant through the integrated function's
      /// evaluated points (the same points as returned by intermed_fnc()), so
      /// that the curvature of each piece is adjusted for agreement with the
      /// Runge-Kutta estimate (returned by intermed_int()) of the integral
      /// over the same piece.
      ///
      /// The coefficients of a quadratic interpolant between two points can be
      /// determined if the area under the interpolant be known. For each piece
      /// between subsequent points
      /// \f[
      ///    p_i = (x_i,y_i)
      /// \f]
      /// \f[
      ///    p_{i+1} = (x_{i+1},y_{i+1})
      /// \f]
      /// returned by intermed_fnc(), there are two corresponding points
      /// \f$(x_i,Y_i)\f$ and \f$(x_{i+1},Y_{i+1})\f$ returned by
      /// intermed_int(). (The notation used here is different from that used
      /// in the documentation for intermed_int().)  The difference,
      /// \f[
      ///    A_i = Y_{i+1} - Y_i
      /// \f]
      /// is the Runge-Kutta estimate of the area under the function that
      /// passes through \f$p_i\f$ and \f$p_{i+1}\f$.
      ///
      /// The linear interpolant between the points has the form
      /// \f[
      ///    \lambda_i(x) = y_i + m_i [x - x_i],
      /// \f]
      /// where the slope is
      /// \f[
      ///    m_i = \frac{\Delta y_i}{\Delta x_i};
      /// \f]
      /// the area under the line is
      /// \f[
      ///    a_i = \Delta x_i \; \left[ y_i + \frac{\Delta y_i}{2} \right];
      /// \f]
      /// and where
      /// \f[
      ///    \Delta x_i = x_{i+1} - x_i
      /// \f]
      /// \f[
      ///    \Delta y_i = y_{i+1} - y_i.
      /// \f]
      ///
      /// The quadratic interpolant that not only passes through \f$p_i\f$ and
      /// \f$p_{i+1}\f$ but also covers the area \f$A_i\f$ is of the form
      /// \f[
      ///    q_i(x) = y_i + [m_i + c_i[x - x_{i+1}]][x - x_i],
      /// \f]
      /// where
      /// \f[
      ///    c_i = -\frac{6 [A_i - a_i]}{[\Delta x_i]^3}.
      /// \f]
      sparse_table<X> make_fnc_interp() const
      {
         if (dl.size() < 2) {
            throw "Must have at least two control points.";
         }
         // Each subdomain is the x region between subsequent control points.
         X const        a0 = 0.5 * (dl[0].first + dl[1].first);
         unsigned const nd = dl.size() - 1; // number of deltas
         using namespace std;
         vector<pair<X, GiNaC::ex>> vf(nd);
         for (unsigned i = 0; i < nd; ++i) {
            unsigned const j = i + 1;
            // horizontal geometry
            X const &x1 = dl[i].first; // left edge of piece
            X const &x2 = dl[j].first; // right edge of piece
            X const  dx = x2 - x1;     // width of piece
            // vertical geometry
            DYDX const &y1 = dl[i].second; // func val at left
            DYDX const &y2 = dl[j].second; // func val at right
            DYDX const  dy = y2 - y1;      // change in func val
            // areas
            Y const a1 = (y1 + 0.5 * dy) * dx;        // under linear interp
            Y const a2 = yl[j].second - yl[i].second; // of R-K estimate
            Y const da = a2 - a1;                     // difference
            // calculation of coefficients
            auto const dx3 = dx * dx * dx;
            auto const c   = -6.0 * da / dx3; // quadratic coef
            auto const m   = dy / dx;         // slope of linear interp
            // construction of expression for piece
            auto const &x = sparse_table<X>::x; // independent var's symbol
            vf[i].first   = dx;                 // width of piece
            vf[i].second  = y1 + (m + c * (x - x2)) * (x - x1); // expression
         }
         return sparse_table<X>(a0, move(vf));
      }

      /// Construct the cubic interpolant through the integral's partial-sum
      /// points (the same points as returned by intermed_int()), so that the
      /// curvature of each piece is adjusted for agreement with the derivative
      /// values (returned by intermed_fnc()) at the endpoints of the same
      /// piece.
      ///
      /// The coefficients of a cubic interpolant between two points can be
      /// determined if the derivative be known at each piece's endpoint.  For
      /// each piece between subsequent points
      /// \f[
      ///    p_i = (x_i,y_i)
      /// \f]
      /// \f[
      ///    p_{i+1} = (x_{i+1},y_{i+1})
      /// \f]
      /// returned by intermed_int(), there are two corresponding derivative
      /// constraints \f$(x_i,y'_i)\f$ and \f$(x_{i+1},y'_{i+1})\f$ returned by
      /// intermed_fnc(). (The notation used here is different from that used
      /// in the documentation for intermed_fnc().)
      ///
      /// The linear interpolant between the points has the form
      /// \f[
      ///    \lambda_i(x) = y_i + m_i [x - x_i],
      /// \f]
      /// where the slope is
      /// \f[
      ///    m_i = \frac{\Delta y_i}{\Delta x_i},
      /// \f]
      /// and where
      /// \f[
      ///    \Delta x_i = x_{i+1} - x_i
      /// \f]
      /// \f[
      ///    \Delta y_i = y_{i+1} - y_i.
      /// \f]
      ///
      /// The cubic interpolant that not only passes through \f$p_i\f$ and
      /// \f$p_{i+1}\f$ but also obeys the constraints on the derivative
      /// returned by intermed_fnc() has the form
      /// \f[
      ///    \kappa_i(x) = y_i + y'_i [x - x_i] + \frac{3m_i - 2y'_i -
      ///    y'_{i+1}}{\Delta x_i} [x - x_i]^2 + \frac{y'_i + y'_{i+1} -
      ///    2m_i}{[\Delta x_i]^2} [x - x_i]^3.
      /// \f]
      sparse_table<X> make_int_interp() const
      {
         if (dl.size() < 2) {
            throw "Must have at least two control points.";
         }
         // Each subdomain is the x region between subsequent control points.
         X const        a0 = 0.5 * (dl[0].first + dl[1].first);
         unsigned const nd = dl.size() - 1; // number of deltas
         using namespace std;
         vector<pair<X, GiNaC::ex>> vf(nd);
         for (unsigned i = 0; i < nd; ++i) {
            unsigned const j = i + 1;
            // horizontal geometry
            X const &x1 = yl[i].first; // left edge of piece
            X const &x2 = yl[j].first; // right edge of piece
            X const  dx = x2 - x1;     // width of piece
            // vertical geometry
            Y const &y1 = yl[i].second; // func val at left
            Y const &y2 = yl[j].second; // func val at right
            Y const  dy = y2 - y1;      // change in func val
            // derivatives
            DYDX const &yp1 = dl[i].second; // derivative at left edge
            DYDX const &yp2 = dl[j].second; // derivative at right edge
            // calculation of coefficients
            auto const dx2 = dx * dx;
            auto const m   = dy / dx; // slope of linear interp
            auto const c = (3.0 * m - 2.0 * yp1 - yp2) / dx; // quadratic coef
            auto const d = (yp1 + yp2 - 2.0 * m) / dx2;      // cubic coef
            // construction of expression for piece
            auto const &x  = sparse_table<X>::x; // independent var's symbol
            auto const  xt = x - x1;             // change of variable
            vf[i].first    = dx;                 // width of piece
            vf[i].second   = y1 + yp1 * xt + c * pow(xt, 2) + d * pow(xt, 3);
         }
         return sparse_table<X>(a0, move(vf));
      }
   };

   /// Short alias for Runge-Kutta solver for ordinary double-precision values.
   using rk_quadd = rk_quad<double, double>;
}

#endif // ndef NUMERIC_RK_HPP

