
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   interpolant.hpp
///
/// \brief  Definition for each of num::make_const_interp() and
///         num_make_linear_interp().

#ifndef NUMERIC_INTERPOLANT_HPP
#define NUMERIC_INTERPOLANT_HPP

#include <algorithm> // for sort()
#include <iostream>  // for cerr, endl
#include <limits>    // for numeric_limits::epsilon()
#include <memory>    // for unique_ptr
#include <string>    // for string
#include <utility>   // for pair

#include <ilist.hpp>          // for ipoint, ilist
#include <integral-stats.hpp> // for integral_stats
#include <interval.hpp>       // for interval and subinterval_stack
#include <sparse-table.hpp>   // for sparse_table

namespace num
{
   /// Construct a (\ref sparse_table) piecewise-constant interpolant from the
   /// first two space-delimited columns in an ASCII file.  Each line of the
   /// file must consist either
   ///
   /// - of only white space, optionally followed by the '#' character and a
   ///   subsequent comment or
   ///
   /// - of white space followed by at least two space-delimited floating-point
   ///   numbers (after which everything else on the line is ignored).
   ///
   /// \tparam X  Type of first column, representing the x coordinate.
   /// \tparam Y  Type of second column, representing the Y coordinate.
   template <typename X = double, typename Y = double>
   sparse_table<X> make_const_interp(
         /** Name of ASCII file.     */ std::string  file,
         /** Unit multiplying first col. */ X const &xu = 1,
         /** Unit multiplying secnd col. */ Y const &yu = 1)
   {
      using V    = ilist<X, Y>;
      V const cp = get_points(file, xu, yu); // control points
      if (cp.size() < 1) {
         throw "Must have at least one control point.";
      }
      // For linear interpolation, each subdomain is just the x region between
      // subsequent control points.  For constant interpolation, however, each
      // subdomain is the x region between subsequent *midpoints* between
      // control points, each between a pair of control points.  However, in
      // general only each of the first and last points lies in the center of
      // its subdomain.
      V const cpmid = midpoints(cp); // control midpoints
      X const a0    = cp[0].first;
      std::vector<std::pair<X, GiNaC::ex>> vf(cp.size());
      vf[0].first  = 2.0 * (cpmid[0].first - cp[0].first);
      vf[0].second = cp[0].second;
      for (unsigned i = 1; i < cpmid.size(); ++i) {
         vf[i].first  = cpmid[i].first - cpmid[i - 1].first;
         vf[i].second = cp[i].second;
      }
      unsigned const i = cpmid.size();
      vf[i].first      = 2.0 * (cp[i].first - cpmid[i - 1].first);
      vf[i].second     = cp[i].second;
      return sparse_table<X>(a0, move(vf));
   }

   /// Construct a (\ref sparse_table) piecewise-linear interpolant from a set
   /// of ordered pairs.
   ///
   /// \tparam X  Type of first element of each ordered pair.
   /// \tparam Y  Type of second element of each ordered pair.
   template <typename X = double, typename Y = double>
   sparse_table<X> make_linear_interp(ilist<X, Y> cp)
   {
      if (cp.size() == 0) {
         return sparse_table<X>();
      } else if (cp.size() == 1) {
         throw "Must have at least two control points.";
      }
      std::sort(cp.begin(), cp.end());
      // For linear interpolation, each subdomain is just the x region between
      // subsequent control points.
      X const        a0      = 0.5 * (cp[0].first + cp[1].first);
      unsigned const ndeltas = cp.size() - 1;
      using namespace std;
      vector<pair<X, GiNaC::ex>> vf(ndeltas);
      for (unsigned i = 0; i < ndeltas; ++i) {
         unsigned const j  = i + 1;
         X const &      x1 = cp[i].first;
         X const &      x2 = cp[j].first;
         Y const &      y1 = cp[i].second;
         Y const &      y2 = cp[j].second;
         X const        dx = x2 - x1;
         auto const     c1 = (y2 - y1) / dx;
         Y const        c0 = y1 - c1 * x1;
         vf[i].first       = dx;
         vf[i].second      = c0 + c1 * sparse_table<X>::x;
      }
      return sparse_table<X>(a0, move(vf));
   }

   /// Construct a (\ref sparse_table) piecewise-linear interpolant from the
   /// first two space-delimited columns in an ASCII file.  Each line of the
   /// file must consist either
   ///
   /// - of only white space, optionally followed by the '#' character and a
   ///   subsequent comment or
   ///
   /// - of white space followed by at least two space-delimited floating-point
   ///   numbers (after which everything else on the line is ignored).
   ///
   /// \tparam X  Type of first column, representing the x coordinate.
   /// \tparam Y  Type of second column, representing the Y coordinate.
   template <typename X = double, typename Y = double>
   sparse_table<X> make_linear_interp(
         /** Name of ASCII file.     */ std::string  file = "",
         /** Unit multiplying first col. */ X const &xu   = 1,
         /** Unit multiplying secnd col. */ Y const &yu   = 1)
   {
      using namespace std;
      if (file == "") {
         return sparse_table<X>();
      }
      return make_linear_interp(get_points(file, xu, yu));
   }

   /// Initialize control points from copy of initial subinterval_stack.
   /// \tparam X  Type of independent variable.
   /// \tparam Y  Type of dependent variable.
   template <typename X, typename Y>
   void init_from_stack(
         /** Copy of stack.                  */ subinterval_stack<X, Y> s,
         /** Initially empty list of points. */ ilist<X, Y> &           d)
   {
      using interval = interval<X, Y>;
      // Default ordering of stack is not what we want here.
      auto comp = [](interval const &i1, interval const &i2) {
         return i1.a < i2.a;
      };
      // Work from right-most to left-most subinterval.
      std::sort(s.begin(), s.end(), comp);
      while (s.size() > 1) {
         // Push right-most end of every subinterval except for left-most
         // subinterval.
         interval const r = *s.rbegin();
         s.pop_back();
         d.push_back({r.b, r.fb});
      }
      // For left-most subinterval, push both ends.
      interval const r = *s.rbegin();
      s.pop_back();
      d.push_back({r.b, r.fb});
      d.push_back({r.a, r.fa});
   }

   /// Construct a (\ref sparse_table) piecewise-linear interpolant for a
   /// continuous function over the specified interval of its domain. The
   /// initial number of evenly spaced samples should be sufficient to allow
   /// recursive subdivision to produce an interpolant with the desired
   /// fractional error tolerance.  The numerical integral of the function is
   /// computed as a by-product and---if the user supply a pointer in the final
   /// argument---is then stored at the location indicated by the supplied
   /// pointer.
   ///
   /// \tparam X   Type of independent variable.
   /// \tparam Y   Type of dependent variable.
   /// \param  f   Function to approximate via interpolation.
   /// \param  aa  Left edge of domain.
   /// \param  bb  Right edge of domain.
   /// \param  t   Fractional tolerance of approximation.
   /// \param  n   Initial number of evenly spaced samples of function.
   /// \param  i   If non-null, pointer to storage integral.
   template <typename X, typename Y>
   sparse_table<X> make_linear_interp(
         std::function<Y(X)> f, X aa, X bb, double t = 1.0E-06,
         unsigned n = 16, decltype(X() * Y()) *i = nullptr)
   {
      double constexpr eps     = std::numeric_limits<double>::epsilon();
      double constexpr min_tol = 1000.0 * eps;
      double tol               = t;
      if (tol <= 0.0) {
         throw "tolerance not positive";
      } else if (tol < min_tol) {
         tol = min_tol;
      }
      double sign = 1.0;
      if (aa > bb) {
         std::swap(aa, bb);
         sign = -1.0;
      }
      subinterval_stack<X, Y> s(n, aa, bb, f); // Stack of intervals.
      ilist<X, Y>             d;               // Control points.
      init_from_stack(s, d);                   // Add initial n points to d.
      using A = decltype(X() * Y());
      integral_stats<A> stats(0.0 * aa * f(aa));
      while (s.size()) {
         using interval   = interval<X, Y>;
         interval const r = *s.rbegin();
         s.pop_back();
         X const midp  = 0.5 * (r.a + r.b);   // midpoint of interval
         Y const fmid  = f(midp);             // function value at midpoint
         X const len   = r.b - r.a;           // length of interval
         Y const mean  = 0.5 * (r.fa + r.fb); // mean of function values
         Y const rmean = 0.5 * (mean + fmid); // refined mean
         Y const u0    = fabs(rmean);
         Y const u1    = fabs(mean - rmean);
         Y const u2    = fabs(mean + rmean);
         Y const u3    = u0 * tol;
         A const ds    = rmean * len;
         if (u1 <= u3) {
            // Estimated error sufficiently small.  Stop refining estimate.
            stats.add(ds, u1 * len);
            d.push_back({midp, fmid});
         } else if (u1 <= u2 * tol) {
            // Calculated error too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
            d.push_back({midp, fmid});
         } else if (len <= fabs(midp) * tol) {
            // Length of interval too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
            d.push_back({midp, fmid});
         } else if (fabs(ds) <= fabs(stats.area()) * tol) {
            // Increment to integral too small.  Stop refining estimate.
            stats.add(ds, u1 * len);
            d.push_back({midp, fmid});
         } else {
            // Continue subdividing.
            s.push_back(interval{r.a, midp, r.fa, fmid});
            s.push_back(interval{midp, r.b, fmid, r.fb});
         }
      }
      A const farea = fabs(stats.area());
      A const sigma = stats.stdev(); // statistical error
      A const derr  = farea * t;     // desired error
      A const rerr  = farea * eps;   // round-off error
      A       eerr; // greater of statistical and round-off error
      if (sigma < rerr) {
         eerr = rerr;
      } else {
         eerr = sigma;
      }
      if (eerr > derr) {
         std::cerr << "integral: WARNING: Estimated error " << eerr / farea
                   << " is greater than tolerance " << t << "." << std::endl;
      }
      if (i) {
         *i = sign * stats.area();
      }
      return make_linear_interp(d);
   }

   /// Construct a (\ref sparse_table) piecewise-linear interpolant for a
   /// continuous function over the specified interval of its domain. The
   /// initial number of evenly spaced samples should be sufficient to allow
   /// recursive subdivision to produce an interpolant with the desired
   /// fractional error tolerance.  The numerical integral of the function is
   /// computed as a by-product and---if the user supply a pointer in the final
   /// argument---is then stored at the location indicated by the supplied
   /// pointer.
   ///
   /// \tparam X   Type of independent variable.
   /// \tparam Y   Type of dependent variable.
   /// \param  f   Function to approximate via interpolation.
   /// \param  aa  Left edge of domain.
   /// \param  bb  Right edge of domain.
   /// \param  t   Fractional tolerance of approximation.
   /// \param  n   Initial number of evenly spaced samples of function.
   /// \param  i   If non-null, pointer to storage integral.
   template <typename X = double, typename Y = double>
   sparse_table<X> make_linear_interp(
         Y (*f)(X), X aa, X bb, double t = 1.0E-06, unsigned n = 16,
         decltype(X() * Y()) *i = nullptr)
   {
      return make_linear_interp(std::function<Y(X)>(f), aa, bb, t, n, i);
   }
}

#endif // ndef NUMERIC_INTERPOLANT_HPP

