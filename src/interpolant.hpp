
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   interpolant.hpp
///
/// \brief  Source code for num::interpolant, num::ipoint, num::ilist,
///         num::interpolantd, and num::ilistd.

#ifndef NUMERIC_INTERPOLANT_HPP
#define NUMERIC_INTERPOLANT_HPP

#include <algorithm> // for sort()
#include <fstream>   // for ifstream
#include <iostream>  // for cerr, endl
#include <limits>    // for numeric_limits::epsilon()
#include <memory>    // for unique_ptr
#include <sstream>   // for istringstream
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

#include <integral-stats.hpp> // for integral_stats
#include <interval.hpp>       // for interval and subinterval_stack
#include <sparse-table.hpp>   // for sparse_table
#include <cpoly.hpp>          // for cpoly

namespace num
{
   /// Namespace for utility functions used by numeric.
   namespace util
   {
      /// Extract point from space-delimited ASCII line.
      /// \tparam X  Type of first column, representing X coordinate.
      /// \tparam Y  Type of second column, representing Y coordinate.
      template <typename X, typename Y>
      std::pair<X, Y>
      get_point(/** Line of non-blank ASCII text. */ std::string line,
                /** Unit to multiply against first column. */ X xu,
                /** Unit to multiply against secnd column. */ Y yu)
      {
         std::istringstream iss(line);
         double x, y;
         if (!(iss >> x)) {
            throw "error reading x";
         }
         if (!(iss >> y)) {
            throw "error reading y";
         }
         return std::pair<X, Y>(x * xu, y * yu);
      }

      /// Extract a set of points from a space-delimited ASCII file.  Each line
      /// of the file must consist either
      ///
      /// - of only white space, optionally followed by the '#' character and a
      ///   subsequent comment or
      ///
      /// - of white space followed by at least two space-delimited
      ///   floating-point numbers (after which everything else on the line is
      ///   ignored).
      ///
      /// \tparam X  Type of first column, representing X coordinate.
      /// \tparam Y  Type of second column, representing Y coordinate.
      template <typename X, typename Y>
      std::vector<std::pair<X, Y>>
      get_points(/** Name of ASCII file. */ std::string file,
                 /** Unit to multiply against first column. */ X xu,
                 /** Unit to multiply against secnd column. */ Y yu)
      {
         using namespace std;
         ifstream is(file);
         if (!is) {
            throw "Failed to open '" + file + "'.";
         }
         string line;
         vector<pair<X, Y>> points;
         while (getline(is, line)) {
            size_t const p = line.find_first_not_of(" \f\n\r\t");
            if (p == string::npos || line[p] == '#') {
               continue;
            }
            points.push_back(get_point(line, xu, yu));
         }
         auto comp = [](pair<X, Y> const &a, pair<X, Y> const &b) {
            return a.first < b.first;
         };
         sort(points.begin(), points.end(), comp);
         return points;
      }

      /// Return the list of midpoints, each between a subsequent pair of
      /// points in a sorted list.
      ///
      /// \tparam V  Type of vector of points.
      /// \return    Midpoints sorted by x coordinate.
      template <typename V>
      V midpoints(/** Points sorted by x coordinate. */ V const &v)
      {
         if (v.size() < 2) {
            return V();
         }
         V rv(v.size() - 1); // Return value.
         for (unsigned j = 1; j < v.size(); ++j) {
            unsigned const i = j - 1;
            auto const &xi = v[i].first, &yi = v[i].second;
            auto const &xj = v[j].first, &yj = v[j].second;
            rv[i] = {0.5 * (xi + xj), 0.5 * (yi + yj)};
         }
         return rv;
      }
   }

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
   sparse_table<X, cpoly<0, X, Y>>
   make_const_interp(/** Name of ASCII file. */ std::string file,
                     /** Unit multiplying first col. */ X const &xu = 1,
                     /** Unit multiplying secnd col. */ Y const &yu = 1)
   {
      using namespace std;
      using V = vector<pair<X, Y>>;
      using C = cpoly<0, X, Y>;
      V const cp = util::get_points(file, xu, yu); // control points
      if (cp.size() < 1) {
         throw "Must have at least one control point.";
      }
      // For linear interpolation, each subdomain is just the x region between
      // subsequent control points.  For constant interpolation, however, each
      // subdomain is the x region between subsequent *midpoints* between
      // control points, each between a pair of control points.  However, in
      // general only each of the first and last points lies in the center of
      // its subdomain.
      V const cpmid = util::midpoints(cp); // control midpoints
      X const a0 = cp[0].first;
      vector<pair<X, C>> vf(cp.size());
      vf[0].first = 2.0 * (cpmid[0].first - cp[0].first);
      vf[0].second = cp[0].second;
      for (unsigned i = 1; i < cpmid.size(); ++i) {
         vf[i].first = cpmid[i].first - cpmid[i - 1].first;
         vf[i].second = cp[i].second;
      }
      unsigned const i = cpmid.size();
      vf[i].first = 2.0 * (cp[i].first - cpmid[i - 1].first);
      vf[i].second = cp[i].second;
      return sparse_table<X, C>(a0, move(vf));
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
   sparse_table<X, cpoly<1, X, Y>>
   make_linear_interp(/** Name of ASCII file. */ std::string file,
                      /** Unit multiplying first col. */ X const &xu = 1,
                      /** Unit multiplying secnd col. */ Y const &yu = 1)
   {
      using namespace std;
      using V = vector<pair<X, Y>>;
      using C = cpoly<1, X, Y>;
      V const cp = util::get_points(file, xu, yu); // control points
      if (cp.size() < 2) {
         throw "Must have at least two control points.";
      }
      // For linear interpolation, each subdomain is just the x region between
      // subsequent control points.
      X const a0 = 0.5 * (cp[0].first + cp[1].first);
      unsigned const ndeltas = cp.size() - 1;
      vector<pair<X, C>> vf(ndeltas);
      for (unsigned i = 0; i < ndeltas; ++i) {
         unsigned const j = i + 1;
         X const delta = cp[j].first - cp[i].first;
         Y const c0 = 0.5 * (cp[i].second + cp[j].second);
         auto const c1 = (cp[j].second - cp[i].second) / delta;
         vf[i].first = delta;
         vf[i].second.template set_coef<0>(c0);
         vf[i].second.template set_coef<1>(c1);
      }
      return sparse_table<X, C>(a0, move(vf));
   }

   template <char TI, char D, char M, char C, char TE>
   class statdim;

   /// Point used as one of the constraints of a linear interpolant.
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   struct ipoint : public std::pair<I, D> {
      using std::pair<I, D>::pair;
   };

   /// Send point to output stream.
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   std::ostream &operator<<(std::ostream &os, ipoint<I, D> const &p)
   {
      return os << p.first << " " << p.second;
   }

   /// List of points used to constrain a linear interpolant.
   ///
   /// An instance of ilist can be used to initialize an instance of
   /// interpolant.
   ///
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   struct ilist : public std::vector<ipoint<I, D>> {
      using std::vector<ipoint<I, D>>::vector;
   };

   /// Send ilist to output stream.
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   std::ostream &operator<<(std::ostream &os, ilist<I, D> const &list)
   {
      for (auto const &p : list) {
         os << p << "\n";
      }
      return os;
   }

   /// Linear interpolant.
   ///
   /// This original C++ implementation is due entirely to Thomas E. Vaughan,
   /// who wrote it in 2016.
   ///
   /// Class interpolant is designed with the idea that either of its template
   /// type parameters may be a double or else a dimval.
   ///
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   class interpolant
   {
      using point = ipoint<I, D>;    ///< Type of points defining interpolant.
      using list = ilist<I, D>;      ///< Type of list of points.
      using A = decltype(I() * D()); ///< Type of area under interpolant.

      list d_; ///< x-y points between which to interpolate.

      /// Value of definite integral between first and last point defining
      /// interpolant.
      A integral_;

      /// Compare points so that they can be sorted.
      /// \return True only if p1 lie to the left of p2.
      static bool cmpr_pts(point const &p1, point const &p2)
      {
         return p1.first < p2.first;
      }

      /// Sort points.
      void sort() { std::sort(d_.begin(), d_.end(), cmpr_pts); }

      /// Open file.
      /// \param fname  Name of file.
      /// \return       Pointer to input stream from file.
      static std::unique_ptr<std::ifstream> ifstr(std::string fname)
      {
         using namespace std;
         unique_ptr<ifstream> p(new ifstream(fname));
         if (!(*p)) {
            throw "failed to open '" + fname + "'";
         }
         return p;
      }

      /// Extract point from space-delimited ASCII line.
      /// \param line    Line of text from input ASCII file.
      /// \param x_unit  Factor to multiply against first column.
      /// \param y_unit  Factor to multiply against second column.
      static point get_point(std::string line, I x_unit, D y_unit)
      {
         using namespace std;
         istringstream iss(line);
         double x, y;
         if (!(iss >> x)) {
            throw "error reading independent variable";
         }
         if (!(iss >> y)) {
            throw "error reading dependent variable";
         }
         return point(x * x_unit, y * y_unit);
      }

      /// combine() needs to be friendly with another kind of interpolant.
      template <typename OI, typename OD>
      friend class interpolant;

      /// Combine the present instance I1 with another interpolant I2 such that
      /// there is a point in the combined interpolant I3 at every unique x
      /// coordinate across I1 and I2.  Each new y coordinate is the
      /// combination of a y coordinate, in either I1 or I2, and the
      /// corresponding interpolated y coordinate in the other.  The
      /// combination is represented by a passed-in function.
      ///
      /// combine() is used to support multiplication and division of
      /// interpolants.
      ///
      /// \tparam OD    Type of y coordinate in other interpolant.
      /// \param  oint  Reference to other interpolant.
      /// \param  f     Combining function.
      /// \return       Combined interpolant.
      template <typename OD>
      auto combine(interpolant<I, OD> const &oint,
                   std::function<decltype(D() * OD())(D, OD)> f) const
            -> interpolant<I, decltype(D() * OD())>
      {
         using PD = decltype(D() * OD());
         ilist<I, PD> pl;
         ilist<I, D> const &l1 = d_;
         ilist<I, OD> const &l2 = oint.d_;
         unsigned const s1 = l1.size();
         unsigned const s2 = l2.size();
         unsigned n1 = 0;
         unsigned n2 = 0;
         I ind;
         if (s1 == 0 || s2 == 0) {
            return pl;
         }
         while (n1 < s1 || n2 < s2) {
            if (n1 < s1 && n2 < s2) {
               I const i1 = l1[n1].first;
               I const i2 = l2[n2].first;
               D const d1 = l1[n1].second;
               OD const d2 = l2[n2].second;
               if (i1 < i2) {
                  pl.push_back({i1, f(d1, oint(i1))});
                  ++n1;
               } else {
                  pl.push_back({i2, f((*this)(i2), d2)});
                  ++n2;
               }
            } else if (n1 < s1) {
               I const i1 = l1[n1].first;
               D const d1 = l1[n1].second;
               pl.push_back({i1, f(d1, oint(i1))});
               ++n1;
            } else {
               I const i2 = l2[n2].first;
               OD const d2 = l2[n2].second;
               pl.push_back({i2, f((*this)(i2), d2)});
               ++n2;
            }
         }
         return pl;
      }

      /// Initialize control points from copy of initial subinterval_stack.
      void init_from_stack(subinterval_stack<I, D> s)
      {
         using interval = interval<I, D>;
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
            d_.push_back({r.b, r.fb});
         }
         // For left-most subinterval, push both ends.
         interval const r = *s.rbegin();
         s.pop_back();
         d_.push_back({r.b, r.fb});
         d_.push_back({r.a, r.fa});
      }

      /// Initialize from a continuous function and from an interval on its
      /// domain.  The initial number of evenly spaced samples should be
      /// sufficient to allow recursive subdivision to produce an interpolant
      /// with the desired fractional tolerance as an approximation of the
      /// function.
      ///
      /// \tparam A1  Type of value indicating left edge of domain.  An
      ///             instance of A1 must be convertible to type I.
      ///
      /// \tparam A2  Type of value indicating right edge of domain.  An
      ///             instance of A2 must be convertible to type I.
      ///
      /// \param  f   Function to approximate via interpolation.
      /// \param  aa  Left edge of domain.
      /// \param  bb  Right edge of domain.
      /// \param  t   Fractional tolerance of approximation.
      /// \param  n   Initial number of evenly spaced samples of function.
      template <typename A1, typename A2>
      void init(std::function<D(I)> f, A1 aa, A2 bb, double t = 1.0E-06,
                unsigned n = 16)
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
         I a = aa;
         I b = bb;
         if (a > b) {
            std::swap(a, b);
            sign = -1.0;
         }
         subinterval_stack<I, D> s(n, a, b, f);
         init_from_stack(s); // Add initial n points to d_.
         integral_stats<A> stats(0.0 * aa * f(aa));
         while (s.size()) {
            using interval = interval<I, D>;
            interval const r = *s.rbegin();
            s.pop_back();
            I const midp = 0.5 * (r.a + r.b);    // midpoint of interval
            D const fmid = f(midp);              // function value at midpoint
            I const len = r.b - r.a;             // length of interval
            D const mean = 0.5 * (r.fa + r.fb);  // mean of function values
            D const rmean = 0.5 * (mean + fmid); // refined mean
            D const u0 = fabs(rmean);
            D const u1 = fabs(mean - rmean);
            D const u2 = fabs(mean + rmean);
            D const u3 = u0 * tol;
            A const ds = rmean * len;
            if (u1 <= u3) {
               // Estimated error sufficiently small.  Stop refining estimate.
               stats.add(ds, u1 * len);
               d_.push_back({midp, fmid});
            } else if (u1 <= u2 * tol) {
               // Calculated error too small.  Stop refining estimate.
               stats.add(ds, u1 * len);
               d_.push_back({midp, fmid});
            } else if (len <= fabs(midp) * tol) {
               // Length of interval too small.  Stop refining estimate.
               stats.add(ds, u1 * len);
               d_.push_back({midp, fmid});
            } else if (fabs(ds) <= fabs(stats.area()) * tol) {
               // Increment to integral too small.  Stop refining estimate.
               stats.add(ds, u1 * len);
               d_.push_back({midp, fmid});
            } else {
               // Continue subdividing.
               s.push_back(interval{r.a, midp, r.fa, fmid});
               s.push_back(interval{midp, r.b, fmid, r.fb});
            }
         }
         A const farea = fabs(stats.area());
         A const sigma = stats.stdev(); // statistical error
         A const derr = farea * t;      // desired error
         A const rerr = farea * eps;    // round-off error
         A eerr; // greater of statistical and round-off error
         if (sigma < rerr) {
            eerr = rerr;
         } else {
            eerr = sigma;
         }
         if (eerr > derr) {
            std::cerr << "integral: WARNING: Estimated error " << eerr / farea
                      << " is greater than tolerance " << t << "."
                      << std::endl;
         }
         integral_ = sign * stats.area();
         sort();
      }

   public:
      /// Initialize ifrom an ilist.
      /// \param d  Data stored in ilist.
      interpolant(list d = list()) : d_(std::move(d))
      {
         if (d_.size() == 0) {
            integral_ = 0.0 * A();
            return;
         }
         if (d_.size() == 1) {
            integral_ = 0.0 * d[0].first * d[0].second;
         }
         sort();
         integral_ = integral(d_.begin()->first, d_.rbegin()->first);
      }

      /// Initialize from the first two columns of a space-delimited ASCII
      /// file.
      ///
      /// \param fname   Name of input file.
      ///
      /// \param x_unit  Factor multiplied against each number in file's first
      ///                space-delimited column. Product is stored as x value.
      ///
      /// \param y_unit  Factor multiplied against each number in file's second
      ///                space-delimited column. Product is stored as y value.
      interpolant(std::string fname, I x_unit = 1.0, D y_unit = 1.0)
      {
         using namespace std;
         unique_ptr<ifstream> is = ifstr(fname);
         string line;
         while (getline(*is, line)) {
            size_t const p = line.find_first_not_of(" \f\n\r\t");
            if (p == string::npos || line[p] == '#') {
               continue;
            }
            d_.push_back(get_point(line, x_unit, y_unit));
         }
         if (d_.size() < 2) {
            integral_ = 0.0 * x_unit * y_unit;
            return;
         }
         sort();
         integral_ = integral(d_.begin()->first, d_.rbegin()->first);
      }

      /// Initialize from continuous function and interval of its domain. The
      /// initial number of evenly spaced samples should be sufficient to allow
      /// recursive subdivision to produce an interpolant with the desired
      /// fractional tolerance as an approximation of the function.
      ///
      /// \tparam A1  Type of value indicating left edge of domain.  An
      ///             instance of A1 must be convertible to type I.
      ///
      /// \tparam A2  Type of value indicating right edge of domain.  An
      ///             instance of A2 must be convertible to type I.
      ///
      /// \param  f   Function to approximate via interpolation.
      /// \param  aa  Left edge of domain.
      /// \param  bb  Right edge of domain.
      /// \param  t   Fractional tolerance of approximation.
      /// \param  n   Initial number of evenly spaced samples of function.
      template <typename A1, typename A2>
      interpolant(std::function<D(I)> f, A1 aa, A2 bb, double t = 1.0E-06,
                  unsigned n = 16)
      {
         init(f, aa, bb, t, n);
      }

      /// Initialize from continuous function and interval of its domain. The
      /// initial number of evenly spaced samples should be sufficient to allow
      /// recursive subdivision to produce an interpolant with the desired
      /// fractional tolerance as an approximation of the function.
      ///
      /// \tparam A1  Type of value indicating left edge of domain.  An
      ///             instance of A1 must be convertible to type I.
      ///
      /// \tparam A2  Type of value indicating right edge of domain.  An
      ///             instance of A2 must be convertible to type I.
      ///
      /// \param  f   Function to approximate via interpolation.
      /// \param  aa  Left edge of domain.
      /// \param  bb  Right edge of domain.
      /// \param  t   Fractional tolerance of approximation.
      /// \param  n   Initial number of evenly spaced samples of function.
      template <typename A1, typename A2>
      interpolant(D (*f)(I), A1 aa, A2 bb, double t = 1.0E-06, unsigned n = 16)
      {
         init(std::function<D(I)>(f), aa, bb, t, n);
      }

      /// \return  Control points for interpolant.
      list const &points() const { return d_; }

      /// Interpolate.
      /// \param x  Value of independent variable.
      /// \return   Interpolated value of dependent variable.
      D operator()(I x) const
      {
         if (d_.size() == 0) {
            throw "interpolating on empty interpolant";
         }
         if (x <= d_.begin()->first) {
            return d_.begin()->second;
         }
         if (x >= d_.rbegin()->first) {
            return d_.rbegin()->second;
         }
         using namespace std;
         point const xp{x, D()}; // Dummy y-coord used for search.
         auto const j = upper_bound(d_.begin(), d_.end(), xp, cmpr_pts);
         auto const i = j - 1;
         I const &xi = i->first;
         I const &xj = j->first;
         D const &yi = i->second;
         D const &yj = j->second;
         return yi + (yj - yi) * ((x - xi) / (xj - xi));
      }

      /// Integrate interpolant.
      /// \param x1  Beginning of interval of integration.
      /// \param x2  End of interval of integration.
      /// \return    Numerical result of definite integral.
      A integral(I x1, I x2) const
      {
         if (d_.size() == 0) {
            throw "integrating on empty interpolant";
         }
         double sign = 1.0;
         if (x1 > x2) {
            std::swap(x1, x2);
            sign = -1.0;
         }
         using R = decltype(I() * D());
         R sum(0.0 * x1 * d_[0].second);
         point const x1p{x1, D()}; // Dummy y-coord for search.
         auto i = upper_bound(d_.begin(), d_.end(), x1p, cmpr_pts);
         if (i == d_.end()) {
            // We need only calculate one piece to the right of all points.
            return sign * (x2 - x1) * d_.rbegin()->second;
         } else {
            // Handle first piece up to x coord of stored point.
            D const y1 = (*this)(x1);
            if (x2 <= i->first) {
               // We need calculate only one piece.
               return 0.5 * sign * (x2 - x1) * (y1 + (*this)(x2));
            }
            sum += 0.5 * (i->first - x1) * (y1 + i->second);
         }
         auto j = i + 1;
         while (j != d_.end() && x2 >= j->first) {
            // Handle piece between current point and next point.
            sum += 0.5 * (j->first - i->first) * (i->second + j->second);
            ++j;
            i = j - 1;
         }
         // Handle last piece.
         sum += 0.5 * (x2 - i->first) * (i->second + (*this)(x2));
         return sign * sum;
      }

      /// \return  Definite integral of interpolant between first and last
      ///          control point.
      A integral() const { return integral_; }

      /// Multiply the present instance I1 by another interpolant I2 such that
      /// there is a point in the product interpolant I3 at every unique x
      /// coordinate across I1 and I2.  Each new y coordinate is the
      /// product of a y coordinate, in either I1 or I2, and the
      /// corresponding interpolated y coordinate in the other.
      ///
      /// \tparam OD    Type of y coordinate in other interpolant.
      /// \param  oint  Reference to other interpolant.
      /// \return       Product interpolant.
      template <typename OD>
      auto operator*(interpolant<I, OD> const &oint) const
            -> interpolant<I, decltype(D() * OD())>
      {
         using D3 = decltype(D() * OD());
         std::function<D3(D, OD)> f = [](D y1, OD y2) { return y1 * y2; };
         return combine(oint, f);
      }

      /// Divide the present instance I1 by another interpolant I2 such that
      /// there is a point in the product interpolant I3 at every unique x
      /// coordinate across I1 and I2.  Each new y coordinate is the quotient
      /// of a y coordinate, in either I1 or I2, and the corresponding
      /// interpolated y coordinate in the other.
      ///
      /// \tparam OD    Type of y coordinate in other interpolant.
      /// \param  oint  Reference to other interpolant.
      /// \return       Quotient interpolant.
      template <typename OD>
      auto operator/(interpolant<I, OD> const &oint) const
            -> interpolant<I, decltype(D() * OD())>
      {
         using D3 = decltype(D() * OD());
         std::function<D3(D, OD)> f = [](D y1, OD y2) { return y1 / y2; };
         return combine(oint, f);
      }

      /// Modify every y-value by scalar multiplication.
      /// \return  Reference to present instance.
      interpolant &operator*=(double s /**< Scale factor. */)
      {
         for (unsigned i = 0; i < d_.size(); ++i) {
            d_[i].second *= s;
         }
         return *this;
      }

      /// Modify every y-value by scalar division.
      /// \return  Reference to present instance.
      interpolant &operator/=(double s /**< Inverse scale factor. */)
      {
         for (unsigned i = 0; i < d_.size(); ++i) {
            d_[i].second /= s;
         }
         return *this;
      }

      /// Multiply every y-value of interpolant by a scale factor on the
      /// right, and return the resultant interpolant.
      ///
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      auto operator*(double s) const -> interpolant<I, decltype(D() * s)>
      {
         using ND = decltype(D() * s);
         ilist<I, ND> nl(d_.size());
         for (unsigned i = 0; i < d_.size(); ++i) {
            nl[i] = ipoint<I, ND>(d_[i].first, d_[i].second * s);
         }
         return nl;
      }

      /// Multiply every y-value of interpolant by a scale factor on the left,
      /// and return the resultant interpolant.  Assume that multiplication of
      /// y-values is commutative.
      ///
      /// \param  s  Scale factor for y-values.
      /// \param  i  Interpolant to scale.
      /// \return    Scaled interpolant.
      friend auto operator*(double s, interpolant const &i)
            -> interpolant<I, decltype(D() * s)>
      {
         return i * s;
      }

      /// Multiply every y-value of interpolant by a statdim scale factor on
      /// the right, and return the resultant interpolant.
      ///
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      template <char TI, char ID, char M, char C, char TE>
      auto operator*(statdim<TI, ID, M, C, TE> s) const
            -> interpolant<I, decltype(D() * s)>
      {
         using ND = decltype(D() * s);
         ilist<I, ND> nl(d_.size());
         for (unsigned i = 0; i < d_.size(); ++i) {
            nl[i] = ipoint<I, ND>(d_[i].first, d_[i].second * s);
         }
         return nl;
      }

      /// Multiply every y-value of interpolant by a statdim scale factor on
      /// the left, and return the resultant interpolant.  Assume that
      /// multiplication of y-values is commutative.
      ///
      /// \param  s  Scale factor for y-values.
      /// \param  i  Interpolant to scale.
      /// \return    Scaled interpolant.
      template <char TI, char ID, char M, char C, char TE>
      friend auto operator*(statdim<TI, ID, M, C, TE> s, interpolant const &i)
            -> interpolant<I, decltype(D() * s)>
      {
         return i * s;
      }

      /// Divide every y-value of interpolant by a scale factor, and return the
      /// resultant interpolant.
      ///
      /// \tparam Y  Type of scale factor for y-values.
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      template <typename Y>
      auto operator/(Y s) const -> interpolant<I, decltype(D() / s)>
      {
         using ND = decltype(D() / s);
         ilist<I, ND> nl(d_.size());
         for (unsigned i = 0; i < d_.size(); ++i) {
            nl[i] = ipoint<I, ND>(d_[i].first, d_[i].second / s);
         }
         return nl;
      }
   };

   /// Alias for interpolant and maps from double to double.
   using interpolantd = interpolant<double, double>;

   /// Alias for initializer for interpolantd.
   using ilistd = ilist<double, double>;
}

#endif // ndef NUMERIC_INTERPOLANT_HPP

