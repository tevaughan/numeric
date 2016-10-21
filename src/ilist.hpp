
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   ilist.hpp
///
/// \brief  Definition for each of num::ipoint, num::ilist, corresponding
///         output-stream operators, num::get_point(), num::get_points(), and
///         num::midpoints().

#ifndef NUMERIC_ILIST_HPP
#define NUMERIC_ILIST_HPP

#include <fstream>  // for ifstream
#include <iostream> // for ostream
#include <sstream>  // for istringstream
#include <vector>   // for vector

namespace num
{
   /// Point used as one of the constraints of a linear interpolant.
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   using ipoint = std::pair<I, D>;

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
   using ilist = std::vector<ipoint<I, D>>;

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

   /// Extract point from space-delimited ASCII line.
   /// \tparam X  Type of first column, representing X coordinate.
   /// \tparam Y  Type of second column, representing Y coordinate.
   template <typename X, typename Y>
   ipoint<X, Y> get_point(
         /** Line of non-blank ASCII text. */ std::string line,
         /** Unit to multiply against first column. */ X  xu,
         /** Unit to multiply against secnd column. */ Y  yu)
   {
      std::istringstream iss(line);
      double             x, y;
      if (!(iss >> x)) {
         throw "error reading x";
      }
      if (!(iss >> y)) {
         throw "error reading y";
      }
      return ipoint<X, Y>(x * xu, y * yu);
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
   ilist<X, Y> get_points(
         /** Name of ASCII file.          */ std::string file,
         /** Unit to multiply against first column. */ X xu,
         /** Unit to multiply against secnd column. */ Y yu)
   {
      std::ifstream is(file);
      if (!is) {
         throw "Failed to open '" + file + "'.";
      }
      std::string line;
      ilist<X, Y> points;
      while (getline(is, line)) {
         size_t const p = line.find_first_not_of(" \f\n\r\t");
         if (p == std::string::npos || line[p] == '#') {
            continue;
         }
         points.push_back(get_point(line, xu, yu));
      }
      auto comp = [](ipoint<X, Y> const &a, ipoint<X, Y> const &b) {
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
         unsigned const i  = j - 1;
         auto const &   xi = v[i].first, &yi = v[i].second;
         auto const &   xj = v[j].first, &yj = v[j].second;
         rv[i] = {0.5 * (xi + xj), 0.5 * (yi + yj)};
      }
      return rv;
   }
}

#endif // ndef NUMERIC_ILIST_HPP

