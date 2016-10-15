
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   ilist.hpp
/// \brief  Definition of ipoint, ilist, and output-stream operator for each.

#ifndef NUMERIC_ILIST_HPP
#define NUMERIC_ILIST_HPP

#include <iostream> // for ostream
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
}

#endif // ndef NUMERIC_ILIST_HPP

