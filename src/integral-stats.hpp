
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

/// \file   integral-stats.hpp
///
/// \brief  Source for num::integral_stats, which is used by num::integral()
///         and by num::interpolant.

#ifndef NUMERIC_INTEGRAL_STATS_HPP
#define NUMERIC_INTEGRAL_STATS_HPP

namespace num
{
   /// Summary of statistics on trapezoids in integral.
   /// \tparam I  Type of integral, corresponding to "area" under curve.
   template <typename I>
   class integral_stats
   {
      /// Type of square deviation in area under trapezoid.
      using S = decltype(I() * I());

      unsigned num_; ///< Number of trapezoids.
      I area_;       ///< Total area of trapezoids.
      S sqdv_;       ///< Sum of estimated square deviations in area.

   public:
      /// Construct null statistical summary.
      integral_stats(I zero) : num_(0), area_(zero), sqdv_(zero * zero) {}

      /// Add a trapezoidal area and an estimated error in area.
      /// \param a  Trapezoidal area.
      /// \param d  Estimated error in area.
      void add(I const &a, I const &d)
      {
         ++num_;
         area_ += a;
         sqdv_ += d * d;
      }

      /// Total area of trapezoids.
      I const &area() const { return area_; }

      /// Standard deviation of estimated total area in trapezoids.
      I stdev() const { return sqrt(sqdv_ / num_); }
   };
}

#endif // ndef NUMERIC_INTEGRAL_STATS_HPP

