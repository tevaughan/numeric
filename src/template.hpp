
/// \file   template.hpp
/// \brief  Template-related utilities.

#ifndef NUMERIC_TEMPLATE_HPP
#define NUMERIC_TEMPLATE_HPP

namespace num
{
   /// Type of ratio of two tings.
   /// \tparam Y  Type of numerator.
   /// \tparam X  Type of denominator.
   template <typename Y, typename X>
   using RAT = decltype(Y() / X());

   /// Type of product of two things.
   /// \tparam X  Type of left factor.
   /// \tparam Y  Type of right factor.
   template <typename X, typename Y>
   using PRD = decltype(X() * Y());
}

#endif // ndef NUMERIC_TEMPLATE_HPP

