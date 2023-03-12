/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2022 by the ExaDG authors
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *  ______________________________________________________________________
 */

#ifndef INCLUDE_EXADG_UTILITIES_ENUM_PATTERNS_H_
#define INCLUDE_EXADG_UTILITIES_ENUM_PATTERNS_H_

// inspired by following code which uses different backend for enums
// https://github.com/peterrum/dealii-parameter-handler-enum/blob/master/include/deal.II/base/patterns_enum.h

#include <deal.II/base/patterns.h>

#include <exadg/utilities/enum_utilities.h>

#include <string>

namespace dealii
{
namespace Patterns
{
namespace Tools
{
/**
 * Converter class for structs of enum (scoped or unscoped) type.
 * dealii::Patterns::Selection is automatically generated from the given enum type and
 * automatically checked. Outputting possible parameters is also possible.
 */
template<class T>
struct Convert<T, typename std::enable_if<ExaDG::Utilities::is_enum<T>()>::type>
{
  /**
   * Convert to pattern.
   */
  static std::unique_ptr<Patterns::PatternBase>
  to_pattern()
  {
    return std::make_unique<Patterns::Selection>(ExaDG::Utilities::serialized_string<T>());
  }

  /**
   * Convert value to string.
   */
  static std::string
  to_string(const T & t, const Patterns::PatternBase & = *Convert<T>::to_pattern())
  {
    return ExaDG::Utilities::enum_to_string(t);
  }

  /**
   * Convert string to value.
   */
  static T
  to_value(const std::string & s, const Patterns::PatternBase & = *Convert<T>::to_pattern())
  {
    T value;
    ExaDG::Utilities::string_to_enum(value, s);
    return value;
  }
};

} // namespace Tools
} // namespace Patterns
} // namespace dealii

#endif /*INCLUDE_EXADG_UTILITIES_ENUM_PATTERNS_H_*/
