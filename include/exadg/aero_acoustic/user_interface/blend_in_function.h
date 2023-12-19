/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2023 by the ExaDG authors
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

#ifndef INCLUDE_EXADG_AERO_ACOUSTIC_USER_INTERFACE_BLEND_IN_FUNCTION_H_
#define INCLUDE_EXADG_AERO_ACOUSTIC_USER_INTERFACE_BLEND_IN_FUNCTION_H_

namespace ExaDG
{
namespace AeroAcoustic
{
class BlendInFunction
{
public:
  BlendInFunction() : start_time(-1.0), end_time(-1.0)
  {
  }

  virtual ~BlendInFunction() = default;

  void
  setup(double const start_time_in, double const end_time_in)
  {
    start_time = start_time_in;
    end_time   = end_time_in;
  }

  bool
  is_finished(double const time) const
  {
    if(time < end_time)
      return false;
    return true;
  }

  double
  get_scaling_factor(double const time) const
  {
    return value(std::make_pair(start_time, end_time), time);
  }

private:
  virtual double
  value(std::pair<double, double> const & start_end_time, double const time) const = 0;

  double start_time;
  double end_time;
};

} // namespace AeroAcoustic
} // namespace ExaDG


#endif /* INCLUDE_EXADG_AERO_ACOUSTIC_USER_INTERFACE_BLEND_IN_FUNCTION_H_ */
