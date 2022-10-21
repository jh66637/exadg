/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2021 by the ExaDG authors
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

// ExaDG
#include <exadg/postprocessor/time_control_statistics.h>
#include <exadg/utilities/print_functions.h>

namespace ExaDG
{
TimeControlDataStatistics::TimeControlDataStatistics()
  : write_preliminary_results_every_nth_time_step(numbers::invalid_timestep)
{
}

void
TimeControlDataStatistics::print(dealii::ConditionalOStream & pcout, bool const unsteady) const
{
  time_control_data.print(pcout, unsteady);
  if(Utilities::is_valid_timestep(write_preliminary_results_every_nth_time_step))
    print_parameter(pcout,
                    "Write preliminary results every nth time step",
                    write_preliminary_results_every_nth_time_step);
}

void
TimeControlStatistics::setup(TimeControlDataStatistics const & time_control_data_statistics_in)
{
  time_control_data_statistics = time_control_data_statistics_in;
  time_control.setup(time_control_data_statistics.time_control_data);
}

bool
TimeControlStatistics::needs_evaluation(double const           time,
                                        types::time_step const time_step_number) const
{
  if(time_control.needs_evaluation(time, time_step_number))
    return true;
  if(write_preliminary_results(time_step_number))
    return true;

  return false;
}

unsigned int
TimeControlStatistics::get_counter() const
{
  return time_control.get_counter();
}

bool
TimeControlStatistics::write_preliminary_results(types::time_step const time_step_number) const
{
  if(Utilities::is_valid_timestep(
       time_control_data_statistics.write_preliminary_results_every_nth_time_step))
    if((time_step_number - 1) %
           (time_control_data_statistics.write_preliminary_results_every_nth_time_step) ==
         0 ||
       time_control.reached_end_time())
      return true;

  return false;
}


} // namespace ExaDG
