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

#ifndef INCLUDE_EXADG_TIME_INTEGRATION_TIME_INTEGRATION_CONSTANTS_BASE_H_
#define INCLUDE_EXADG_TIME_INTEGRATION_TIME_INTEGRATION_CONSTANTS_BASE_H_

#include <vector>

namespace ExaDG
{
class TimeIntegratorConstantsBase
{
public:
  TimeIntegratorConstantsBase(unsigned int const order, bool const start_with_low_order)
    : order(order), start_with_low_order(start_with_low_order)
  {
  }

  virtual ~TimeIntegratorConstantsBase()
  {
  }

  /*
   *  This function updates the time integrator constants. The argument time_steps is only used in
   * case of adaptive time stepping.
   */
  void
  update(unsigned int const          current_order,
         bool const                  adaptive_time_stepping,
         std::vector<double> const & time_steps)
  {
    // when starting the time integrator with a low order method, ensure that
    // the time integrator constants are set properly
    unsigned int const update_order =
      (current_order <= order && start_with_low_order == true) ? current_order : order;

    if(adaptive_time_stepping)
      set_adaptive_time_step(update_order, time_steps);
    else
      set_constant_time_step(update_order);
  }

  unsigned int
  get_order() const
  {
    return order;
  }

  /*
   *  This function prints the time integrator constants
   */
  virtual void
  print(dealii::ConditionalOStream & pcout) const = 0;

protected:
  /**
   * Can be used to fill the unused components of the vector with zeros.
   * This is needed since current_order might be smaller than order, e.g.,
   * when using start_with_low_order = true
   */
  void
  zero_out_unused_constants(unsigned int const first_idx, std::vector<double> & constants)
  {
    for(unsigned int i = first_idx; i < constants.size(); ++i)
      constants[i] = 0.0;
  }


  // order of time integrator
  unsigned int const order;

  // use a low order time integration scheme to start the time integrator?
  bool const start_with_low_order;

private:
  /*
   *  This function calculates the time integrator constants in case of constant time step sizes.
   */
  virtual void
  set_constant_time_step(unsigned int const current_order) = 0;

  /*
   *  This function calculates time integrator constants in case of varying time step sizes
   * (adaptive time stepping).
   */
  virtual void
  set_adaptive_time_step(unsigned int const          current_order,
                         std::vector<double> const & time_steps) = 0;
};
} // namespace ExaDG


#endif /* INCLUDE_EXADG_TIME_INTEGRATION_TIME_INTEGRATION_CONSTANTS_BASE_H_ */
