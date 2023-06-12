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

#ifndef INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_ABM_BASE_H_
#define INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_ABM_BASE_H_

// ExaDG
#include <exadg/time_integration/ab_constants.h>
#include <exadg/time_integration/am_constants.h>
#include <exadg/time_integration/bdf_constants.h>
#include <exadg/time_integration/time_int_multistep_base.h>

namespace ExaDG
{
class TimeIntABMBase : public TimeIntMultistepBase
{
public:
  TimeIntABMBase(double const        start_time_,
                 double const        end_time_,
                 unsigned int const  max_number_of_time_steps_,
                 unsigned int const  order_ab,
                 unsigned int const  order_am,
                 bool const          start_with_low_order_,
                 bool const          adaptive_time_stepping_,
                 RestartData const & restart_data_,
                 MPI_Comm const &    mpi_comm_,
                 bool const          is_test_)
    : TimeIntMultistepBase(start_time_,
                           end_time_,
                           max_number_of_time_steps_,
                           order_,
                           start_with_low_order_,
                           adaptive_time_stepping_,
                           restart_data_,
                           mpi_comm_,
                           is_test_),
      ab(order_ab, start_with_low_order_),
      am(order_am, start_with_low_order_),
      bdf(order_, start_with_low_order_)
  {
  }

protected:
  void
  update_time_integrator_constants() override
  {
    ab.update(time_step_number, adaptive_time_stepping, time_steps);
    am.update(time_step_number, adaptive_time_stepping, time_steps);
    bdf.update(time_step_number, adaptive_time_stepping, time_steps);
  }

  ABTimeIntegratorConstants ab;

  AMTimeIntegratorConstants am;
  // needed for ALE
  BDFTimeIntegratorConstants bdf;
};

} // namespace ExaDG

#endif /* INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_AMB_BASE_H_ */
