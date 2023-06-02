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

#ifndef INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_BDF_BASE_H_
#define INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_BDF_BASE_H_

// ExaDG
#include <exadg/time_integration/bdf_constants.h>
#include <exadg/time_integration/extrapolation_scheme.h>
#include <exadg/time_integration/time_int_multistep_base.h>

namespace ExaDG
{
class TimeIntBDFBase : public TimeIntMultistepBase
{
public:
  TimeIntBDFBase(double const        start_time_,
                 double const        end_time_,
                 unsigned int const  max_number_of_time_steps_,
                 unsigned const      order_,
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
      bdf(order_, start_with_low_order_),
      extra(order_, start_with_low_order_)
  {
  }

  double
  get_scaling_factor_time_derivative_term() const
  {
    return bdf.get_gamma0() / time_steps[0];
  }


  BDFTimeIntegratorConstants const &
  get_bdf_constants() const
  {
    return bdf;
  }

protected:
  void
  update_time_integrator_constants() override
  {
    bdf.update(time_step_number, adaptive_time_stepping, time_steps);
    extra.update(time_step_number, adaptive_time_stepping, time_steps);

    // use this function to check the correctness of the time integrator constants
    //  std::cout << std::endl << "Time step " << time_step_number << std::endl << std::endl;
    //  std::cout << "Coefficients BDF time integration scheme:" << std::endl;
    //  bdf.print();
    //  std::cout << "Coefficients extrapolation scheme:" << std::endl;
    //  extra.print();
  }

  /*
   * Time integration constants. The extrapolation scheme is not necessarily used for a BDF time
   * integration scheme with fully implicit time stepping, implying a violation of the Liskov
   * substitution principle (OO software design principle). However, it does not appear to be
   * reasonable to complicate the inheritance due to this fact.
   */
  BDFTimeIntegratorConstants bdf;
  ExtrapolationConstants     extra;
};
} // namespace ExaDG

#endif /* INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_BDF_BASE_H_ */
