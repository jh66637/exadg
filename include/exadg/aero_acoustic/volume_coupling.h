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

#ifndef INCLUDE_EXADG_AERO_ACOUSTIC_VOLUME_COUPLING_H_
#define INCLUDE_EXADG_AERO_ACOUSTIC_VOLUME_COUPLING_H_

#include <exadg/aero_acoustic/calculators/source_term_calculator.h>
#include <exadg/aero_acoustic/single_field_solvers/acoustics.h>
#include <exadg/aero_acoustic/single_field_solvers/fluid.h>
#include <exadg/aero_acoustic/user_interface/parameters.h>

namespace ExaDG
{
namespace AeroAcoustic
{
/**
 * A class that handles the volume coupling between fluid and acoustic.
 */
template<int dim, typename Number>
class VolumeCoupling
{
  using VectorType = dealii::LinearAlgebra::distributed::Vector<Number>;

public:
  void
  setup(Parameters const &                           parameters_in,
        std::shared_ptr<SolverAcoustic<dim, Number>> acoustic_solver_in,
        std::shared_ptr<SolverFluid<dim, Number>>    fluid_solver_in,
        std::shared_ptr<FieldFunctions<dim>>         field_functions_in)
  {
    parameters      = parameters_in;
    acoustic_solver = acoustic_solver_in;
    fluid_solver    = fluid_solver_in;

    acoustic_solver_in->pde_operator->initialize_dof_vector_pressure(source_term_acoustic);
    fluid_solver_in->pde_operator->initialize_vector_pressure(source_term_fluid);

    // setup the transfer operator
    if(parameters.fluid_to_acoustic_coupling_strategy ==
         FluidToAcousticCouplingStrategy::ConservativeInterpolation or
       parameters.acoustic_to_fluid_coupling_strategy ==
         FluidToAcousticCouplingStrategy::NonNestedGridTransfer)
    {
      non_nested_grid_transfer.reinit(fluid_solver_in->pde_operator->get_dof_handler_p(),
                                      acoustic_solver_in->pde_operator->get_dof_handler_p(),
                                      *fluid_solver_in->pde_operator->get_mapping(),
                                      *acoustic_solver_in->pde_operator->get_mapping());
    }
    else
    {
      AssertThrow(false, dealii::ExcMessage("FluidToAcousticCouplingStrategy not implemented."));
    }

    // setup aeroacoustic source term calculator
    SourceTermCalculatorData<dim> data;
    data.dof_index_pressure  = fluid_solver_in->pde_operator->get_dof_index_pressure();
    data.dof_index_velocity  = fluid_solver_in->pde_operator->get_dof_index_velocity();
    data.quad_index          = fluid_solver_in->pde_operator->get_quad_index_pressure();
    data.density             = parameters_in.density;
    data.consider_convection = parameters_in.source_term_with_convection;
    data.blend_in            = parameters.blend_in_source_term;
    data.blend_in_function   = field_functions_in->source_term_blend_in;

    source_term_calculator.setup(fluid_solver_in->pde_operator->get_matrix_free(), data);
  }

  void
  fluid_to_acoustic()
  {
    if(parameters.fluid_to_acoustic_coupling_strategy ==
       FluidToAcousticCouplingStrategy::ConservativeInterpolation)
    {
      source_term_calculator.evaluate_integrate(source_term_fluid,
                                                fluid_solver->time_integrator->get_velocity(),
                                                fluid_solver->time_integrator->get_pressure(),
                                                fluid_solver->get_pressure_time_derivative(),
                                                fluid_solver->time_integrator->get_time());

      non_nested_grid_transfer.restrict_and_add(source_term_acoustic, source_term_fluid);
    }
    else
    {
      AssertThrow(false, dealii::ExcMessage("FluidToAcousticCouplingStrategy not implemented."));
    }

    acoustic_solver->pde_operator->set_aero_acoustic_source_term(source_term_acoustic);
  }

  void
  acoustic_to_fluid()
  {
    if(parameters.acoustic_to_fluid_coupling_strategy ==
       FluidToAcousticCouplingStrategy::NonNestedGridTransfer)
    {
      // TODO: we somehow have to transfer u_a to the fluid mesh. for this we need a different
      // grid transfer class and different integration rules because we are working with vecotrs
      // Hence we should split the class!!
      non_nested_grid_transfer.prolongate_and_add(acoustic_particle_velocity_on_fluid_mesh,
                                                  acoustic_solver->time_integrator->get_velocity());

      feedback_term_calculator.evaluate_integrate(feedback_term_fluid,
                                                  fluid_solver->time_integrator->get_velocity(),
                                                  acoustic_particle_velocity_on_fluid_mesh);
    }
    else
    {
      AssertThrow(false, dealii::ExcMessage("FluidToAcousticCouplingStrategy not implemented."));
    }

    fluid_solver->pde_operator->set_aero_acoustic_feedback_term(feedback_term_fluid);
  }

private:
  Parameters parameters;

  // Single field solvers
  std::shared_ptr<SolverAcoustic<dim, Number>> acoustic_solver;
  std::shared_ptr<SolverFluid<dim, Number>>    fluid_solver;

  // Field functions
  std::shared_ptr<FieldFunctions<dim>> field_functions;

  // Transfer operator
  dealii::MGTwoLevelTransferNonNested<dim, VectorType> non_nested_grid_transfer;

  // Class that knows how to compute the source term
  SourceTermCalculator<dim, Number> source_term_calculator;

  // Class that knows how to compute the feedback term
  SourceTermCalculator<dim, Number> feedback_term_calculator;

  // Aeroacoustic source and feedback terms defined on the acoustic mesh
  VectorType source_term_acoustic;
  VectorType feedback_term_acoustic;

  // Aeroacoustic source and feedback terms defined on the fluid mesh
  VectorType source_term_fluid;
  VectorType feedback_term_fluid;
};
} // namespace AeroAcoustic
} // namespace ExaDG

#endif /*INCLUDE_EXADG_AERO_ACOUSTIC_VOLUME_COUPLING_H_*/
