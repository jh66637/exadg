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

#ifndef EXADG_AERO_ACOUSTIC_USER_INTERFACE_INPUT_PARAMETERS_H_
#define EXADG_AERO_ACOUSTIC_USER_INTERFACE_INPUT_PARAMETERS_H_

// deal.II
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>

namespace ExaDG
{
namespace AeroAcoustic
{
enum class FluidToAcousticCouplingStrategy
{
  Undefined,
  NonNestedMGRestriction
};

class Parameters
{
public:
  Parameters()
    : density(-1.0),
      source_term_with_convection(false),
      fluid_to_acoustic_coupling_strategy(FluidToAcousticCouplingStrategy::Undefined),
      use_analytical_source_term(false),
      use_analytical_cfd_solution(false)
  {
  }

  void
  check() const
  {
    AssertThrow(density >= 0.0, dealii::ExcMessage("Density has to be set."));

    AssertThrow(fluid_to_acoustic_coupling_strategy != FluidToAcousticCouplingStrategy::Undefined,
                dealii::ExcMessage("Coupling strategy has to be set."));

    AssertThrow(not(use_analytical_source_term and use_analytical_cfd_solution),
                dealii::ExcMessage(
                  "only use_analytical_source_term OR use_analytical_cfd_solution can be true."));
  }

  void
  print(dealii::ConditionalOStream const & pcout, std::string const & name) const
  {
    pcout << std::endl << name << std::endl << std::endl;
    print_parameter(pcout, "Density", density);
    print_parameter(pcout, "Source term has convective part", source_term_with_convection);
    print_parameter(pcout, "Fluid to acoustic coupling", fluid_to_acoustic_coupling_strategy);
    if(use_analytical_source_term)
      print_parameter(pcout, "Use analytical source term", use_analytical_source_term);
    if(use_analytical_cfd_solution)
      print_parameter(pcout, "Use analytical CFD solution", use_analytical_cfd_solution);
  }

  void
  add_parameters(dealii::ParameterHandler & prm, std::string const & subsection_name)
  {
    prm.enter_subsection(subsection_name);
    {
      prm.add_parameter(
        "Density", density, "Mean density of underlying fluid.", dealii::Patterns::Double(), true);

      prm.add_parameter("SourceTermWithConvection",
                        source_term_with_convection,
                        "Source term includes convective part.",
                        dealii::Patterns::Bool(),
                        true);

      prm.add_parameter("FluidToAcousticCouplingStrategy",
                        fluid_to_acoustic_coupling_strategy,
                        "Volume coupling strategy from the fluid to the acoustic field.",
                        Patterns::Enum<FluidToAcousticCouplingStrategy>(),
                        true);

      prm.add_parameter("UseAnalyticalSourceTerm",
                        use_analytical_source_term,
                        "Use an analytical representation of the aeroacoustic source term.",
                        dealii::Patterns::Bool());

      prm.add_parameter("UseAnalyticalCFDSolution",
                        use_analytical_cfd_solution,
                        "Compute the aeroacoustic source term from an analytical CFD solution.",
                        dealii::Patterns::Bool());
    }
    prm.leave_subsection();
  }

  // mean density of underlying fluid
  double density;

  // The aero-acoustic source term is the material derivative of the
  // pressure. Sometimes, it is sufficient to neglect the convective
  // part of the material derivative.
  bool source_term_with_convection;

  // Strategy to couple from fluid to acoustic
  FluidToAcousticCouplingStrategy fluid_to_acoustic_coupling_strategy;

  // In case an analytical representation of the aeroacoustic source
  // term exists we can use it to benchmark the aeroacoustic solver.
  bool use_analytical_source_term;

  // In case an analytical solution of the CFD is knwon, the aeroacoustic
  // source terms can be computed from this solution and it is not necessary
  // to solve the CFD.
  bool use_analytical_cfd_solution;
};

} // namespace AeroAcoustic
} // namespace ExaDG

#endif /* EXADG_AERO_ACOUSTIC_USER_INTERFACE_INPUT_PARAMETERS_H_ */
