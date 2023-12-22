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
    : source_term_with_convection(false),
      fluid_to_acoustic_coupling_strategy(FluidToAcousticCouplingStrategy::Undefined),
      blend_in_source_term(false),
      blend_in_duration(-1.0)
  {
  }

  void
  check() const
  {
    AssertThrow(fluid_to_acoustic_coupling_strategy != FluidToAcousticCouplingStrategy::Undefined,
                dealii::ExcMessage("Coupling strategy has to be set."));

    if(blend_in_source_term)
      AssertThrow(blend_in_duration > 0.0, dealii::ExcMessage("Blend in duration has to be set."));
  }

  void
  print(dealii::ConditionalOStream const & pcout, std::string const & name) const
  {
    pcout << std::endl << name << std::endl << std::endl;
    print_parameter(pcout, "Source term has convective part", source_term_with_convection);
    print_parameter(pcout, "Fluid to acoustic coupling", fluid_to_acoustic_coupling_strategy);
    print_parameter(pcout, "Blend in source term", blend_in_source_term);
    if(blend_in_source_term)
      print_parameter(pcout, "Blend in duration", blend_in_duration);
  }

  void
  add_parameters(dealii::ParameterHandler & prm, std::string const & subsection_name)
  {
    prm.enter_subsection(subsection_name);
    {
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

      prm.add_parameter("BlendInSourceTerm",
                        blend_in_source_term,
                        "Blend in the aeroacoustic source term.",
                        dealii::Patterns::Bool(),
                        true);

      prm.add_parameter("BlendInDuration",
                        blend_in_duration,
                        "Blend in the aeroacoustic source term.",
                        dealii::Patterns::Double());
    }
    prm.leave_subsection();
  }

  // The aero-acoustic source term is the material derivative of the
  // pressure. Sometimes, it is sufficient to neglect the convective
  // part of the material derivative.
  bool source_term_with_convection;

  // Strategy to couple from fluid to acoustic
  FluidToAcousticCouplingStrategy fluid_to_acoustic_coupling_strategy;

  // Blend in the aeroacoustic source term.
  bool blend_in_source_term;

  // Duration of source term blend in
  double blend_in_duration;
};

} // namespace AeroAcoustic
} // namespace ExaDG

#endif /* EXADG_AERO_ACOUSTIC_USER_INTERFACE_INPUT_PARAMETERS_H_ */
