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

#ifndef INCLUDE_EXADG_POSTPROCESSOR_OUTPUT_PARAMETERS_H_
#define INCLUDE_EXADG_POSTPROCESSOR_OUTPUT_PARAMETERS_H_

#include <exadg/utilities/parameter_handler.h>

namespace ExaDG
{
struct OutputParameters
{
  OutputParameters() : directory("output/"), filename("solution"), write(false)
  {
  }

  void
  add_parameters(dealii::ParameterHandler & prm, std::string const & subsection_name = "Output")
  {
    // clang-format off
    prm.enter_subsection(subsection_name);
      prm.add_parameter("OutputDirectory",  directory, "Directory where output is written.");
      prm.add_parameter("OutputName",       filename,  "Name of output files.");
      prm.add_parameter("WriteOutput",      write,     "Decides whether output is written.");
    prm.leave_subsection();
    // clang-format on
  }

  std::string directory;
  std::string filename;
  bool        write;
};
} // namespace ExaDG



#endif /* INCLUDE_EXADG_POSTPROCESSOR_OUTPUT_PARAMETERS_H_ */
