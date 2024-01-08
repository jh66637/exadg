/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2024 by the ExaDG authors
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

/* #ifndef EXADG_AERO_ACOUSTIC_USER_INTERFACE_FIELD_FUNCTIONS_H_ */
/* #define EXADG_AERO_ACOUSTIC_USER_INTERFACE_FIELD_FUNCTIONS_H_ */

namespace ExaDG
{
namespace AeroAcoustic
{
template<int dim>
struct FieldFunctions
{
  /*
   * The function analytical_source_term to prescribe the aeroacoustic
   * source term analytically. If the analytical source term is used,
   * the incompressible Navier-Stokes equations do not have to be solved.
   */
  std::shared_ptr<dealii::Function<dim>> analytical_source_term;

  /*
   * Use the analytical CFD solutions to compute aeroacoustic source terms.
   * If the analytical CFD solution is known, the incompressible Navier-Stokes
   * equations do not have to be solved.
   */
  std::shared_ptr<dealii::Function<dim>> analytical_cfd_solution_velocity;
  std::shared_ptr<dealii::Function<dim>> analytical_cfd_solution_pressure;
};

} // namespace AeroAcoustic
} // namespace ExaDG

/* #endif /\* EXADG_AERO_ACOUSTIC_USER_INTERFACE_FIELD_FUNCTIONS_H_ *\/ */
