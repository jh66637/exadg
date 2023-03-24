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


#include <exadg/postprocessor/pointwise_output_generator_base.h>
#include <exadg/utilities/print_functions.h>

namespace ExaDG
{
template<int dim>
PointwiseOutputDataBase<dim>::PointwiseOutputDataBase()
  : directory("output/"),
    filename("name"),
    update_points_before_evaluation(false),
    restarted_simulation(false)
{
}

template<int dim>
void
PointwiseOutputDataBase<dim>::print(dealii::ConditionalOStream & pcout) const
{
  if(time_control_data.is_active && evaluation_points.size() > 0)
  {
    pcout << std::endl << "Pointwise output" << std::endl;

    // this class makes only sense for the unsteady case
    time_control_data.print(pcout, true /*unsteady*/);

    print_parameter(pcout, "Output directory", directory);
    print_parameter(pcout, "Name of output file", filename);
    print_parameter(pcout, "Update Points before evaluation", update_points_before_evaluation);
  }
}

template struct PointwiseOutputDataBase<2>;
template struct PointwiseOutputDataBase<3>;
} // namespace ExaDG
