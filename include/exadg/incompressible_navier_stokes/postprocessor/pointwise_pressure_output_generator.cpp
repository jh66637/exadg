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

#include <deal.II/numerics/vector_tools.h>

#include <exadg/utilities/print_functions.h>

#include <exadg/incompressible_navier_stokes/postprocessor/pointwise_output_generator.h>

namespace ExaDG
{
namespace IncNS
{
template<int dim>
PointwisePressureOutputData<dim>::PointwisePressureOutputData()
  : write_pressure(false), write_velocity(false)
{
}

template<int dim>
void
PointwisePressureOutputData<dim>::print(dealii::ConditionalOStream & pcout) const
{
  PointwiseOutputDataBase<dim>::print(pcout);

  if(this->time_control_data.is_active && this->evaluation_points.size() > 0)
  {
    print_parameter(pcout, "Write pressure", write_pressure);
  }
}

template class PointwisePressureOutputData<2>;
template class PointwisePressureOutputData<3>;

template<int dim, typename Number>
PointwisePressureOutputGenerator<dim, Number>::PointwisePressureOutputGenerator(
  MPI_Comm const & comm)
  : PointwisePressureOutputGeneratorBase<dim, VectorType>(comm)
{
}

template<int dim, typename Number>
void
PointwisePressureOutputGenerator<dim, Number>::setup(
  dealii::DoFHandler<dim> const &  dof_handler_in,
  dealii::Mapping<dim> const &     mapping_in,
  PointwiseOutputData<dim> const & pointwise_output_data_in)
{
  this->setup_base(dof_handler_velocity_in.get_triangulation(),
                   mapping_in,
                   pointwise_output_data_in);

  dof_handler = &dof_handler_in;

  pointwise_output_data = pointwise_output_data_in;

  if(pointwise_output_data.write_pressure)
    this->add_quantity("Pressure", 1);
}

template<int dim, typename Number>
void
PointwisePressureOutputGenerator<dim, Number>::do_evaluate(BlockVectorType const & pressure)
{
  if(pointwise_output_data.write_pressure)
  {
    auto const values = this->template compute_point_values<1>(pressure, *dof_handler);
    this->write_quantity("Pressure", values);
  }
}

template class PointwisePressureOutputGenerator<2, float>;
template class PointwisePressureOutputGenerator<2, double>;

template class PointwisePressureOutputGenerator<3, float>;
template class PointwisePressureOutputGenerator<3, double>;

} // namespace IncNS
} // namespace ExaDG
