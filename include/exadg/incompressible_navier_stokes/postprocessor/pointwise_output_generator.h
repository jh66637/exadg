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

#ifndef INCLUDE_EXADG_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_
#define INCLUDE_EXADG_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_

#include <deal.II/lac/la_parallel_vector.h>

#include <exadg/postprocessor/pointwise_output_generator_base.h>

namespace ExaDG
{
namespace IncNS
{
template<int dim>
struct PointwiseOutputData : public PointwiseOutputDataBase<dim>
{
  PointwiseOutputData() : write_velocity(false), write_pressure(false)
  {
  }

  void
  print(dealii::ConditionalOStream & pcout) const
  {
    PointwiseOutputDataBase<dim>::print(pcout);

    if(this->time_control_data.is_active && this->evaluation_points.size() > 0)
    {
      print_parameter(pcout, "Write velocity", write_velocity);
      print_parameter(pcout, "Write pressure", write_pressure);
    }
  }

  bool write_velocity;
  bool write_pressure;
};

template<int dim, typename VectorType>
class PointwiseOutputGenerator
  : public PointwiseOutputGeneratorBase<dim,
                                        std::vector<VectorType const *>,
                                        typename VectorType::value_type>
{
public:
  PointwiseOutputGenerator(MPI_Comm const & comm)
    : PointwiseOutputGeneratorBase<dim,
                                   std::vector<VectorType const *>,
                                   typename VectorType::value_type>(comm)
  {
  }

  void
  setup(dealii::DoFHandler<dim> const &  dof_handler_velocity_in,
        dealii::DoFHandler<dim> const &  dof_handler_pressure_in,
        dealii::Mapping<dim> const &     mapping_in,
        PointwiseOutputData<dim> const & pointwise_output_data_in)
  {
    this->setup_base(dof_handler_pressure_in.get_triangulation(),
                     mapping_in,
                     pointwise_output_data_in);

    dof_handler_velocity = &dof_handler_velocity_in;
    dof_handler_pressure = &dof_handler_pressure_in;

    pointwise_output_data = pointwise_output_data_in;

    if(pointwise_output_data.write_velocity)
      this->add_quantity("Velocity", dim);
    if(pointwise_output_data.write_pressure)
      this->add_quantity("Pressure", 1);
  }

private:
  void
  do_evaluate(std::vector<VectorType const *> const & solution) final
  {
    if(pointwise_output_data.write_velocity)
    {
      auto const values =
        this->template compute_point_values<dim>(*solution[0], *dof_handler_velocity);
      this->write_quantity("Velocity", values, 0);
    }
    if(pointwise_output_data.write_pressure)
    {
      auto const values =
        this->template compute_point_values<1>(*solution[1], *dof_handler_pressure);
      this->write_quantity("Pressure", values);
    }
  }


  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler_velocity;
  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler_pressure;

  PointwiseOutputData<dim> pointwise_output_data;
};

} // namespace IncNS
} // namespace ExaDG


#endif /* INCLUDE_EXADG_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_ \
        */
