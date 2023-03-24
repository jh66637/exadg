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

#ifndef INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_
#define INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_

#include <deal.II/lac/la_parallel_vector.h>

#include <exadg/postprocessor/pointwise_output_generator_base.h>

namespace ExaDG
{
namespace CompNS
{
template<int dim>
struct PointwiseOutputData : public PointwiseOutputDataBase<dim>
{
  PointwiseOutputData() : write_rho(false), write_rho_u(false), write_rho_E(false)
  {
  }

  void
  print(dealii::ConditionalOStream & pcout) const
  {
    PointwiseOutputDataBase<dim>::print(pcout);

    if(this->time_control_data.is_active and this->evaluation_points.size() > 0)
    {
      print_parameter(pcout, "Write rho", write_rho);
      print_parameter(pcout, "Write rho_u", write_rho_u);
      print_parameter(pcout, "Write rho_E", write_rho_E);
    }
  }

  bool write_rho;
  bool write_rho_u;
  bool write_rho_E;
};

template<int dim, typename VectorType>
class PointwiseOutputGenerator : public PointwiseOutputGeneratorBase<dim, VectorType>
{
public:
  PointwiseOutputGenerator(MPI_Comm const & comm)
    : PointwiseOutputGeneratorBase<dim, VectorType>(comm)
  {
  }

  void
  setup(dealii::DoFHandler<dim> const &  dof_handler_in,
        dealii::Mapping<dim> const &     mapping_in,
        PointwiseOutputData<dim> const & pointwise_output_data_in)
  {
    this->setup_base(dof_handler_in.get_triangulation(), mapping_in, pointwise_output_data_in);

    dof_handler           = &dof_handler_in;
    pointwise_output_data = pointwise_output_data_in;

    if(pointwise_output_data.write_rho)
      this->add_quantity("Rho", 1);
    if(pointwise_output_data.write_rho_u)
      this->add_quantity("Rho_U", dim);
    if(pointwise_output_data.write_rho_E)
      this->add_quantity("Rho_E", 1);
  }

private:
  void
  do_evaluate(VectorType const & solution) final
  {
    if(pointwise_output_data.write_rho or pointwise_output_data.write_rho_u or
       pointwise_output_data.write_rho_E)
    {
      auto const values = this->template compute_point_values<dim + 2>(solution, *dof_handler);
      if(pointwise_output_data.write_rho)
        this->write_quantity("Rho", values, 0);
      if(pointwise_output_data.write_rho_u)
        this->write_quantity("Rho_U", values, 1);
      if(pointwise_output_data.write_rho_E)
        this->write_quantity("Rho_E", values, 1 + dim);
    }
  }

  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler;
  PointwiseOutputData<dim>                            pointwise_output_data;
};

} // namespace CompNS
} // namespace ExaDG


#endif /* INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_ \
        */
