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
// Container that collects pressure and velocity. This is only used to pass both quantities to the
// base class
template<typename Number>
struct SolutionContainer
{
  using VectorType = dealii::LinearAlgebra::distributed::Vector<Number>;
  using value_type = Number;

  SolutionContainer(VectorType const & velocity, VectorType const & pressure)
    : velocity(velocity), pressure(pressure)
  {
  }

  VectorType const & velocity;
  VectorType const & pressure;
};


template<int dim>
struct PointwiseOutputData : public PointwiseOutputDataBase<dim>
{
  PointwiseOutputData();

  void
  print(dealii::ConditionalOStream & pcout) const;

  bool write_pressure;
  bool write_velocity;
};

template<int dim, typename Number>
class PointwiseOutputGenerator : public PointwiseOutputGeneratorBase<dim, SolutionContainer<Number>>
{
  using VectorType = typename SolutionContainer<Number>::VectorType;

public:
  PointwiseOutputGenerator(MPI_Comm const & comm);

  void
  setup(dealii::DoFHandler<dim> const &  dof_handler_velocity_in,
        dealii::DoFHandler<dim> const &  dof_handler_pressure_in,
        dealii::Mapping<dim> const &     mapping_in,
        PointwiseOutputData<dim> const & pointwise_output_data_in);

private:
  void
  do_evaluate(SolutionContainer<Number> const & solution) final;

  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler_pressure;
  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler_velocity;

  PointwiseOutputData<dim> pointwise_output_data;
};

} // namespace IncNS
} // namespace ExaDG


#endif /* INCLUDE_EXADG_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_H_ \
        */
