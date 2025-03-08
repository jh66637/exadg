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

#pragma once

#include <deal.II/lac/la_parallel_block_vector.h>

#include <exadg/matrix_free/integrators.h>
#include <exadg/operators/mapping_flags.h>

#include <exadg/acoustic_conservation_equations/user_interface/parameters.h>

namespace ExaDG::Acoustics
{
namespace PML::Utilities
{
template<int dim>
unsigned int
categorize_pml_cells(dealii::Triangulation<dim> const & tria, std::vector<unsigned int> & dst)
{
  unsigned int n_pml_cells = 0;
  dst.resize(tria.n_active_cells());
  for(const auto & cell : tria.active_cell_iterators())
  {
    if(cell->is_locally_owned())
    {
      AssertIndexRange(cell->active_cell_index(), tria.n_active_cells());

      if(cell->material_id() == numbers::pml_material_id)
      {
        dst[cell->active_cell_index()] = numbers::pml_material_id;
        ++n_pml_cells;
      }
      else
      {
        dst[cell->active_cell_index()] = 0;
      }
    }
  }
  return n_pml_cells;
}

} // namespace PML::Utilities

namespace Operators
{
template<int dim, typename Number>
class PMLKernel
{
public:
  static MappingFlags
  get_mapping_flags()
  {
    MappingFlags flags;
    flags.cells =
      dealii::update_JxW_values | dealii::update_gradients | dealii::update_quadrature_points;
    return flags;
  }
};
} // namespace Operators

template<unsigned int dim>
struct PMLOperatorData
{
  unsigned int dof_index_pressure;
  unsigned int dof_index_velocity;
  unsigned int quad_index;

  unsigned int block_index_pressure;
  unsigned int block_index_velocity;
  unsigned int block_index_auxiliary;

  // TODO: do we really need it or can we always have the same damping? but then we
  // have to know the normal direction of the pml?
  std::shared_ptr<dealii::Function<dim>> pml_damping;
};

template<int dim, typename Number>
class PMLOperator
{
public:
  using This = PMLOperator<dim, Number>;

  using BlockVectorType = dealii::LinearAlgebra::distributed::BlockVector<Number>;
  using VectorType      = dealii::LinearAlgebra::distributed::Vector<Number>;

  using scalar = dealii::VectorizedArray<Number>;
  using vector = dealii::Tensor<1, dim, dealii::VectorizedArray<Number>>;
  using tensor = dealii::Tensor<2, dim, dealii::VectorizedArray<Number>>;

  using Range = std::pair<unsigned int, unsigned int>;

  using CellIntegratorU = CellIntegrator<dim, dim, Number>;
  using CellIntegratorP = CellIntegrator<dim, 1, Number>;

  void
  initialize(dealii::MatrixFree<dim, Number> const & matrix_free_in,
             PMLOperatorData<dim> const              data_in)
  {
    AssertThrow(data_in.pml_damping, dealii::ExcMessage("No PML damping function provided"));
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;
  }

  void
  evaluate(BlockVectorType & dst, BlockVectorType const & src) const
  {
    dst.zero_out_ghost_values();
    matrix_free->cell_loop(&This::cell_loop, this, dst, src, true);
  }

  void
  evaluate_add(BlockVectorType & dst, BlockVectorType const & src) const
  {
    dst.zero_out_ghost_values();
    matrix_free->cell_loop(&This::cell_loop, this, dst, src, false);
  }

private:
  static inline DEAL_II_ALWAYS_INLINE //
    vector
    multiply_component_wise(vector const & vec1, vector const & vec2)
  {
    vector result;
    for(unsigned int d = 0; d < dim; ++d)
      result[d] = vec1[d] * vec2[d];
    return result;
  }

  static inline DEAL_II_ALWAYS_INLINE //
    vector
    diago(tensor const & tensor_in)
  {
    vector result;
    for(unsigned int d = 0; d < dim; ++d)
      result[d] = tensor_in[d][d];
    return result;
  }

  void
  cell_loop(dealii::MatrixFree<dim, Number> const & matrix_free_in,
            BlockVectorType &                       dst,
            BlockVectorType const &                 src,
            Range const &                           cell_range) const
  {
    CellIntegratorU auxiliary(matrix_free_in, data.dof_index_velocity, data.quad_index);
    CellIntegratorU velocity(matrix_free_in, data.dof_index_velocity, data.quad_index);
    CellIntegratorP pressure(matrix_free_in, data.dof_index_pressure, data.quad_index);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if(matrix_free_in.get_cell_category(cell) != numbers::pml_material_id)
        continue;

      pressure.reinit(cell);

      velocity.reinit(cell);
      velocity.gather_evaluate(src.block(data.block_index_velocity),
                               dealii::EvaluationFlags::values |
                                 dealii::EvaluationFlags::gradients);

      auxiliary.reinit(cell);
      auxiliary.gather_evaluate(src.block(data.block_index_auxiliary),
                                dealii::EvaluationFlags::values);

      for(unsigned int q = 0; q < auxiliary.n_q_points; ++q)
      {
        // sigma_d is the diagonal of the damping tensor. since we restrict ourselfs to
        // pmls in normal direction of the coordinate system the off-diagonal entries are
        // always 0.
        // TODO: can we somehow use a vector that only holds required PML values? currently we
        // use the same amout of storage for PML auxiliary and velocity. (will not make anything
        // faster but still...)
        // TODO: optimize performance by precomputing these values
        // TODO: maybe Patrick knows right away how we can generalize the PML. We have to
        // look at the derivation of the PML formulation in detail.
        vector sigma_d =
          FunctionEvaluator<1, dim, Number>::value(*data.pml_damping,
                                                   auxiliary.quadrature_point(q),
                                                   std::numeric_limits<double>::max());

        // (q, -sigma_d * auxiliary)
        pressure.submit_value(-sigma_d * auxiliary.get_value(q), q);

        // (w[i], sigma_d[i] * u[i])
        velocity.submit_value(multiply_component_wise(sigma_d, velocity.get_value(q)), q);

        // (b, sigma * a - D * u)
        auxiliary.submit_value(multiply_component_wise(sigma_d, auxiliary.get_value(q)) -
                                 diago(velocity.get_gradient(q)),
                               q);
      }

      pressure.integrate_scatter(dealii::EvaluationFlags::values,
                                 dst.block(data.block_index_pressure));
      velocity.integrate_scatter(dealii::EvaluationFlags::values,
                                 dst.block(data.block_index_velocity));
      auxiliary.integrate_scatter(dealii::EvaluationFlags::values,
                                  dst.block(data.block_index_auxiliary));
    }
  }

  dealii::MatrixFree<dim, Number> const * matrix_free = nullptr;

  PMLOperatorData<dim> data;

  Operators::PMLKernel<dim, Number> kernel;
};

} // namespace ExaDG::Acoustics
