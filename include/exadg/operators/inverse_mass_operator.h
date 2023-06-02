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

#ifndef INCLUDE_OPERATORS_INVERSEMASSMATRIX_H_
#define INCLUDE_OPERATORS_INVERSEMASSMATRIX_H_

// deal.II
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/operators.h>

// ExaDG
#include <exadg/matrix_free/integrators.h>
#include <exadg/operators/mass_operator.h>
#include <exadg/solvers_and_preconditioners/preconditioners/block_jacobi_preconditioner.h>

namespace ExaDG
{
struct InverseMassOperatorData
{
  unsigned int dof_index               = 0;
  unsigned int quad_index              = 0;
  unsigned int only_eval_cell_category = dealii::numbers::invalid_unsigned_int;

  // only relevant if an explicit matrix-free inverse mass operator is not available
  bool implement_block_diagonal_preconditioner_matrix_free = true;

  // only relevant if elementwise mass operators are inverted by elementwise
  // iterative solvers with matrix-free implementation
  SolverData solver_data_block_diagonal = SolverData(1000, 1e-12, 1e-20);
};

template<int dim, int n_components, typename Number>
class InverseMassOperator
{
private:
  typedef dealii::LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef InverseMassOperator<dim, n_components, Number> This;

  typedef CellIntegrator<dim, n_components, Number> Integrator;

  // use a template parameter of -1 to select the precompiled version of this operator
  typedef dealii::MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, n_components, Number>
    ExplicitMatrixFreeInverseMass;

  typedef std::pair<unsigned int, unsigned int> Range;

public:
  InverseMassOperator()
    : matrix_free(nullptr),
      dof_index(0),
      quad_index(0),
      only_eval_cell_category(dealii::numbers::invalid_unsigned_int),
      explicit_matrix_free_inverse_mass_available(false)
  {
  }

  void
  initialize(dealii::MatrixFree<dim, Number> const & matrix_free_in,
             InverseMassOperatorData const           inverse_mass_operator_data)
  {
    this->matrix_free       = &matrix_free_in;
    dof_index               = inverse_mass_operator_data.dof_index;
    quad_index              = inverse_mass_operator_data.quad_index;
    only_eval_cell_category = inverse_mass_operator_data.only_eval_cell_category;
    if(only_eval_cell_category != dealii::numbers::invalid_unsigned_int)
      AssertThrow(
        explicit_matrix_free_inverse_mass_available,
        dealii::ExcMessage(
          "Only evaluation certain cell category currently only supported if explicit_matrix_free_inverse_mass_available."));

    dealii::FiniteElement<dim> const & fe = matrix_free->get_dof_handler(dof_index).get_fe();

    if(fe.conforms(dealii::FiniteElementData<dim>::L2))
    {
      // this checks if we have a tensor-product element
      if(fe.base_element(0).dofs_per_cell == dealii::Utilities::pow(fe.degree + 1, dim))
        explicit_matrix_free_inverse_mass_available = true;
    }
    else
      AssertThrow(false, dealii::ExcMessage("InverseMassOperator only implemented for DG!"));

    if(not(explicit_matrix_free_inverse_mass_available))
    {
      // initialize mass operator
      dealii::AffineConstraints<Number> constraint;
      constraint.clear();
      constraint.close();

      MassOperatorData<dim> mass_operator_data;
      mass_operator_data.dof_index  = dof_index;
      mass_operator_data.quad_index = quad_index;
      mass_operator_data.implement_block_diagonal_preconditioner_matrix_free =
        inverse_mass_operator_data.implement_block_diagonal_preconditioner_matrix_free;
      mass_operator_data.solver_block_diagonal         = Elementwise::Solver::CG;
      mass_operator_data.preconditioner_block_diagonal = Elementwise::Preconditioner::None;
      mass_operator_data.solver_data_block_diagonal =
        inverse_mass_operator_data.solver_data_block_diagonal;

      mass_operator.initialize(*matrix_free, constraint, mass_operator_data);

      elementwise_inverse_mass =
        std::make_shared<BlockJacobiPreconditioner<MassOperator<dim, n_components, Number>>>(
          mass_operator);
    }
  }

  void
  apply(VectorType & dst, VectorType const & src) const
  {
    dst.zero_out_ghost_values();

    if(explicit_matrix_free_inverse_mass_available)
    {
      matrix_free->cell_loop(&This::cell_loop, this, dst, src);
    }
    else
    {
      elementwise_inverse_mass->vmult(dst, src);
    }
  }

private:
  void
  cell_loop(dealii::MatrixFree<dim, Number> const &,
            VectorType &       dst,
            VectorType const & src,
            Range const &      cell_range) const
  {
    Integrator                    integrator(*matrix_free, dof_index, quad_index);
    ExplicitMatrixFreeInverseMass inverse(integrator);

    bool const check_cell_category =
      (only_eval_cell_category != dealii::numbers::invalid_unsigned_int);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if(check_cell_category)
        if(matrix_free->get_cell_category(cell) != only_eval_cell_category)
          continue;

      integrator.reinit(cell);
      integrator.read_dof_values(src, 0);

      inverse.apply(integrator.begin_dof_values(), integrator.begin_dof_values());

      integrator.set_dof_values(dst, 0);
    }
  }

  dealii::MatrixFree<dim, Number> const * matrix_free;

  unsigned int dof_index, quad_index, only_eval_cell_category;

  bool explicit_matrix_free_inverse_mass_available;

  MassOperator<dim, n_components, Number> mass_operator;

  std::shared_ptr<BlockJacobiPreconditioner<MassOperator<dim, n_components, Number>>>
    elementwise_inverse_mass;
};

} // namespace ExaDG


#endif /* INCLUDE_OPERATORS_INVERSEMASSMATRIX_H_ */
