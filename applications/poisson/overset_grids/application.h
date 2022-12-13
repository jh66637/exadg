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
#ifndef APPLICATIONS_POISSON_TEST_CASES_TEMPLATE_H_
#define APPLICATIONS_POISSON_TEST_CASES_TEMPLATE_H_

namespace ExaDG
{
namespace Poisson
{
template<int dim, int n_components, typename Number>
class Application : public ApplicationBase<dim, n_components, Number>
{
private:
  static unsigned int const rank =
    (n_components == 1) ? 0 : ((n_components == dim) ? 1 : dealii::numbers::invalid_unsigned_int);

public:
  Application(std::string input_file, MPI_Comm const & comm)
    : ApplicationBase<dim, n_components, Number>(input_file, comm)
  {
  }

  void
  add_parameters(dealii::ParameterHandler & prm)
  {
    ApplicationBase<dim, n_components, Number>::add_parameters(prm);

    // clang-format off
    prm.enter_subsection("Application");
      prm.add_parameter("AdditionalRefinementDomain2", additional_refinement_domain_2, "Additional refinement Domain 2.", dealii::Patterns::Integer(0));
    prm.leave_subsection();
    // clang-format on
  }


  void
  set_parameters() final
  {
    Parameters & p = this->param;

    // MATHEMATICAL MODEL
    p.right_hand_side = true;

    // SPATIAL DISCRETIZATION
    p.grid.triangulation_type = TriangulationType::Distributed;
    p.grid.mapping_degree     = 3;
    p.spatial_discretization  = SpatialDiscretization::DG;
    p.IP_factor               = 1.0e0;

    // SOLVER
    p.solver                      = Solver::CG;
    p.solver_data.abs_tol         = 1.e-20;
    p.solver_data.rel_tol         = 1.e-10;
    p.solver_data.max_iter        = 1e4;
    p.compute_performance_metrics = true;
    p.preconditioner              = Preconditioner::Multigrid;
    p.multigrid_data.type         = MultigridType::cphMG;
    p.multigrid_data.p_sequence   = PSequenceType::Bisect;
    // MG smoother
    p.multigrid_data.smoother_data.smoother        = MultigridSmoother::Chebyshev;
    p.multigrid_data.smoother_data.iterations      = 5;
    p.multigrid_data.smoother_data.smoothing_range = 20;
    // MG coarse grid solver
    p.multigrid_data.coarse_problem.solver              = MultigridCoarseGridSolver::CG;
    p.multigrid_data.coarse_problem.preconditioner      = MultigridCoarseGridPreconditioner::AMG;
    p.multigrid_data.coarse_problem.solver_data.rel_tol = 1.e-3;
  }

  void
  create_grid() final
  {
    dealii::Triangulation<dim> domain1;
    {
      double const       right = 1.0;
      dealii::Point<dim> p1, p2;
      p1[0] = 0.0;
      p1[1] = 0.5;
      p2[0] = right;
      p2[1] = 1.5;
      dealii::GridGenerator::subdivided_hyper_rectangle(domain1, {3, 3}, p1, p2);
    }

    dealii::Triangulation<dim> domain2;
    {
      dealii::Triangulation<dim> domain2_temp;

      double const       left = 0.5;
      dealii::Point<dim> p1, p2;
      p1[0] = left;
      p1[1] = 0.0;
      p2[0] = left + 1.0;
      p2[1] = 1.0;
      dealii::GridGenerator::subdivided_hyper_rectangle(domain2_temp, {2, 2}, p1, p2);
      domain2_temp.refine_global(additional_refinement_domain_2);
      dealii::GridGenerator::flatten_triangulation(domain2_temp, domain2);
    }

    dealii::GridGenerator::merge_triangulations(
      domain1, domain2, *this->grid->triangulation, 0.0, true, true);

    this->grid->triangulation->refine_global(this->param.grid.n_refine_global);
  }

  void
  set_boundary_descriptor() final
  {
    this->boundary_descriptor->dirichlet_bc.insert(
      std::make_pair(0, new dealii::Functions::ZeroFunction<dim>(dim)));

    dealii::types::boundary_id const bnd_id_max =
      std::numeric_limits<dealii::types::boundary_id>::max();

    this->boundary_descriptor->overset_bc.insert(std::make_pair(bnd_id_max - 1, bnd_id_max - 2));
    this->boundary_descriptor->overset_bc.insert(std::make_pair(bnd_id_max - 2, bnd_id_max - 1));
  }

  void
  set_field_functions() final
  {
    // these lines show exemplarily how the field functions are filled
    this->field_functions->initial_solution.reset(
      new dealii::Functions::ZeroFunction<dim>(n_components));
    this->field_functions->right_hand_side.reset(
      new dealii::Functions::ConstantFunction<dim>(1.0, n_components));
  }

  std::shared_ptr<PostProcessorBase<dim, Number>>
  create_postprocessor() final
  {
    PostProcessorData<dim> pp_data;
    pp_data.output_data.time_control_data.is_active = this->output_parameters.write;
    pp_data.output_data.directory                   = this->output_parameters.directory + "vtu/";
    pp_data.output_data.filename                    = this->output_parameters.filename;
    pp_data.output_data.write_higher_order          = true;
    pp_data.output_data.degree                      = this->param.degree;

    std::shared_ptr<PostProcessorBase<dim, Number>> pp;
    pp.reset(new PostProcessor<dim, Number>(pp_data, this->mpi_comm));

    return pp;
  }

  int additional_refinement_domain_2 = 2;
};

} // namespace Poisson

} // namespace ExaDG

#include <exadg/poisson/user_interface/implement_get_application.h>

#endif /* APPLICATIONS_POISSON_TEST_CASES_TEMPLATE_H_ */
