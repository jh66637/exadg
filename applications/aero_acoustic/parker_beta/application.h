/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2023 by the ExaDG authors
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

#ifndef APPLICATIONS_AERO_ACOUSTIC_TEST_CASES_PARKER_H_
#define APPLICATIONS_AERO_ACOUSTIC_TEST_CASES_PARKER_H_

#include <exadg/functions_and_boundary_conditions/linear_interpolation.h>

#include "grid.h"

namespace ExaDG
{
double const rampup_duration = 0.1;
double const start_acoustic  = 0.1;
double const end_time        = start_acoustic + rampup_duration + 0.5;
double const dt_max          = 1e-5;

bool const         CONSIDER_BACKCOUPLING = false;
unsigned int const REFINEMENTS_BND_LAYER = 1;

namespace AcousticsAeroAcoustic
{
using namespace Acoustics;

template<int dim, typename Number>
class Application : public ApplicationBase<dim, Number>
{
public:
  Application(std::string input_file_in, MPI_Comm const & comm)
    : ApplicationBase<dim, Number>(input_file_in, comm)
  {
  }

  void
  add_parameters(dealii::ParameterHandler & prm) final
  {
    ApplicationBase<dim, Number>::add_parameters(prm);
    prm.enter_subsection("Application");
    prm.add_parameter("BulkVelocity", bulk_velocity);
    prm.add_parameter("SpeedOfSound", this->param.speed_of_sound, "Speed of sound.");
    prm.add_parameter("CFLAcoustics", this->param.cfl, "Courant Number.");
    prm.leave_subsection();
  }

private:
  void
  set_parameters() final
  {
    this->param.formulation               = Formulation::SkewSymmetric;
    this->param.aero_acoustic_source_term = true;

    this->param.calculation_of_time_step_size = TimeStepCalculation::CFL;
    this->param.order_time_integrator         = 2;
    this->param.start_with_low_order          = true;
    this->param.adaptive_time_stepping        = true;

    this->param.start_time = start_acoustic;
    this->param.end_time   = end_time;

    this->param.solver_info_data.interval_time = (this->param.end_time - this->param.start_time);

    this->param.grid.triangulation_type = TriangulationType::Distributed;

    this->param.mapping_degree = 2;
    this->param.degree_p       = this->param.degree_u;
    this->param.degree_u       = this->param.degree_p;
  }

  void
  create_grid(Grid<dim> & grid, std::shared_ptr<dealii::Mapping<dim>> & mapping) final
  {
    auto const lambda_create_triangulation =
      [&](dealii::Triangulation<dim, dim> &                        tria,
          std::vector<dealii::GridTools::PeriodicFacePair<
            typename dealii::Triangulation<dim>::cell_iterator>> & periodic_face_pairs,
          unsigned int const                                       global_refinements,
          std::vector<unsigned int> const &                        vector_local_refinements)
    {
      (void)periodic_face_pairs;
      (void)vector_local_refinements;

      Parker::create_triangulation_fluid(tria, Parker::channel_width, global_refinements, 0);
    };

    GridUtilities::create_triangulation<dim>(
      grid, this->mpi_comm, this->param.grid, lambda_create_triangulation, {});

    GridUtilities::create_mapping(mapping,
                                  this->param.grid.element_type,
                                  this->param.mapping_degree);
  }

  void
  set_boundary_descriptor() final
  {
    auto Y = std::make_shared<dealii::Functions::ConstantFunction<dim>>(0.01);
    this->boundary_descriptor->admittance_bc.insert(std::make_pair(0, Y));
    this->boundary_descriptor->admittance_bc.insert(std::make_pair(3, Y));

    auto ABC = std::make_shared<dealii::Functions::ConstantFunction<dim>>(1.0);
    this->boundary_descriptor->admittance_bc.insert(std::make_pair(1, ABC));
    this->boundary_descriptor->admittance_bc.insert(std::make_pair(2, ABC));
  }

  void
  set_field_functions() final
  {
    this->field_functions->initial_solution_velocity.reset(
      new dealii::Functions::ZeroFunction<dim>(dim));
    this->field_functions->initial_solution_pressure.reset(
      new dealii::Functions::ZeroFunction<dim>(1));

    this->field_functions->right_hand_side.reset(new dealii::Functions::ZeroFunction<dim>(1));
  }

  std::shared_ptr<PostProcessorBase<dim, Number>>
  create_postprocessor() final
  {
    PostProcessorData<dim> pp_data;

    // write output for visualization of results
    pp_data.output_data.time_control_data.is_active        = this->output_parameters.write;
    pp_data.output_data.time_control_data.start_time       = 0.0;
    pp_data.output_data.time_control_data.trigger_interval = this->param.end_time / 20.0;
    pp_data.output_data.directory          = this->output_parameters.directory + "vtu/";
    pp_data.output_data.filename           = this->output_parameters.filename + "_acoustic";
    pp_data.output_data.write_velocity     = true;
    pp_data.output_data.write_pressure     = true;
    pp_data.output_data.write_processor_id = false;
    pp_data.output_data.write_boundary_IDs = true;
    pp_data.output_data.write_higher_order = true;
    pp_data.output_data.degree             = this->param.degree_u;


    // record pressure and velocity in the wake of the obstacle
    pp_data.pointwise_output_data.write_pressure                     = true;
    pp_data.pointwise_output_data.write_velocity                     = false;
    pp_data.pointwise_output_data.time_control_data.is_active        = true;
    pp_data.pointwise_output_data.time_control_data.start_time       = start_acoustic;
    pp_data.pointwise_output_data.time_control_data.end_time         = end_time;
    pp_data.pointwise_output_data.time_control_data.trigger_interval = dt_max;

    pp_data.pointwise_output_data.directory =
      this->output_parameters.directory + "pointwise_output/";
    pp_data.pointwise_output_data.filename = this->output_parameters.filename + "_acoustic";
    pp_data.pointwise_output_data.update_points_before_evaluation = false;

    pp_data.pointwise_output_data.evaluation_points.push_back(
      {0.5 * Parker::plate_length, 0.49 * Parker::channel_height, 0.5 * Parker::channel_width});

    std::shared_ptr<PostProcessorBase<dim, Number>> pp;
    pp.reset(new PostProcessor<dim, Number>(pp_data, this->mpi_comm));

    return pp;
  }

  double bulk_velocity = 32.0;
};
} // namespace AcousticsAeroAcoustic

namespace FluidAeroAcoustic
{
using namespace ExaDG;
using namespace IncNS;

template<int dim>
class UniformInflow : public dealii::Function<dim>
{
public:
  explicit UniformInflow(double const bulk_velocity_in)
    : dealii::Function<dim>(dim, 0.0), bulk_velocity(bulk_velocity_in)
  {
  }

  double
  value(dealii::Point<dim> const &, unsigned int const component = 0) const final
  {
    // flow in x-direction
    if(component == 0)
      return bulk_velocity;
    else
      return 0.0;
  }

private:
  double const bulk_velocity;
};

template<int dim, typename Number>
class Application : public ApplicationBase<dim, Number>
{
public:
  Application(std::string input_file_in, MPI_Comm const & comm)
    : ApplicationBase<dim, Number>(input_file_in, comm)
  {
  }

  void
  add_parameters(dealii::ParameterHandler & prm) final
  {
    ApplicationBase<dim, Number>::add_parameters(prm);

    prm.enter_subsection("Application");
    prm.add_parameter("BulkVelocity", bulk_velocity);
    prm.add_parameter("CFLFluid", this->param.cfl, "Courant Number.");
    prm.add_parameter("TemporalDiscretizationFluid",
                      this->param.temporal_discretization,
                      "Temporal discretization of the fluid.");
    prm.leave_subsection();
  }

private:
  void
  set_parameters() final
  {
    // MATHEMATICAL MODEL
    this->param.problem_type                = ProblemType::Unsteady;
    this->param.equation_type               = EquationType::NavierStokes;
    this->param.formulation_viscous_term    = FormulationViscousTerm::LaplaceFormulation;
    this->param.formulation_convective_term = FormulationConvectiveTerm::DivergenceFormulation;
    this->param.right_hand_side             = CONSIDER_BACKCOUPLING;



    // PHYSICAL QUANTITIES
    this->param.start_time = 0.0;
    this->param.end_time   = end_time;
    this->param.viscosity  = 1.58e-5;

    // TEMPORAL DISCRETIZATION
    this->param.solver_type                     = SolverType::Unsteady;
    this->param.treatment_of_convective_term    = TreatmentOfConvectiveTerm::Explicit;
    this->param.calculation_of_time_step_size   = TimeStepCalculation::CFL;
    this->param.order_time_integrator           = 2;
    this->param.adaptive_time_stepping          = true;
    this->param.start_with_low_order            = true;
    this->param.cfl_exponent_fe_degree_velocity = 1.5;

    this->param.max_velocity = bulk_velocity;

    // output of solver information
    this->param.solver_info_data.interval_time =
      (this->param.end_time - this->param.start_time) / 100;

    // SPATIAL DISCRETIZATION
    this->param.grid.triangulation_type           = TriangulationType::Distributed;
    this->param.grid.create_coarse_triangulations = REFINEMENTS_BND_LAYER > 0;

    this->param.mapping_degree = 2;
    this->param.degree_p       = DegreePressure::MixedOrder;

    // convective term
    if(this->param.formulation_convective_term == FormulationConvectiveTerm::DivergenceFormulation)
      this->param.upwind_factor = 0.5; // allows using larger CFL values for explicit formulations

    // divergence penalty
    this->param.use_divergence_penalty                     = true;
    this->param.divergence_penalty_factor                  = 1.0e0;
    this->param.use_continuity_penalty                     = true;
    this->param.continuity_penalty_factor                  = this->param.divergence_penalty_factor;
    this->param.continuity_penalty_components              = ContinuityPenaltyComponents::Normal;
    this->param.continuity_penalty_use_boundary_data       = true;
    this->param.apply_penalty_terms_in_postprocessing_step = true;

    // viscous term
    this->param.IP_formulation_viscous = InteriorPenaltyFormulation::SIPG;

    // NUMERICAL PARAMETERS
    this->param.implement_block_diagonal_preconditioner_matrix_free = false;
    this->param.use_cell_based_face_loops                           = false;
    this->param.quad_rule_linearization = QuadratureRuleLinearization::Overintegration32k;

    // PROJECTION METHODS

    // pressure Poisson equation
    this->param.solver_pressure_poisson              = SolverPressurePoisson::CG;
    this->param.solver_data_pressure_poisson         = SolverData(1000, ABS_TOL, REL_TOL, 30);
    this->param.preconditioner_pressure_poisson      = PreconditionerPressurePoisson::Multigrid;
    this->param.multigrid_data_pressure_poisson.type = MultigridType::cphMG;

    this->param.multigrid_data_pressure_poisson.smoother_data.smoother =
      MultigridSmoother::Chebyshev;
    this->param.multigrid_data_pressure_poisson.p_sequence = PSequenceType::Bisect;

    this->param.multigrid_data_pressure_poisson.smoother_data.iterations = 5;
    this->param.multigrid_data_pressure_poisson.coarse_problem.solver =
      MultigridCoarseGridSolver::CG;
    this->param.multigrid_data_pressure_poisson.coarse_problem.preconditioner =
      MultigridCoarseGridPreconditioner::AMG;
    this->param.update_preconditioner_pressure_poisson = false;

    // projection step
    this->param.solver_projection                = SolverProjection::CG;
    this->param.solver_data_projection           = SolverData(1000, ABS_TOL, REL_TOL);
    this->param.preconditioner_projection        = PreconditionerProjection::InverseMassMatrix;
    this->param.update_preconditioner_projection = false;

    // HIGH-ORDER DUAL SPLITTING SCHEME

    // formulations
    this->param.order_extrapolation_pressure_nbc =
      this->param.order_time_integrator <= 2 ? this->param.order_time_integrator : 2;
    this->param.use_outflow_bc_convective_term = true;

    // PRESSURE-CORRECTION SCHEME

    // formulation
    this->param.order_pressure_extrapolation = 1;
    this->param.rotational_formulation       = true;

    // momentum step
    if(this->param.temporal_discretization == TemporalDiscretization::BDFPressureCorrection)
    {
      // Newton solver
      this->param.newton_solver_data_momentum = Newton::SolverData(100, ABS_TOL, REL_TOL);

      // linear solver
      this->param.solver_momentum = SolverMomentum::GMRES;
      if(this->param.treatment_of_convective_term == TreatmentOfConvectiveTerm::Implicit)
        this->param.solver_data_momentum = SolverData(1e4, ABS_TOL_LINEAR, REL_TOL_LINEAR, 100);
      else
        this->param.solver_data_momentum = SolverData(1e4, ABS_TOL, REL_TOL, 100);

      this->param.preconditioner_momentum = MomentumPreconditioner::InverseMassMatrix;
    }

    // COUPLED NAVIER-STOKES SOLVER
    this->param.use_scaling_continuity = false;

    // nonlinear solver (Newton solver)
    this->param.newton_solver_data_coupled = Newton::SolverData(100, ABS_TOL, REL_TOL);

    // linear solver
    this->param.solver_coupled = SolverCoupled::GMRES;
    if(this->param.treatment_of_convective_term == TreatmentOfConvectiveTerm::Implicit)
      this->param.solver_data_coupled = SolverData(1e3, ABS_TOL_LINEAR, REL_TOL_LINEAR, 100);
    else
      this->param.solver_data_coupled = SolverData(1e3, ABS_TOL, REL_TOL, 100);

    // preconditioning linear solver
    this->param.preconditioner_coupled = PreconditionerCoupled::BlockTriangular;

    // preconditioner velocity/momentum block
    this->param.preconditioner_velocity_block = MomentumPreconditioner::InverseMassMatrix;

    // preconditioner Schur-complement block
    this->param.preconditioner_pressure_block      = SchurComplementPreconditioner::CahouetChabard;
    this->param.multigrid_data_pressure_block.type = MultigridType::cphMG;
    this->param.multigrid_data_pressure_block.coarse_problem.solver = MultigridCoarseGridSolver::CG;
    this->param.multigrid_data_pressure_block.coarse_problem.preconditioner =
      MultigridCoarseGridPreconditioner::AMG;
  }



  void
  create_grid(Grid<dim> &                                       grid,
              std::shared_ptr<dealii::Mapping<dim>> &           mapping,
              std::shared_ptr<MultigridMappings<dim, Number>> & multigrid_mappings) final
  {
    auto const lambda_create_triangulation =
      [&](dealii::Triangulation<dim, dim> &                        tria,
          std::vector<dealii::GridTools::PeriodicFacePair<
            typename dealii::Triangulation<dim>::cell_iterator>> & periodic_face_pairs,
          unsigned int const                                       global_refinements,
          std::vector<unsigned int> const &                        vector_local_refinements)
    {
      (void)periodic_face_pairs;
      (void)vector_local_refinements;

      Parker::create_triangulation_fluid(tria,
                                         Parker::channel_width,
                                         global_refinements,
                                         REFINEMENTS_BND_LAYER);
    };

    GridUtilities::create_triangulation_with_multigrid<dim>(grid,
                                                            this->mpi_comm,
                                                            this->param.grid,
                                                            this->param.involves_h_multigrid(),
                                                            lambda_create_triangulation,
                                                            {});

    // mappings
    GridUtilities::create_mapping_with_multigrid(mapping,
                                                 multigrid_mappings,
                                                 this->param.grid.element_type,
                                                 this->param.mapping_degree,
                                                 this->param.mapping_degree_coarse_grids,
                                                 this->param.involves_h_multigrid());
  }

  void
  set_boundary_descriptor() final
  {
    // obstacle
    this->boundary_descriptor->pressure->neumann_bc.insert(0);
    this->boundary_descriptor->velocity->dirichlet_bc.insert(
      std::make_pair(0, std::make_shared<dealii::Functions::ZeroFunction<dim>>(dim)));

    // slip
    this->boundary_descriptor->velocity->symmetry_bc.insert(
      std::make_pair(3, new dealii::Functions::ZeroFunction<dim>(dim)));
    this->boundary_descriptor->pressure->neumann_bc.insert(3);

    // inflow
    this->boundary_descriptor->pressure->neumann_bc.insert(1);
    this->boundary_descriptor->velocity->dirichlet_bc.insert(
      std::make_pair(1, std::make_shared<UniformInflow<dim>>(bulk_velocity)));

    // outflow
    this->boundary_descriptor->pressure->dirichlet_bc.insert(
      std::make_pair(2, std::make_shared<dealii::Functions::ZeroFunction<dim>>(1)));
    this->boundary_descriptor->velocity->neumann_bc.insert(
      std::make_pair(2, std::make_shared<dealii::Functions::ZeroFunction<dim>>(dim)));
  }

  void
  set_field_functions() final
  {
    this->field_functions->initial_solution_velocity.reset(
      new dealii::Functions::ZeroFunction<dim>(dim));
    this->field_functions->initial_solution_pressure.reset(
      new dealii::Functions::ZeroFunction<dim>(1));
    this->field_functions->analytical_solution_pressure.reset(
      new dealii::Functions::ZeroFunction<dim>(1));
    this->field_functions->right_hand_side.reset(new dealii::Functions::ZeroFunction<dim>(dim));
  }

  std::shared_ptr<IncNS::PostProcessorBase<dim, Number>>
  create_postprocessor() final
  {
    PostProcessorData<dim> pp_data;

    // write output for visualization of results
    pp_data.output_data.time_control_data.is_active        = this->output_parameters.write;
    pp_data.output_data.time_control_data.start_time       = 0.0;
    pp_data.output_data.time_control_data.trigger_interval = this->param.end_time / 20.0;
    pp_data.output_data.directory          = this->output_parameters.directory + "vtu/";
    pp_data.output_data.filename           = this->output_parameters.filename + "_fluid";
    pp_data.output_data.write_processor_id = false;
    pp_data.output_data.write_higher_order = true;
    pp_data.output_data.degree             = this->param.degree_u;
    pp_data.output_data.write_boundary_IDs = true;


    std::shared_ptr<PostProcessorBase<dim, Number>> pp;
    pp.reset(new PostProcessor<dim, Number>(pp_data, this->mpi_comm));

    return pp;
  }

  double bulk_velocity = 32.0;

  // solver tolerances
  double const ABS_TOL = 1.e-12;
  double const REL_TOL = 1.e-6;

  double const ABS_TOL_LINEAR = 1.e-12;
  double const REL_TOL_LINEAR = 1.e-2;
};

} // namespace FluidAeroAcoustic


namespace AeroAcoustic
{
template<int dim>
class BlendInFunction : public Utilities::SpatialAwareFunction<dim>
{
public:
  BlendInFunction(double const blend_in_start, double const blend_in_end)
    : Utilities::SpatialAwareFunction<dim>(1, 0.0), start(blend_in_start), end(blend_in_end)
  {
  }

  double
  value(dealii::Point<dim> const &, unsigned int const) const final
  {
    AssertThrow(false, dealii::ExcMessage("Should not end up here"));
    return {};
  }

  bool
  varies_in_space(double const) const final
  {
    return false;
  }

  double
  compute_time_factor(double const time) const final
  {
    double const pi = dealii::numbers::PI;
    double const T  = end - start;
    double const t  = time - start;

    return 0.5 * (1.0 - std::cos(pi * std::min(t / T, 1.0)));
  }

private:
  double const start;
  double const end;
};


template<int dim, typename Number>
class Application : public ApplicationBase<dim, Number>
{
public:
  Application(std::string input_file, MPI_Comm const & comm)
    : ApplicationBase<dim, Number>(input_file, comm)
  {
  }

private:
  void
  set_single_field_solvers(std::string input_file, MPI_Comm const & comm) final
  {
    this->acoustic =
      std::make_shared<AcousticsAeroAcoustic::Application<dim, Number>>(input_file, comm);
    this->fluid = std::make_shared<FluidAeroAcoustic::Application<dim, Number>>(input_file, comm);
  }


  void
  set_field_functions() final
  {
    this->field_functions->source_term_blend_in =
      std::make_shared<BlendInFunction<dim>>(start_acoustic, rampup_duration);
  }
};
} // namespace AeroAcoustic

} // namespace ExaDG
#include <exadg/aero_acoustic/user_interface/implement_get_application.h>

#endif /*APPLICATIONS_AERO_ACOUSTIC_TEST_CASES_PARKER_H_*/
