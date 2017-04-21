/*
 * DGNavierStokesDualSplitting.h
 *
 *  Created on: Jun 27, 2016
 *      Author: fehn
 */

#ifndef INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_DG_NAVIER_STOKES_DUAL_SPLITTING_H_
#define INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_DG_NAVIER_STOKES_DUAL_SPLITTING_H_

#include "../../incompressible_navier_stokes/preconditioners/multigrid_preconditioner_navier_stokes.h"
#include "../../incompressible_navier_stokes/spatial_discretization/dg_navier_stokes_projection_methods.h"
#include "../../incompressible_navier_stokes/spatial_discretization/helmholtz_operator.h"
#include "../../incompressible_navier_stokes/spatial_discretization/pressure_neumann_bc_convective_term.h"
#include "../../incompressible_navier_stokes/spatial_discretization/pressure_neumann_bc_viscous_term.h"
#include "../../incompressible_navier_stokes/spatial_discretization/velocity_divergence_convective_term.h"
#include "../include/solvers_and_preconditioners/jacobi_preconditioner.h"
#include "solvers_and_preconditioners/inverse_mass_matrix_preconditioner.h"
#include "solvers_and_preconditioners/iterative_solvers.h"
#include "solvers_and_preconditioners/newton_solver.h"

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
class DGNavierStokesDualSplitting : public DGNavierStokesProjectionMethods<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>
{
public:
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::value_type value_type;

  DGNavierStokesDualSplitting(parallel::distributed::Triangulation<dim> const &triangulation,
                              InputParametersNavierStokes<dim> const          &parameter)
    :
    DGNavierStokesProjectionMethods<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>(triangulation,parameter),
    fe_param(parameter),
    sum_alphai_ui(nullptr),
    evaluation_time(0.0),
    scaling_factor_time_derivative_term(1.0)
  {}

  virtual ~DGNavierStokesDualSplitting()
  {}

  void setup_solvers(double const &time_step_size,
                     double const &scaling_factor_time_derivative_term);

  /*
   *  implicit solution of convective step
   */
  void solve_nonlinear_convective_problem (parallel::distributed::Vector<value_type>       &dst,
                                           parallel::distributed::Vector<value_type> const &sum_alphai_ui,
                                           double const                                    &eval_time,
                                           double const                                    &scaling_factor_mass_matrix_term,
                                           unsigned int                                    &newton_iterations,
                                           unsigned int                                    &linear_iterations);

  /*
   *  The implementation of the Newton solver requires that the underlying operator
   *  implements a function called "evaluate_nonlinear_residual"
   */
  void evaluate_nonlinear_residual (parallel::distributed::Vector<value_type>       &dst,
                                    parallel::distributed::Vector<value_type> const &src);

  /*
   * The implementation of the Newton solver requires that the underlying operator
   * implements a function called "initialize_vector_for_newton_solver"
   */
  void initialize_vector_for_newton_solver(parallel::distributed::Vector<value_type> &src) const
  {
    this->initialize_vector_velocity(src);
  }

  /*
   * The implementation of the Newton solver requires that the underlying operator
   * implements a function called "set_solution_linearization"
   */
  void set_solution_linearization(parallel::distributed::Vector<value_type> const &solution_linearization)
  {
    velocity_linearization = solution_linearization;
  }

  /*
   *  To solve the linearized convective problem, the underlying operator
   *  has to implement a function called "vmult"
   */
  void vmult (parallel::distributed::Vector<value_type>       &dst,
              parallel::distributed::Vector<value_type> const &src) const;

  /*
   *  This function calculates the matrix vector product for the linearized convective problem.
   */
  void apply_linearized_convective_problem (parallel::distributed::Vector<value_type>       &dst,
                                            parallel::distributed::Vector<value_type> const &src) const;

  // convective term
  void evaluate_convective_term_and_apply_inverse_mass_matrix(parallel::distributed::Vector<value_type>       &dst,
                                                              parallel::distributed::Vector<value_type> const &src,
                                                              value_type const                                evaluation_time) const;

  // body forces
  void  evaluate_body_force_and_apply_inverse_mass_matrix (parallel::distributed::Vector<value_type>  &dst,
                                                           const value_type                           evaluation_time) const;

  // rhs pressure: velocity divergence
  void apply_velocity_divergence_term(parallel::distributed::Vector<value_type>        &dst,
                                      const parallel::distributed::Vector<value_type>  &src) const;

  void rhs_velocity_divergence_term(parallel::distributed::Vector<value_type>        &dst,
                                    double const                                     &evaluation_time) const;

  void rhs_ppe_div_term_body_forces_add(parallel::distributed::Vector<value_type> &dst,
                                        double const                              &eval_time);

  void rhs_ppe_div_term_convective_term_add(parallel::distributed::Vector<value_type>       &dst,
                                            parallel::distributed::Vector<value_type> const &src);

  // rhs pressure
  void rhs_ppe_nbc_add (parallel::distributed::Vector<value_type> &dst,
                        double const                              &evaluation_time);

  // rhs pressure: Neumann BC convective term
  void rhs_ppe_convective_add (parallel::distributed::Vector<value_type>       &dst,
                               const parallel::distributed::Vector<value_type> &src) const;

  // rhs pressure: Neumann BC viscous term
  void rhs_ppe_viscous_add(parallel::distributed::Vector<value_type>       &dst,
                           const parallel::distributed::Vector<value_type> &src) const;

  // viscous step
  unsigned int solve_viscous (parallel::distributed::Vector<value_type>       &dst,
                              parallel::distributed::Vector<value_type> const &src,
                              double const                                    &scaling_factor_time_derivative_term);

  FEParameters<dim> const & get_fe_parameters() const
  {
    return this->fe_param;
  }

protected:
  FEParameters<dim> fe_param;
  HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule,value_type> helmholtz_operator;
  std_cxx11::shared_ptr<PreconditionerBase<value_type> > helmholtz_preconditioner;

private:
  // Helmholtz solver
  std_cxx11::shared_ptr<IterativeSolverBase<parallel::distributed::Vector<value_type> > > helmholtz_solver;

  // Implicit solution of convective step
  parallel::distributed::Vector<value_type> velocity_linearization;
  parallel::distributed::Vector<value_type> temp;
  parallel::distributed::Vector<value_type> const * sum_alphai_ui;

  // implicit solution of convective step
  std_cxx11::shared_ptr<InverseMassMatrixPreconditioner<dim,fe_degree,value_type> > preconditioner_convective_problem;
  std_cxx11::shared_ptr<IterativeSolverBase<parallel::distributed::Vector<value_type> > > linear_solver;
  std_cxx11::shared_ptr<NewtonSolver<parallel::distributed::Vector<value_type>,
                                     DGNavierStokesDualSplitting<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>,
                                     DGNavierStokesDualSplitting<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>,
                                     IterativeSolverBase<parallel::distributed::Vector<value_type> > > > newton_solver;


  // rhs pressure Poisson equation: Velocity divergence term
  VelocityDivergenceConvectiveTerm<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule, value_type> velocity_divergence_convective_term;

  // pressure Neumann BC
  PressureNeumannBCConvectiveTerm<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule, value_type> pressure_nbc_convective_term;
  PressureNeumannBCViscousTerm<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule, value_type> pressure_nbc_viscous_term;

  double evaluation_time;
  double scaling_factor_time_derivative_term;

  // setup of solvers
  void setup_convective_solver();
  virtual void setup_pressure_poisson_solver(double const time_step_size);
  void setup_helmholtz_solver(double const &scaling_factor_time_derivative_term);
  virtual void setup_helmholtz_preconditioner();

  // rhs pressure: boundary conditions for intermediate velocity u_hat
  // body force term
  void local_rhs_ppe_div_term_body_forces (const MatrixFree<dim,value_type>                &data,
                                           parallel::distributed::Vector<value_type>       &dst,
                                           const parallel::distributed::Vector<value_type> &src,
                                           const std::pair<unsigned int,unsigned int>      &cell_range) const;

  void local_rhs_ppe_div_term_body_forces_face (const MatrixFree<dim,value_type>                 &data,
                                                parallel::distributed::Vector<value_type>        &dst,
                                                const parallel::distributed::Vector<value_type>  &src,
                                                const std::pair<unsigned int,unsigned int>       &face_range) const;

  void local_rhs_ppe_div_term_body_forces_boundary_face(const MatrixFree<dim,value_type>                 &data,
                                                        parallel::distributed::Vector<value_type>        &dst,
                                                        const parallel::distributed::Vector<value_type>  &src,
                                                        const std::pair<unsigned int,unsigned int>       &face_range) const;

  // rhs pressure: NBC term
  void local_rhs_ppe_nbc_add (const MatrixFree<dim,value_type>                &data,
                              parallel::distributed::Vector<value_type>       &dst,
                              const parallel::distributed::Vector<value_type> &src,
                              const std::pair<unsigned int,unsigned int>      &cell_range) const;

  void local_rhs_ppe_nbc_add_face (const MatrixFree<dim,value_type>                 &data,
                                   parallel::distributed::Vector<value_type>        &dst,
                                   const parallel::distributed::Vector<value_type>  &src,
                                   const std::pair<unsigned int,unsigned int>       &face_range) const;

  void local_rhs_ppe_nbc_add_boundary_face(const MatrixFree<dim,value_type>                 &data,
                                           parallel::distributed::Vector<value_type>        &dst,
                                           const parallel::distributed::Vector<value_type>  &src,
                                           const std::pair<unsigned int,unsigned int>       &face_range) const;

};

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
setup_solvers (double const &time_step_size,
               double const &scaling_factor_time_derivative_term)
{
  ConditionalOStream pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
  pcout << std::endl << "Setup solvers ..." << std::endl;

  // initialize vectors that are needed by the nonlinear solver
  if(this->param.equation_type == EquationType::NavierStokes &&
     this->param.treatment_of_convective_term == TreatmentOfConvectiveTerm::Implicit)
  {
    setup_convective_solver();
  }

  this->setup_pressure_poisson_solver(time_step_size);

  this->setup_projection_solver();

  setup_helmholtz_solver(scaling_factor_time_derivative_term);

  pcout << std::endl << "... done!" << std::endl;
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
setup_convective_solver ()
{
  this->initialize_vector_velocity(temp);
  this->initialize_vector_velocity(velocity_linearization);

  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::QuadratureSelector QuadratureSelector;

  // preconditioner implicit convective step
  preconditioner_convective_problem.reset(new InverseMassMatrixPreconditioner<dim,fe_degree,value_type>
     (this->data,
      static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::velocity),
      static_cast<typename std::underlying_type<QuadratureSelector>::type>(QuadratureSelector::velocity)));

  // linear solver (GMRES)
  GMRESSolverData solver_data;
  solver_data.max_iter = this->param.max_iter_linear_convective;
  solver_data.solver_tolerance_abs = this->param.abs_tol_linear_convective;
  solver_data.solver_tolerance_rel = this->param.rel_tol_linear_convective;
  solver_data.right_preconditioning = this->param.use_right_preconditioning_convective;
  solver_data.max_n_tmp_vectors = this->param.max_n_tmp_vectors_convective;

  // always use inverse mass matrix preconditioner
  solver_data.use_preconditioner = true;

  // setup linear solver
  linear_solver.reset(new GMRESSolver<DGNavierStokesDualSplitting<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>,
                                      InverseMassMatrixPreconditioner<dim,fe_degree,value_type>,
                                      parallel::distributed::Vector<value_type> >
      (*this,*preconditioner_convective_problem,solver_data));

  // setup Newton solver
  newton_solver.reset(new NewtonSolver<parallel::distributed::Vector<value_type>,
                                       DGNavierStokesDualSplitting<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>,
                                       DGNavierStokesDualSplitting<dim,fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>,
                                       IterativeSolverBase<parallel::distributed::Vector<value_type> > >
     (this->param.newton_solver_data_convective,*this,*this,*linear_solver));
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
setup_pressure_poisson_solver (const double time_step_size)
{
  // Call setup function of base class
  DGNavierStokesProjectionMethods<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::setup_pressure_poisson_solver(time_step_size);

  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;

  // RHS PPE: Velocity divergence term
  VelocityDivergenceConvectiveTermData<dim> velocity_divergence_convective_data;
  velocity_divergence_convective_data.dof_index_velocity = static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::velocity);
  velocity_divergence_convective_data.dof_index_pressure =  static_cast<typename std::underlying_type<DofHandlerSelector>::type >(DofHandlerSelector::pressure);
  velocity_divergence_convective_data.bc = this->boundary_descriptor_velocity;

  velocity_divergence_convective_term.initialize(this->data,velocity_divergence_convective_data);

  // Pressure NBC: Convective term
  PressureNeumannBCConvectiveTermData<dim> pressure_nbc_convective_data;
  pressure_nbc_convective_data.dof_index_velocity = static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::velocity);
  pressure_nbc_convective_data.dof_index_pressure =  static_cast<typename std::underlying_type<DofHandlerSelector>::type >(DofHandlerSelector::pressure);
  pressure_nbc_convective_data.bc = this->boundary_descriptor_pressure;

  pressure_nbc_convective_term.initialize(this->data,pressure_nbc_convective_data);

  // Pressure NBC: Viscous term
  PressureNeumannBCViscousTermData<dim> pressure_nbc_viscous_data;
  pressure_nbc_viscous_data.dof_index_velocity = static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::velocity);
  pressure_nbc_viscous_data.dof_index_pressure =  static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::pressure);
  pressure_nbc_viscous_data.bc = this->boundary_descriptor_pressure;

  pressure_nbc_viscous_term.initialize(this->data,pressure_nbc_viscous_data,this->viscous_operator);
}


template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
setup_helmholtz_solver (double const &scaling_factor_time_derivative_term)
{
  // 1. Setup Helmholtz operator
  HelmholtzOperatorData<dim> helmholtz_operator_data;

  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;

  helmholtz_operator_data.dof_index = static_cast<typename std::underlying_type<DofHandlerSelector>::type >(DofHandlerSelector::velocity);

  // always unsteady problem
  helmholtz_operator_data.unsteady_problem = true;

  helmholtz_operator.initialize(this->data,helmholtz_operator_data,this->mass_matrix_operator,this->viscous_operator);

  // set scaling factor time derivative term!
  helmholtz_operator.set_scaling_factor_time_derivative_term(scaling_factor_time_derivative_term);

  // 2. Setup Helmholtz preconditioner
  setup_helmholtz_preconditioner();

  // 3. Setup Helmholtz solver
  if(this->param.solver_viscous == SolverViscous::PCG)
  {
    // setup solver data
    CGSolverData solver_data;
    // use default value of max_iter
    solver_data.solver_tolerance_abs = this->param.abs_tol_viscous;
    solver_data.solver_tolerance_rel = this->param.rel_tol_viscous;
    // default value of use_preconditioner = false
    if(this->param.preconditioner_viscous == PreconditionerViscous::PointJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::BlockJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::InverseMassMatrix ||
       this->param.preconditioner_viscous == PreconditionerViscous::GeometricMultigrid)
    {
      solver_data.use_preconditioner = true;
    }
    solver_data.update_preconditioner = this->param.update_preconditioner_viscous;

    // setup helmholtz solver
    helmholtz_solver.reset(new CGSolver<HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule,value_type>,
                                        PreconditionerBase<value_type>,
                                        parallel::distributed::Vector<value_type> >
       (helmholtz_operator,
        *helmholtz_preconditioner,
        solver_data));
  }
  else if(this->param.solver_viscous == SolverViscous::GMRES)
  {
    // setup solver data
    GMRESSolverData solver_data;
    // use default value of max_iter
    solver_data.solver_tolerance_abs = this->param.abs_tol_viscous;
    solver_data.solver_tolerance_rel = this->param.rel_tol_viscous;
    // use default value of right_preconditioning
    // use default value of max_n_tmp_vectors
    // use default value of compute_eigenvalues
    solver_data.update_preconditioner = this->param.update_preconditioner_viscous;

    // default value of use_preconditioner = false
    if(this->param.preconditioner_viscous == PreconditionerViscous::PointJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::BlockJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::InverseMassMatrix ||
       this->param.preconditioner_viscous == PreconditionerViscous::GeometricMultigrid)
    {
      solver_data.use_preconditioner = true;
    }

    // setup helmholtz solver
    helmholtz_solver.reset(new GMRESSolver<HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule,value_type>,
                                           PreconditionerBase<value_type>,
                                           parallel::distributed::Vector<value_type> >
       (helmholtz_operator,
        *helmholtz_preconditioner,
        solver_data));
  }
  else if(this->param.solver_viscous == SolverViscous::FGMRES)
  {
    FGMRESSolverData solver_data;
    // use default value of max_iter
    solver_data.solver_tolerance_abs = this->param.abs_tol_viscous;
    solver_data.solver_tolerance_rel = this->param.rel_tol_viscous;
    // use default value of max_n_tmp_vectors
    solver_data.update_preconditioner = this->param.update_preconditioner_viscous;

    if(this->param.preconditioner_viscous == PreconditionerViscous::PointJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::BlockJacobi ||
       this->param.preconditioner_viscous == PreconditionerViscous::InverseMassMatrix ||
       this->param.preconditioner_viscous == PreconditionerViscous::GeometricMultigrid)
    {
      solver_data.use_preconditioner = true;
    }

    helmholtz_solver.reset(new FGMRESSolver<HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule,value_type>,
                                            PreconditionerBase<value_type>,
                                            parallel::distributed::Vector<value_type>  >
        (helmholtz_operator,
         *helmholtz_preconditioner,
         solver_data));
  }
  else
  {
    AssertThrow(this->param.solver_viscous == SolverViscous::PCG ||
                this->param.solver_viscous == SolverViscous::GMRES ||
                this->param.solver_viscous == SolverViscous::FGMRES,
                ExcMessage("Specified Viscous Solver not implemented - possibilities are PCG and GMRES"));
  }
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
setup_helmholtz_preconditioner ()
{
  if(this->param.preconditioner_viscous == PreconditionerViscous::InverseMassMatrix)
  {
    typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;
    typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::QuadratureSelector QuadratureSelector;

    helmholtz_preconditioner.reset(new InverseMassMatrixPreconditioner<dim,fe_degree,value_type>
       (this->data,
        static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::velocity),
        static_cast<typename std::underlying_type<QuadratureSelector>::type>(QuadratureSelector::velocity)));
  }
  else if(this->param.preconditioner_viscous == PreconditionerViscous::PointJacobi)
  {
    helmholtz_preconditioner.reset(new JacobiPreconditioner<value_type,
        HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule, value_type> >(helmholtz_operator));
  }
  else if(this->param.preconditioner_viscous == PreconditionerViscous::BlockJacobi)
  {
    helmholtz_preconditioner.reset(new BlockJacobiPreconditioner<value_type,
        HelmholtzOperator<dim,fe_degree,fe_degree_xwall,xwall_quad_rule, value_type> >(helmholtz_operator));
  }
  else if(this->param.preconditioner_viscous == PreconditionerViscous::GeometricMultigrid)
  {
    MultigridData mg_data;
    mg_data = this->param.multigrid_data_viscous;

    // use single precision for multigrid
    typedef float Number;

    typedef MyMultigridPreconditionerVelocityDiffusion<dim,value_type,
        HelmholtzOperator<dim, fe_degree, fe_degree_xwall, xwall_quad_rule, Number>,
        HelmholtzOperator<dim, fe_degree, fe_degree_xwall, xwall_quad_rule, value_type> > MULTIGRID;

    helmholtz_preconditioner.reset(new MULTIGRID());

    std_cxx11::shared_ptr<MULTIGRID> mg_preconditioner = std::dynamic_pointer_cast<MULTIGRID>(helmholtz_preconditioner);

    mg_preconditioner->initialize(mg_data,
                                  this->dof_handler_u,
                                  this->mapping,
                                  helmholtz_operator,
                                  this->periodic_face_pairs);
  }
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
solve_nonlinear_convective_problem (parallel::distributed::Vector<value_type>       &dst,
                                    parallel::distributed::Vector<value_type> const &sum_alphai_ui,
                                    double const                                    &eval_time,
                                    double const                                    &scaling_factor_mass_matrix_term,
                                    unsigned int                                    &newton_iterations,
                                    unsigned int                                    &linear_iterations)
{
  // Set sum_alphai_ui, this variable is used when evaluating the nonlinear residual
  this->sum_alphai_ui = &sum_alphai_ui;

  // set evaluation time for both the linear and the nonlinear operator (=DGNavierStokesDualSplitting)
  evaluation_time = eval_time;

  // set scaling_factor_time_derivative term for both the linear and the nonlinear operator
  // (=DGNavierStokesDualSplitting)
  scaling_factor_time_derivative_term = scaling_factor_mass_matrix_term;

  // solve nonlinear problem
  newton_solver->solve(dst,newton_iterations,linear_iterations);

  // Reset sum_alphai_ui
  this->sum_alphai_ui = nullptr;
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
evaluate_nonlinear_residual(parallel::distributed::Vector<value_type>       &dst,
                            const parallel::distributed::Vector<value_type> &src)
{
  if(this->param.right_hand_side == true)
  {
    this->body_force_operator.evaluate(dst,evaluation_time);
    // shift body force term to the left-hand side of the equation
    dst *= -1.0;
  }
  else // right_hand_side == false
  {
    // set dst to zero. This is necessary since the subsequent operators
    // call functions of type ..._add
    dst = 0.0;
  }

  // temp, src, sum_alphai_ui have the same number of blocks
  temp.equ(scaling_factor_time_derivative_term,src);
  temp.add(-1.0,*sum_alphai_ui);

  this->mass_matrix_operator.apply_add(dst,temp);

  this->convective_operator.evaluate_add(dst,src,evaluation_time);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
vmult (parallel::distributed::Vector<value_type>       &dst,
       parallel::distributed::Vector<value_type> const &src) const
{
  apply_linearized_convective_problem(dst,src);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
apply_linearized_convective_problem (parallel::distributed::Vector<value_type>       &dst,
                                     parallel::distributed::Vector<value_type> const &src) const
{
  this->mass_matrix_operator.apply(dst,src);

  dst *= scaling_factor_time_derivative_term;

  this->convective_operator.apply_linearized_add(dst,
                                                 src,
                                                 &velocity_linearization,
                                                 evaluation_time);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
evaluate_convective_term_and_apply_inverse_mass_matrix (parallel::distributed::Vector<value_type>       &dst,
                                                        parallel::distributed::Vector<value_type> const &src,
                                                        value_type const                                evaluation_time) const
{
  this->convective_operator.evaluate(dst,src,evaluation_time);
  this->inverse_mass_matrix_operator->apply(dst,dst);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
evaluate_body_force_and_apply_inverse_mass_matrix(parallel::distributed::Vector<value_type>  &dst,
                                                  const value_type                           evaluation_time) const
{
  this->body_force_operator.evaluate(dst,evaluation_time);

  this->inverse_mass_matrix_operator->apply(dst,dst);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
apply_velocity_divergence_term(parallel::distributed::Vector<value_type>        &dst,
                               const parallel::distributed::Vector<value_type>  &src) const
{
  this->divergence_operator.apply(dst,src);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_velocity_divergence_term(parallel::distributed::Vector<value_type>        &dst,
                             double const                                     &evaluation_time) const
{
  this->divergence_operator.rhs(dst,evaluation_time);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_ppe_div_term_body_forces_add(parallel::distributed::Vector<value_type> &dst,
                                 double const                              &eval_time)
{
  evaluation_time = eval_time;

  parallel::distributed::Vector<value_type> src_dummy;
  this->data.loop (&DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_div_term_body_forces,
                   &DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_div_term_body_forces_face,
                   &DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_div_term_body_forces_boundary_face,
                   this, dst, src_dummy);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_div_term_body_forces (const MatrixFree<dim,value_type>                 &,
                                    parallel::distributed::Vector<value_type>        &,
                                    const parallel::distributed::Vector<value_type>  &,
                                    const std::pair<unsigned int,unsigned int>       &) const
{

}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_div_term_body_forces_face (const MatrixFree<dim,value_type>                &,
                                         parallel::distributed::Vector<value_type>       &,
                                         const parallel::distributed::Vector<value_type> &,
                                         const std::pair<unsigned int,unsigned int>      &) const
{

}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_div_term_body_forces_boundary_face (const MatrixFree<dim,value_type>                 &data,
                                                  parallel::distributed::Vector<value_type>        &dst,
                                                  const parallel::distributed::Vector<value_type>  &,
                                                  const std::pair<unsigned int,unsigned int>       &face_range) const
{
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::QuadratureSelector QuadratureSelector;

  unsigned int dof_index_pressure = static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::pressure);
  unsigned int quad_index_pressure = static_cast<typename std::underlying_type<QuadratureSelector>::type>(QuadratureSelector::pressure);

  FEFaceEvaluation<dim,fe_degree_p,fe_degree_p+1,1,value_type> fe_eval(data,true,dof_index_pressure,quad_index_pressure);

  for(unsigned int face=face_range.first; face<face_range.second; face++)
  {
    fe_eval.reinit (face);

    typename std::map<types::boundary_id,std_cxx11::shared_ptr<Function<dim> > >::iterator it;
    types::boundary_id boundary_id = data.get_boundary_indicator(face);

    for(unsigned int q=0;q<fe_eval.n_q_points;++q)
    {
      it = this->boundary_descriptor_pressure->dirichlet_bc.find(boundary_id);
      if(it != this->boundary_descriptor_pressure->dirichlet_bc.end())
      {
        Point<dim,VectorizedArray<value_type> > q_points = fe_eval.quadrature_point(q);

        // evaluate right-hand side
        Tensor<1,dim,VectorizedArray<value_type> > rhs;
        evaluate_vectorial_function(rhs,this->field_functions->right_hand_side,q_points,evaluation_time);

        VectorizedArray<value_type> flux_times_normal = rhs * fe_eval.get_normal_vector(q);
        // minus sign is introduced here which allows to call a function of type ...add()
        // and avoids a scaling of the resulting vector by the factor -1.0
        fe_eval.submit_value(-flux_times_normal,q);
      }

      it = this->boundary_descriptor_pressure->neumann_bc.find(boundary_id);
      if (it != this->boundary_descriptor_pressure->neumann_bc.end())
      {
        // do nothing on Neumann boundaries
        VectorizedArray<value_type> zero = make_vectorized_array<value_type>(0.0);
        fe_eval.submit_value(zero,q);
      }
    }
    fe_eval.integrate(true,false);
    fe_eval.distribute_local_to_global(dst);
  }
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_ppe_div_term_convective_term_add(parallel::distributed::Vector<value_type>       &dst,
                                     parallel::distributed::Vector<value_type> const &src)
{
  velocity_divergence_convective_term.calculate(dst,src);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_ppe_nbc_add (parallel::distributed::Vector<value_type> &dst,
                 double const                              &eval_time)
{
  evaluation_time = eval_time;

  parallel::distributed::Vector<value_type> src_dummy;
  this->data.loop (&DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_nbc_add,
                   &DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_nbc_add_face,
                   &DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::local_rhs_ppe_nbc_add_boundary_face,
                   this, dst, src_dummy);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_nbc_add (const MatrixFree<dim,value_type>                 &,
                       parallel::distributed::Vector<value_type>        &,
                       const parallel::distributed::Vector<value_type>  &,
                       const std::pair<unsigned int,unsigned int>       &) const
{

}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_nbc_add_face (const MatrixFree<dim,value_type>                &,
                            parallel::distributed::Vector<value_type>       &,
                            const parallel::distributed::Vector<value_type> &,
                            const std::pair<unsigned int,unsigned int>      &) const
{

}


template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
local_rhs_ppe_nbc_add_boundary_face (const MatrixFree<dim,value_type>                 &data,
                                     parallel::distributed::Vector<value_type>        &dst,
                                     const parallel::distributed::Vector<value_type>  &,
                                     const std::pair<unsigned int,unsigned int>       &face_range) const
{
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::DofHandlerSelector DofHandlerSelector;
  typedef typename DGNavierStokesBase<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::QuadratureSelector QuadratureSelector;

  unsigned int dof_index_pressure = static_cast<typename std::underlying_type<DofHandlerSelector>::type>(DofHandlerSelector::pressure);
  unsigned int quad_index_pressure = static_cast<typename std::underlying_type<QuadratureSelector>::type>(QuadratureSelector::pressure);

  FEFaceEvaluation<dim,fe_degree_p,fe_degree_p+1,1,value_type> fe_eval(data,true,dof_index_pressure,quad_index_pressure);

  for(unsigned int face=face_range.first; face<face_range.second; face++)
  {
    fe_eval.reinit (face);

    typename std::map<types::boundary_id,std_cxx11::shared_ptr<Function<dim> > >::iterator it;
    types::boundary_id boundary_id = data.get_boundary_indicator(face);

    for(unsigned int q=0;q<fe_eval.n_q_points;++q)
    {
      it = this->boundary_descriptor_pressure->dirichlet_bc.find(boundary_id);
      if(it != this->boundary_descriptor_pressure->dirichlet_bc.end())
      {
        Point<dim,VectorizedArray<value_type> > q_points = fe_eval.quadrature_point(q);

        // evaluate right-hand side
        Tensor<1,dim,VectorizedArray<value_type> > rhs;
        evaluate_vectorial_function(rhs,this->field_functions->right_hand_side,q_points,evaluation_time);

        // evaluate boundary condition
        Tensor<1,dim,VectorizedArray<value_type> > dudt;
        evaluate_vectorial_function(dudt,it->second,q_points,evaluation_time);

        Tensor<1,dim,VectorizedArray<value_type> > normal = fe_eval.get_normal_vector(q);
        VectorizedArray<value_type> h;

        h = - normal * (dudt - rhs);

        fe_eval.submit_value(h,q);
      }

      it = this->boundary_descriptor_pressure->neumann_bc.find(boundary_id);
      if (it != this->boundary_descriptor_pressure->neumann_bc.end())
      {
        VectorizedArray<value_type> zero = make_vectorized_array<value_type>(0.0);
        fe_eval.submit_value(zero,q);
      }
    }
    fe_eval.integrate(true,false);
    fe_eval.distribute_local_to_global(dst);
  }
}


template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_ppe_convective_add (parallel::distributed::Vector<value_type>       &dst,
                        const parallel::distributed::Vector<value_type> &src) const
{
  pressure_nbc_convective_term.calculate(dst,src);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
void DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
rhs_ppe_viscous_add(parallel::distributed::Vector<value_type>       &dst,
                    const parallel::distributed::Vector<value_type> &src) const
{
  pressure_nbc_viscous_term.calculate(dst,src);
}

template<int dim, int fe_degree, int fe_degree_p, int fe_degree_xwall, int xwall_quad_rule>
unsigned int DGNavierStokesDualSplitting<dim, fe_degree, fe_degree_p, fe_degree_xwall, xwall_quad_rule>::
solve_viscous (parallel::distributed::Vector<value_type>       &dst,
               parallel::distributed::Vector<value_type> const &src,
               double const                                    &factor)
{
  // Update Helmholtz operator
  helmholtz_operator.set_scaling_factor_time_derivative_term(factor);
  // viscous_operator.set_constant_viscosity(viscosity);
  // viscous_operator.set_variable_viscosity(viscosity);

  unsigned int n_iter = helmholtz_solver->solve(dst,src);

  return n_iter;
}

#endif /* INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_DG_NAVIER_STOKES_DUAL_SPLITTING_H_ */
