#ifndef INCLUDE_EXADG_AERO_ACOUSTIC_SOURCE_TERM_CALCULATOR_H_
#define INCLUDE_EXADG_AERO_ACOUSTIC_SOURCE_TERM_CALCULATOR_H_

#include <exadg/matrix_free/integrators.h>
#include <exadg/utilities/lazy_ptr.h>

#include <deal.II/matrix_free/fe_remote_evaluation.h>

template<int dim, typename Number, typename VectorizedArrayType>
dealii::FERemoteEvaluationCommunicator<dim>
compute_remote_communicator_cells_point_to_point_interpolation(
  const dealii::MatrixFree<dim, Number, VectorizedArrayType> & matrix_free_dst,
  const dealii::MatrixFree<dim, Number, VectorizedArrayType> & matrix_free_src,
  const unsigned int                                           quad_no_dst = 0,
  const unsigned int                                           dof_no_dst  = 0,
  const unsigned int                                           dof_no_src  = 0,
  const double                                                 tolerance   = 1e-9)
{
  const auto & dof_handler_src = matrix_free_src.get_dof_handler(dof_no_src);
  const auto & tria_src        = dof_handler_src.get_triangulation();
  const auto & mapping_src     = *matrix_free_src.get_mapping_info().mapping;

  dealii::FERemoteCommunicationObjectEntityBatches<dim> comm_object;

  std::vector<unsigned int> global_quadrature_sizes(matrix_free_dst.n_cell_batches(),
                                                    dealii::numbers::invalid_unsigned_int);

  auto rpe =
    std::make_shared<dealii::Utilities::MPI::RemotePointEvaluation<dim>>(tolerance, false, 0);

  std::vector<std::pair<unsigned int, unsigned int>> cell_batch_id_n_cells;

  // Points that are searched by rpe.
  std::vector<dealii::Point<dim>> points;

  // Temporarily set up FEFaceEvaluation to access the quadrature points
  // at the faces on the non-matching interface.
  dealii::FEEvaluation<dim, -1, 0, 1, Number> phi(matrix_free_dst, dof_no_dst, quad_no_dst);

  std::pair<unsigned int, unsigned int> cell_batch_range{0, matrix_free_dst.n_cell_batches()};

  // Iterate over the boundary faces.
  for(unsigned int cell = 0; cell < matrix_free_dst.n_cell_batches(); ++cell)
  {
    phi.reinit(cell);

    // If @c face is on the current side of the non-matching
    // interface. Add the face batch ID and the number of faces in
    // the batch to the corresponding data structure.
    const unsigned int n_cells = matrix_free_dst.n_active_entries_per_cell_batch(cell);
    cell_batch_id_n_cells.emplace_back(std::make_pair(cell, n_cells));

    // Append the quadrature points to the points we need to search
    // for.
    for(unsigned int v = 0; v < n_cells; ++v)
    {
      for(unsigned int q : phi.quadrature_point_indices())
      {
        const auto         point = phi.quadrature_point(q);
        dealii::Point<dim> temp;
        for(unsigned int i = 0; i < dim; ++i)
          temp[i] = point[i][v];

        points.push_back(temp);
      }
    }

    // Insert the quadrature size into the global vector.
    // First check that each face is only considered once.
    Assert(global_quadrature_sizes[cell] == numbers::invalid_unsigned_int,
           ExcMessage("Quadrature for given face already provided."));

    global_quadrature_sizes[cell] = phi.n_q_points;
  }

  // Reinit RPE and ensure all points are found.
  rpe->reinit(points, tria_src, mapping_src);
  Assert(rpe->all_points_found(), ExcMessage("Not all remote points found."));

  comm_object.batch_id_n_entities = cell_batch_id_n_cells;
  comm_object.rpe                 = rpe;

  dealii::FERemoteEvaluationCommunicator<dim> remote_communicator;

  remote_communicator.reinit_cells(comm_object, cell_batch_range, global_quadrature_sizes);

  return remote_communicator;
}



namespace ExaDG
{
namespace AeroAcoustic
{
struct FeedbackTermCalculatorData
{
  unsigned int dof_index;
  unsigned int quad_index;
  unsigned int dof_index_acoustic;
  double       density;
};

template<int dim, typename Number>
class FeedbackTermCalculator
{
  using This                 = FeedbackTermCalculator<dim, Number>;
  using VectorType           = dealii::LinearAlgebra::distributed::Vector<Number>;
  using CellIntegratorVector = CellIntegrator<dim, dim, Number>;
  using RemoteCellIntegratorVector =
    dealii::FERemoteEvaluation<dim, dim, dealii::VectorizedArray<Number>>;


public:
  FeedbackTermCalculator() : matrix_free(nullptr)
  {
  }

  void
  setup(dealii::MatrixFree<dim, Number> const & matrix_free_fluid,
        dealii::MatrixFree<dim, Number> const & matrix_free_acoustic,
        FeedbackTermCalculatorData const &      data_in)
  {
    matrix_free = &matrix_free_fluid;
    data        = data_in;

    communicator =
      compute_remote_communicator_cells_point_to_point_interpolation(matrix_free_fluid,
                                                                     matrix_free_acoustic,
                                                                     data_in.quad_index,
                                                                     data_in.dof_index,
                                                                     data_in.dof_index_acoustic);

    acoustic_particle_velocity =
      std::make_unique<RemoteCellIntegratorVector>(communicator,
                                                   matrix_free->get_dof_handler(data_in.dof_index));
  }

  void
  evaluate_integrate(VectorType &       dst,
                     VectorType const & velocity_cfd,
                     VectorType const & velocity_acoustic)
  {
    dst.zero_out_ghost_values();

    acoustic_particle_velocity->gather_evaluate(velocity_acoustic, dealii::EvaluationFlags::values);

    matrix_free->cell_loop(&This::compute_feedback_term, this, dst, velocity_cfd, true);
  }

  template<typename T1, typename T2>
  static inline DEAL_II_ALWAYS_INLINE //
    T2
    cross_product(T1 const & omega, T2 const & u_a)
  {
    static_assert(dim == 3 || dim == 2, "feedback term only possible for dimensions 2 and 3");

    if constexpr(dim == 3)
    {
      return dealii::cross_product_3d(omega, u_a);
    }

    if constexpr(dim == 2)
    {
      // vorticity is a scalar (stored in component 0)
      // cross_product_2d() rotates vector clockwise, we need it counterclockwise since
      // [omega,omega,omega]^T x [u1, u2, u3] = [-omega*u2, omega*u1, ...]
      return omega[0] * (-1.0 * dealii::cross_product_2d(u_a));
    }

    return {};
  }


  void
  compute_feedback_term(dealii::MatrixFree<dim, Number> const &       matrix_free_in,
                        VectorType &                                  dst,
                        VectorType const &                            velocity_cfd,
                        std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    // − (∇ × u ic ) × u a, we solve for rho*ua
    CellIntegratorVector feedback_term(matrix_free_in, data.dof_index, data.quad_index);
    auto                 u_acoustic = acoustic_particle_velocity->get_data_accessor();

    double m_rho_inv = -1.0 / data.density;

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      u_acoustic.reinit(cell);

      feedback_term.reinit(cell);
      feedback_term.gather_evaluate(velocity_cfd, dealii::EvaluationFlags::gradients);

      for(unsigned int q = 0; q < feedback_term.n_q_points; ++q)
      {
        auto const u_a   = u_acoustic.get_value(q);
        auto const omega = feedback_term.get_curl(q);
        feedback_term.submit_value(m_rho_inv * cross_product(omega, u_a), q);
      }

      feedback_term.integrate_scatter(dealii::EvaluationFlags::values, dst);
    }
  }

  dealii::MatrixFree<dim, Number> const * matrix_free;

  FeedbackTermCalculatorData data;

  dealii::FERemoteEvaluationCommunicator<dim> communicator;
  std::unique_ptr<RemoteCellIntegratorVector> acoustic_particle_velocity;
};


template<int dim>
struct SourceTermCalculatorData
{
  unsigned int dof_index_pressure;
  unsigned int dof_index_velocity;
  unsigned int quad_index;

  // density of the underlying fluid.
  double density;

  // speed of sound of the fluid.
  double speed_of_sound;

  // use material or partial temporal derivative of pressure as source term.
  bool consider_convection;

  // function if blend in is required.
  bool                                                  blend_in;
  std::shared_ptr<Utilities::SpatialAwareFunction<dim>> blend_in_function;
};

/**
 * A class that knows how to compute the aeroacoustic source term on the fluid mesh.
 * evaluate_integrate() computes and integrates the source term on the fluid mesh.
 *
 * The aeroacoustic source term f is definded as:
 * f = - rho * (dp/dt + u * grad(p)).
 * The scaling factor rho has to be used since the pressure of the incompressible
 * module is a kinematic pressure. Using consider_convection=false
 * f = -rho * (dp/dt).
 */
template<int dim, typename Number>
class SourceTermCalculator
{
  using This       = SourceTermCalculator<dim, Number>;
  using VectorType = dealii::LinearAlgebra::distributed::Vector<Number>;

  using CellIntegratorScalar = CellIntegrator<dim, 1, Number>;
  using CellIntegratorVector = CellIntegrator<dim, dim, Number>;

  using scalar = dealii::VectorizedArray<Number>;
  using qpoint = dealii::Point<dim, dealii::VectorizedArray<Number>>;

public:
  SourceTermCalculator() : matrix_free(nullptr), time(std::numeric_limits<double>::min())
  {
  }

  void
  setup(dealii::MatrixFree<dim, Number> const & matrix_free_in,
        SourceTermCalculatorData<dim> const &   data_in)
  {
    matrix_free = &matrix_free_in;
    data        = data_in;
  }

  void
  evaluate_integrate(VectorType &            dst,
                     dealii::Function<dim> & analytical_source_term,
                     double const            evaluation_time)
  {
    time = evaluation_time;

    dst.zero_out_ghost_values();

    analytical_source_term.set_time(time);

    matrix_free->cell_loop(&This::compute_source_term, this, dst, analytical_source_term, true);
  }


  void
  evaluate_integrate(VectorType &       dst,
                     VectorType const & velocity_cfd_in,
                     VectorType const & pressure_cfd_in,
                     VectorType const & pressure_cfd_time_derivative_in,
                     double const       evaluation_time)
  {
    time = evaluation_time;

    dst.zero_out_ghost_values();

    if(data.consider_convection)
    {
      velocity_cfd.reset(velocity_cfd_in);
      velocity_cfd->update_ghost_values();

      pressure_cfd.reset(pressure_cfd_in);
      pressure_cfd->update_ghost_values();
    }

    matrix_free->cell_loop(
      &This::compute_source_term, this, dst, pressure_cfd_time_derivative_in, true);
  }

private:
  void
  compute_source_term(dealii::MatrixFree<dim, Number> const &       matrix_free_in,
                      VectorType &                                  dst,
                      VectorType const &                            dp_cfd_dt,
                      std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar dpdt(matrix_free_in, data.dof_index_pressure, data.quad_index);
    CellIntegratorScalar p(matrix_free_in, data.dof_index_pressure, data.quad_index);
    CellIntegratorVector u(matrix_free_in, data.dof_index_velocity, data.quad_index);

    auto get_scaling_factor = get_scaling_function();

    Number const m_rho_c2 =
      static_cast<Number>(-data.density / (data.speed_of_sound * data.speed_of_sound));

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      dpdt.reinit(cell);
      dpdt.gather_evaluate(dp_cfd_dt, dealii::EvaluationFlags::values);

      if(data.consider_convection)
      {
        p.reinit(cell);
        p.gather_evaluate(*pressure_cfd, dealii::EvaluationFlags::gradients);
        u.reinit(cell);
        u.gather_evaluate(*velocity_cfd, dealii::EvaluationFlags::values);

        for(unsigned int q = 0; q < dpdt.n_q_points; ++q)
        {
          scalar flux = m_rho_c2 * (dpdt.get_value(q) + u.get_value(q) * p.get_gradient(q));

          if(data.blend_in)
            flux *= get_scaling_factor(dpdt.quadrature_point(q));

          dpdt.submit_value(flux, q);
        }
      }
      else
      {
        for(unsigned int q = 0; q < dpdt.n_q_points; ++q)
        {
          scalar flux = m_rho_c2 * dpdt.get_value(q);

          if(data.blend_in)
            flux *= get_scaling_factor(dpdt.quadrature_point(q));

          dpdt.submit_value(flux, q);
        }
      }

      dpdt.integrate_scatter(dealii::EvaluationFlags::values, dst);
    }
  }

  void
  compute_source_term(dealii::MatrixFree<dim, Number> const &       matrix_free_in,
                      VectorType &                                  dst,
                      dealii::Function<dim> const &                 analytical_source_term,
                      std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar dpdt(matrix_free_in, data.dof_index_pressure, data.quad_index);

    auto get_scaling_factor = get_scaling_function();

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      dpdt.reinit(cell);

      for(unsigned int q = 0; q < dpdt.n_q_points; ++q)
      {
        scalar flux = FunctionEvaluator<0, dim, Number>::value(analytical_source_term,
                                                               dpdt.quadrature_point(q));

        if(data.blend_in)
          flux *= get_scaling_factor(dpdt.quadrature_point(q));

        dpdt.submit_value(flux, q);
      }

      dpdt.integrate_scatter(dealii::EvaluationFlags::values, dst);
    }
  }

  std::function<scalar(qpoint const &)>
  get_scaling_function() const
  {
    // In case we blend in the source term, we check if the scaling is space dependent. Only in that
    // case we have to evaluate the function in every equadrature point. Otherwise the scaling
    // is purely temporal and constant during this function.
    if(data.blend_in)
    {
      AssertThrow(data.blend_in_function != nullptr,
                  dealii::ExcMessage("No blend-in function provided."));
    }

    bool const space_dependent_scaling =
      data.blend_in_function != nullptr ? data.blend_in_function->varies_in_space(time) : false;
    Number const pure_temporal_scaling_factor =
      (not space_dependent_scaling) ? data.blend_in_function->compute_time_factor(time) : 1.0;

    if(space_dependent_scaling)
    {
      return [&](qpoint const & q)
      { return FunctionEvaluator<0, dim, Number>::value(*data.blend_in_function, q, time); };
    }
    else
    {
      // capture scaling factor by copy
      return [pure_temporal_scaling_factor](qpoint const &)
      { return pure_temporal_scaling_factor; };
    }
  }


  dealii::MatrixFree<dim, Number> const * matrix_free;

  SourceTermCalculatorData<dim> data;

  lazy_ptr<VectorType> velocity_cfd;
  lazy_ptr<VectorType> pressure_cfd;

  double time;
};
} // namespace AeroAcoustic
} // namespace ExaDG

#endif /*INCLUDE_EXADG_AERO_ACOUSTIC_SOURCE_TERM_CALCULATOR_H_*/
