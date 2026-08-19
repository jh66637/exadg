#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>


#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/solution_transfer.h>

using namespace dealii;

template <int dim>
class InitialConditionVibratingMembrane : public Function<dim>
{
public:
  InitialConditionVibratingMembrane(const double modes);

  // Function that the gives the initial pressure (comp 0) and velocity (comp
  // 1 to 1 + dim).
  double value(const Point<dim> &p, const unsigned int comp) const final;

  // Function that calculates the duration of one oscillation.
  double get_period_duration(const double speed_of_sound) const;

private:
  const double M;
};

template <int dim>
InitialConditionVibratingMembrane<dim>::InitialConditionVibratingMembrane(
  const double modes)
  : Function<dim>(dim + 1, 0.0)
  , M(modes)
{
  static_assert(dim == 2, "Only implemented for dim==2");
}

template <int dim>
double
InitialConditionVibratingMembrane<dim>::value(const Point<dim>  &p,
                                              const unsigned int comp) const
{
  if (comp == 0)
    {

    // return std::sin(M * numbers::PI * p[0]) * std::sin(M * numbers::PI * p[1]);
  double sigma = 0.01;
  double amplitude = 1.0;
  double dx = p[0] - 0.5;
  double dy = p[1] - 0.5;
  double exponent = -((dx * dx + dy * dy) / (2 * sigma * sigma));
  return amplitude * std::exp(exponent);
    }

  return 0.0;
}

template <int dim>
double InitialConditionVibratingMembrane<dim>::get_period_duration(
  const double speed_of_sound) const
{
  return 2.0 / (M * std::sqrt(dim) * speed_of_sound);
}

namespace HelperFunctions
{
  template <int dim, typename Number, typename VectorType>
  void set_initial_condition(MatrixFree<dim, Number> matrix_free,
                             const Function<dim>    &initial_solution,
                             VectorType             &dst)
  {
    VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                             matrix_free.get_dof_handler(),
                             initial_solution,
                             dst);
  }

  double
  compute_dt_cfl(const double hmin, const unsigned int degree, const double c)
  {
    return hmin / (std::pow(degree, 1.5) * c);
  }

  template <typename VectorType, int dim>
  void write_vtu(const VectorType      &solution,
                 const DoFHandler<dim> &dof_handler,
                 const Mapping<dim>    &mapping,
                 const unsigned int     degree,
                 const std::string     &name_prefix)
  {
    DataOut<dim>          data_out;
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
                             interpretation(dim + 1,
                     DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> names(dim + 1, "U");

    interpretation[0] = DataComponentInterpretation::component_is_scalar;
    names[0]          = "P";

    data_out.add_data_vector(dof_handler, solution, names, interpretation);

    data_out.build_patches(mapping, degree, DataOut<dim>::curved_inner_cells);
    data_out.write_vtu_in_parallel(name_prefix + ".vtu",
                                   dof_handler.get_communicator());
  }
} // namespace HelperFunctions

template <int dim, typename Number>
class BCEvaluationP
{
public:
  BCEvaluationP(const FEFaceEvaluation<dim, -1, 0, 1, Number> &pressure_m)
    : pressure_m(pressure_m)
  {}

  typename FEFaceEvaluation<dim, -1, 0, 1, Number>::value_type
  get_value(const unsigned int q) const
  {
    return -pressure_m.get_value(q);
  }

private:
  const FEFaceEvaluation<dim, -1, 0, 1, Number> &pressure_m;
};

template <int dim, typename Number>
class BCEvaluationU
{
public:
  BCEvaluationU(const FEFaceEvaluation<dim, -1, 0, dim, Number> &velocity_m)
    : velocity_m(velocity_m)
  {}

  typename FEFaceEvaluation<dim, -1, 0, dim, Number>::value_type
  get_value(const unsigned int q) const
  {
    return velocity_m.get_value(q);
  }

private:
  const FEFaceEvaluation<dim, -1, 0, dim, Number> &velocity_m;
};


template <int dim, typename Number>
class AcousticOperator
{
public:
  AcousticOperator(const MatrixFree<dim, Number> &matrix_free)
    : matrix_free(matrix_free)
  {}

  mutable  unsigned int category = -1;

  template <typename VectorType>
  void evaluate_add(VectorType &dst, const VectorType &src, unsigned int category_) const
  {
    category=category_;
    {
      // Perform matrix free loop with point-to-point interpolation at
      // non-matching faces.
      matrix_free.loop(&AcousticOperator::local_apply_cell<VectorType>,
                       &AcousticOperator::local_apply_face<VectorType>,
                       &AcousticOperator::local_apply_boundary_face<VectorType>,
                       this,
                       dst,
                       src,
                       false,
                       MatrixFree<dim, Number>::DataAccessOnFaces::values,
                       MatrixFree<dim, Number>::DataAccessOnFaces::values);
    }
    category=-1;
  }

private:
  template <typename VectorType>
  void local_apply_cell(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, 1, Number>   pressure(matrix_free, 0, 0, 0);
    FEEvaluation<dim, -1, 0, dim, Number> velocity(matrix_free, 0, 0, 1);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        pressure.reinit(cell);

        pressure.gather_evaluate(src, EvaluationFlags::gradients);
        velocity.gather_evaluate(src, EvaluationFlags::gradients);

        for (unsigned int q : pressure.quadrature_point_indices())
          {
            pressure.submit_value(velocity.get_divergence(q), q);
            velocity.submit_value(pressure.get_gradient(q), q);
          }

        pressure.integrate_scatter(EvaluationFlags::values, dst);
        velocity.integrate_scatter(EvaluationFlags::values, dst);
      }
  }

  template <typename VectorType>
  void local_apply_face(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_m(
      matrix_free, true, 0, 0, 0);
    FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_p(
      matrix_free, false, 0, 0, 0);
    FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_m(
      matrix_free, true, 0, 0, 1);
    FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_p(
      matrix_free, false, 0, 0, 1);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        auto const [c_m, c_p] = matrix_free.get_face_category(face);
        if(c_m != category && c_p != category)
          {
            continue;
          }


        velocity_m.reinit(face);
        velocity_p.reinit(face);

        pressure_m.reinit(face);
        pressure_p.reinit(face);

        pressure_m.gather_evaluate(src, EvaluationFlags::values);
        pressure_p.gather_evaluate(src, EvaluationFlags::values);

        velocity_m.gather_evaluate(src, EvaluationFlags::values);
        velocity_p.gather_evaluate(src, EvaluationFlags::values);

        const double tau   = 0.5;
        const double gamma = 0.5;
        for (unsigned int q : pressure_m.quadrature_point_indices())
          {
            const auto n  = pressure_m.normal_vector(q);
            const auto pm = pressure_m.get_value(q);
            const auto um = velocity_m.get_value(q);

            const auto pp = pressure_p.get_value(q);
            const auto up = velocity_p.get_value(q);

            // Compute homogeneous local Lax-Friedrichs fluxes and submit the
            // corrsponding values to the integrators.
            const auto momentum_flux =
              0.5 * (pm + pp) + 0.5 * tau * (um - up) * n;
            velocity_m.submit_value((momentum_flux - pm) * n, q);
            velocity_p.submit_value((momentum_flux - pp) * (-n), q);

            const auto mass_flux =
              0.5 * (um + up) + 0.5 * gamma * (pm - pp) * n;
            pressure_m.submit_value((mass_flux - um) * n, q);
            pressure_p.submit_value((mass_flux - up) * (-n), q);
          }

        if (c_m==category)
          {
            pressure_m.integrate_scatter(EvaluationFlags::values, dst);
            velocity_m.integrate_scatter(EvaluationFlags::values, dst);
          }
        if (c_p==category)
          {
            pressure_p.integrate_scatter(EvaluationFlags::values, dst);
            velocity_p.integrate_scatter(EvaluationFlags::values, dst);
          }
      }
  }

  template <typename VectorType>
  void local_apply_boundary_face(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    // Standard face evaluators.
    FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_m(
      matrix_free, true, 0, 0, 0);
    FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_m(
      matrix_free, true, 0, 0, 1);

    // Classes that return the correct BC values.
    BCEvaluationP pressure_bc(pressure_m);
    BCEvaluationU velocity_bc(velocity_m);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {

        auto const [c_m, _] = matrix_free.get_face_category(face);
        if(c_m != category)
         {
             continue;
         }

        velocity_m.reinit(face);
        pressure_m.reinit(face);

        pressure_m.gather_evaluate(src, EvaluationFlags::values);
        velocity_m.gather_evaluate(src, EvaluationFlags::values);

        // Compute penalty parameters from material parameters.
        const auto tau   = 0.5;
        const auto gamma = 0.5;

        for (unsigned int q : pressure_m.quadrature_point_indices())
          {
            const auto n  = pressure_m.normal_vector(q);
            const auto pm = pressure_m.get_value(q);
            const auto um = velocity_m.get_value(q);

            const auto pp = pressure_bc.get_value(q);
            const auto up = velocity_bc.get_value(q);

            const auto momentum_flux =
              0.5 * (pm + pp) + 0.5 * tau * (um - up) * n;
            velocity_m.submit_value((momentum_flux - pm) * n, q);
            const auto mass_flux =
              0.5 * (um + up) + 0.5 * gamma * (pm - pp) * n;
            pressure_m.submit_value((mass_flux - um) * n, q);
          }

        pressure_m.integrate_scatter(EvaluationFlags::values, dst);
        velocity_m.integrate_scatter(EvaluationFlags::values, dst);
      }
  }

  const MatrixFree<dim, Number> &matrix_free;
};

template <int dim, typename Number>
class InverseMassOperator
{
public:
  InverseMassOperator(const MatrixFree<dim, Number> &matrix_free)
    : matrix_free(matrix_free)
  {}

  mutable unsigned int category=-1;

  // Function to apply the inverse mass operator.
  template <typename VectorType>
  void apply(VectorType &dst, const VectorType &src, unsigned int category_) const
  {
    category = category_;

    dst.zero_out_ghost_values();
    matrix_free.cell_loop(&InverseMassOperator::local_apply_cell<VectorType>,
                          this,
                          dst,
                          src);
    category=-1;
  }

private:
  // Apply the inverse mass operator onto every cell batch.
  template <typename VectorType>
  void local_apply_cell(
    const MatrixFree<dim, Number>               &mf,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, dim + 1, Number> phi(mf);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, dim + 1, Number>
      minv(phi);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        if (matrix_free.get_cell_category(cell)!=category)
          {
            continue;
          }

        phi.reinit(cell);
        phi.read_dof_values(src);
        minv.apply(phi.begin_dof_values(), phi.begin_dof_values());
        phi.set_dof_values(dst);
      }
  }

  const MatrixFree<dim, Number> &matrix_free;
};

// ... (Keep your InitialCondition and HelperFunctions as they are)

template <int dim, typename Number>
class ADEROperator
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

public:
  // dofs are at the standard positions.
  std::vector<double>
  compute_subinterval_weights(
      const Quadrature<1>                &quad,
                              const double a,
                              const double b)const
  {
    assert(a<=b);
    assert(a>=0 -1e-12);
    assert(b<=1 +1e-12);

    const unsigned int n_q_points_1d = quad.size();
    const std::vector<Polynomials::Polynomial<double>> poly_coll =
      Polynomials::generate_complete_Lagrange_basis(quad.get_points());

    QGauss<1> internal_quad(n_q_points_1d + 2);
    double const J = (b - a);

    std::vector<double> sub_weights(n_q_points_1d, 0.0);

    for (unsigned int i = 0; i < n_q_points_1d; ++i)
      {
        for (unsigned int q = 0; q < internal_quad.size(); ++q)
          {
            sub_weights[i] += poly_coll[i].value(a + internal_quad.point(q)[0] * J)*internal_quad.weight(q) * J;
          }
      }
    return sub_weights;
  }


  ADEROperator(const MatrixFree<dim, Number> &mf)
    : matrix_free(mf)
  {
    const unsigned int degree = mf.get_dof_handler().get_fe().degree;
    const unsigned int n_dofs = degree + 1;

    FE_Q<1>          fe_time(degree);
    QGaussLobatto<1> q_int(n_dofs);// lagrange dof positions

    for (unsigned int p = 0; p < n_dofs; ++p) // shape functions
      {
        for (unsigned int l = 0; l < n_dofs; ++l) // time nodes
          {
            if (l == 0)
              {
                I_matrix[l][p] = 0.0; // Integral from 0 to 0 is always 0
                continue;
              }

            // We need to integrate phi_p from 0 to tau_l.
            // We create a "sub-quadrature" scaled to the interval [0, tau_l]
            double tau_l    = q_int.point(l)[0];
            double integral = 0.0;

            for (unsigned int q = 0; q < q_int.size(); ++q)
              {
                // Map Gauss point from [0, 1] to [0, tau_l]
                double weight_mapped = q_int.get_weights()[q] * tau_l;
                double phi_p = fe_time.shape_value(p, q_int.point(q) * tau_l);
                integral += phi_p * weight_mapped;
              }
            I_matrix[l][p] = integral;
          }
      }
  }

  mutable double dt = -100;
  mutable unsigned int category = -1;
 mutable std::vector<VectorType> *st_predictor;
  // The Picard Iteration: Computes the space-time trajectory inside each cell
  void compute_predictor(std::vector<VectorType> &st_predictor_,
                         const VectorType        &u_n,
                         unsigned int category_,
                         const double             dt_
                         ) const
  {
    dt = dt_;
    category=category_;
    st_predictor=&st_predictor_;

    VectorType a=u_n;

    matrix_free.cell_loop(&ADEROperator::cell_local_picard_iteration,
                          this,
                          a,
                          u_n,false);

    category=-1;
    dt=-100;
  }


  void space_time_average(VectorType &dst,
           const std::vector<VectorType> &st_predictor,
           const double sub_time,//[0,dt_max]
           std::vector<std::pair<int, double>> all) const
  {
    double dt_max = all.front().second;
    double dt_min = all.back().second;
    dst = 0.0;

    for (auto [category_, dt_]:all)
      {


    category =category_;
    dt = dt_;

        double sub_interval_begin = sub_time/dt + 1e-6;
        sub_interval_begin -= std::floor(sub_interval_begin) +1e-6;

       double sub_interval_end = sub_interval_begin + dt_min/dt;
        // std::cerr<<category_<<std::endl;
        // std::cerr<<sub_interval_begin<<" ";
        // std::cerr<<sub_interval_end<<std::endl;
         assert(sub_interval_end - 1e-6 <= 1.0);

    QGaussLobatto<1>quad(st_predictor.size());//aka degree+1
    auto quad_sub= compute_subinterval_weights(quad, sub_interval_begin,sub_interval_end);
    for (auto&w:quad_sub)
      {
        w*=dt;
      }

        matrix_free.cell_loop(&ADEROperator::average,
                      this,
                      dst,
                      std::make_pair(quad_sub,st_predictor),
                      false);


    category =-1;
    dt=-100;
      }
  }

  void average(
  const MatrixFree<dim, Number> &mf,
  VectorType &dst,
  const std::pair<std::vector<double>,std::vector<VectorType>> &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, dim+1, Number>  phi1(mf, 0, 0, 0);
    FEEvaluation<dim, -1, 0, dim+1, Number>  phi2(mf, 0, 0, 0);

    for (unsigned int cell = cell_range.first; cell < cell_range.second;++cell)
      {
        if (matrix_free.get_cell_category(cell)!=category)
          {
            continue;
          }

        phi1.reinit(cell);
        phi1.read_dof_values(dst);

        for (unsigned int i=0;i<src.first.size();++i)
          {
            phi2.reinit(cell);
            phi2.read_dof_values(src.second[i]);
            auto const weight=src.first[i];
            for (unsigned int j =0; j< phi1.dofs_per_component; ++j)
              {
                phi1.submit_dof_value(phi1.get_dof_value(j) + weight*phi2.get_dof_value(j),j);
              }
          }

        phi1.set_dof_values(dst);
      }
  }

  void add(
  const MatrixFree<dim, Number> &mf,
  VectorType &dst,
  const std::pair<double,VectorType> &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, dim+1, Number>  phi1(mf, 0, 0, 0);
    FEEvaluation<dim, -1, 0, dim+1, Number>  phi2(mf, 0, 0, 0);

    for (unsigned int cell = cell_range.first; cell < cell_range.second;++cell)
      {
        if (matrix_free.get_cell_category(cell)!=category)
          {
            continue;
          }

        phi1.reinit(cell);
        phi1.read_dof_values(dst);

        phi2.reinit(cell);
        phi2.read_dof_values(src.second);

        for (unsigned int j =0; j< phi1.dofs_per_component; ++j)
          {
            phi1.submit_dof_value(phi1.get_dof_value(j) + src.first*phi2.get_dof_value(j),j);
          }

        phi1.set_dof_values(dst);
      }
  }



  // The Corrector: Uses the predictor to update the solution via fluxes
  void correct(VectorType                    &dst,
             VectorType const &predictor_avg,
             unsigned int category_,
             double const dt_
             ) const
  {
    dt=dt_;
    category=category_;

    InverseMassOperator inverse_mass_operator(matrix_free);
    AcousticOperator<dim, Number> acoustic_operator(matrix_free);

    VectorType rhs = predictor_avg;
    rhs=0.0;
    acoustic_operator.evaluate_add(rhs, predictor_avg, category);
    inverse_mass_operator.apply(rhs, rhs, category);

    matrix_free.cell_loop(&ADEROperator::add,
              this,
              dst,
              std::make_pair(-1.0,rhs),
              false);


    category=-1;
    dt=-100.0;
  }

private:
    void cell_local_picard_iteration(
    const MatrixFree<dim, Number> &mf,
    VectorType &,
    const VectorType                            &u_n,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    const unsigned int degree = 3;

    FEEvaluation<dim, -1, 0, 1, Number>    p_eval(mf, 0, 0, 0);
    FEEvaluation<dim, -1, 0, dim, Number>  u_eval(mf, 0, 0, 1);

    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, 1, Number> minv_p(p_eval);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, dim, Number> minv_u(u_eval);

    AlignedVector<VectorizedArray<Number>> p_0(p_eval.dofs_per_cell);
    AlignedVector<VectorizedArray<Number>> u_0(u_eval.dofs_per_cell);

    std::vector<AlignedVector<VectorizedArray<Number>>> rhs_p(st_predictor->size(), p_0);
    std::vector<AlignedVector<VectorizedArray<Number>>> rhs_u(st_predictor->size(), u_0);

    // for every cell compute the cell local space time predictor using picard
    // iteration
    for (unsigned int cell = cell_range.first; cell < cell_range.second;++cell)
      {
        if (matrix_free.get_cell_category(cell)!=category)
          {
            continue;
          }

        // 1. Get initial state and fill cell local space time predictor buffers
        p_eval.reinit(cell);
        p_eval.read_dof_values(u_n);
        std::copy_n(p_eval.begin_dof_values(),p_eval.dofs_per_cell,p_0.begin());
        std::fill(rhs_p.begin(),rhs_p.end(),p_0);

        u_eval.reinit(cell);
        u_eval.read_dof_values(u_n);
        std::copy_n(u_eval.begin_dof_values(),u_eval.dofs_per_cell,u_0.begin());
        std::fill(rhs_u.begin(),rhs_u.end(),u_0);

        // 2. PICARD ITERATIONS
        for (unsigned int iter = 0; iter < degree; ++iter)
          {
            // 2. Compute RHS for each time node
            // A. For each time node, compute the spatial RHS
            for (unsigned int t = 0; t < st_predictor->size(); ++t)
              {
                // compute spatial integral
                std::copy(rhs_p[t].begin(),rhs_p[t].end(),p_eval.begin_dof_values());
                p_eval.evaluate(EvaluationFlags::gradients);

                std::copy(rhs_u[t].begin(),rhs_u[t].end(),u_eval.begin_dof_values());
                u_eval.evaluate(EvaluationFlags::gradients);

                for (unsigned int q : p_eval.quadrature_point_indices())
                  {
                    p_eval.submit_value(-u_eval.get_divergence(q), q);
                    u_eval.submit_value(-p_eval.get_gradient(q), q);
                  }

                p_eval.integrate(EvaluationFlags::values);
                u_eval.integrate(EvaluationFlags::values);

                minv_p.apply(p_eval.begin_dof_values(),rhs_p[t].begin());
                minv_u.apply(u_eval.begin_dof_values(), rhs_u[t].begin());
              }

            // 3. Update Predictor nodes: Q = U_n + dt * I * RHS
            auto p_node=rhs_p;
            auto u_node=rhs_u;
            for (unsigned int l = 0; l < st_predictor->size(); ++l)
              {
                p_node[l] = p_0;
                u_node[l] = u_0;

                for (unsigned int p = 0; p < st_predictor->size(); ++p)
                  {
                    for (unsigned int i = 0; i < p_eval.dofs_per_cell; ++i)
                      {
                        p_node[l][i] += dt * I_matrix[l][p] * rhs_p[p][i];
                        for (unsigned int d = 0; d < dim; ++d)
                          {
                            u_node[l][d * p_eval.dofs_per_cell + i] +=dt * I_matrix[l][p] *rhs_u[p][d * p_eval.dofs_per_cell + i];
                          }
                      }
                  }
              }
            rhs_p=p_node;
            rhs_u=u_node;
          }

        for (unsigned int l = 0; l < st_predictor->size(); ++l)
          {
            std::copy(rhs_p[l].begin(),rhs_p[l].end(),p_eval.begin_dof_values());
            p_eval.set_dof_values(st_predictor->operator[](l));
            std::copy(rhs_u[l].begin(),rhs_u[l].end(),u_eval.begin_dof_values());
            u_eval.set_dof_values(st_predictor->operator[](l));
          }
      }
  }


  const MatrixFree<dim, Number> &matrix_free;
  double                         I_matrix[4][4];
};

template<int dim>
struct MyErrorEstimator
{
  template <typename Number>
static void estimate(
const DoFHandler<dim> &dof,
const Quadrature<dim>       &quad,
const ReadVector<Number> &solution,
Vector<float>            &error)
{
  FEValues<dim> fe_values(dof.get_fe(),quad,update_gradients|update_JxW_values);
    using grad_t = Tensor<1, dim, Number>;
  std::vector<std::vector<grad_t>> buffer{quad.size(),std::vector<grad_t>{dof.get_fe().n_components()}};//[q][comp]
    for (auto &cell : dof.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values.get_function_gradients(solution,buffer);
            double err=0;
            for (unsigned int c = 0; c < dof.get_fe().n_components(); ++c)
              {
             for (unsigned int q = 0; q < quad.size(); ++q)
              {
                err+=buffer[q][c].norm_square()*fe_values.JxW(q);
              }
              }
           error[cell->active_cell_index()] = std::sqrt(err);
          }
      }


}
};

// --- Updated Time Stepper ---
template <int dim, typename Number>
class ADERTimeStepper
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

public:


  void do_sub_stepping(ADEROperator<dim, Number>&ader,VectorType& solution, std::vector<VectorType>& st_predictor, VectorType& accumulated, std::vector<std::pair<int,double>> all, int index, double& sub_time)
  {

    if (index>=all.size())
      {
        return;
      }

    VectorType avg = accumulated;

    double this_dt = all[index].second;
    int id = all[index].first;
    assert(index > 0);
    double parent_dt = all[index-1].second;

    int n_dt = std::floor((parent_dt/this_dt)+1e-6);
    double sub_dt = parent_dt / n_dt;

    for (int i =0; i <n_dt; ++i)
    {
      ader.compute_predictor(st_predictor, solution, id, sub_dt);

        avg=0.0;
      do_sub_stepping(ader, solution, st_predictor, avg, all, index+1, sub_time);

      if (index == all.size()-1)
      {
          ader.space_time_average(avg, st_predictor, sub_time, all);
          sub_time+=sub_dt;
      }

      ader.correct(solution, avg, id, sub_dt);
      accumulated.sadd(1.0,1.0,avg);
    }
  }

  MatrixFree<dim, Number>* matrix_free__=nullptr;
  Triangulation<dim>*tria__=nullptr;
  DoFHandler<dim>*dof__=nullptr;

  //TODO: we could as well have p refinement here, see step 75
void adapt_resolution(VectorType& solution, std::vector<VectorType>& st_predictor, VectorType& avg)
{
  auto& triangulation = *tria__;
  auto& dof_handler = *dof__;

  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  QGauss<dim>quadrature(dof_handler.get_fe().degree + 1);
  MyErrorEstimator<dim>::estimate(dof_handler,
                                     quadrature,
                                     solution,
                                     estimated_error_per_cell);


  parallel::distributed::SolutionTransfer<dim, VectorType> soltrans(dof_handler);

  parallel::distributed::GridRefinement::
 refine_and_coarsen_fixed_fraction(triangulation,
                                   estimated_error_per_cell,
                                   0.6, // refine the cells that in total make up 10% of the error
                                   0.1 );//coarsen the cells that in total make up 10% of the error

  triangulation.prepare_coarsening_and_refinement();

  soltrans.prepare_for_coarsening_and_refinement(solution);


  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->level() >= 5)
            {
              cell->clear_refine_flag();
            }
        }
    }

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(dof_handler.get_fe());

  auto locally_owned_dofs = dof_handler.locally_owned_dofs();
  auto locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, dof_handler.get_communicator());

  // transfer solution
  soltrans.interpolate(solution);

  // the rest of the vectors are zeroed out anyway
  for (auto& st:st_predictor)
    st.reinit(locally_owned_dofs, locally_relevant_dofs, dof_handler.get_communicator());

  avg.reinit(locally_owned_dofs, locally_relevant_dofs, dof_handler.get_communicator());
}

  void reinit_mf(std::vector<unsigned int>const& cell_categories)
{
  auto& triangulation = *tria__;
  auto& dof_handler = *dof__;


    //reinit matrixfree
    typename MatrixFree<dim, Number>::AdditionalData data;
    data.mapping_update_flags = update_gradients | update_values | update_quadrature_points;
    data.mapping_update_flags_inner_faces = update_values;
    data.mapping_update_flags_boundary_faces = update_quadrature_points | update_values;

  data.cell_vectorization_category = cell_categories;

    data.cell_vectorization_categories_strict = true;
    data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::none;

  AffineConstraints<Number> constraints;
  constraints.close();

    matrix_free__->reinit(
      MappingQ1<dim>(), dof_handler, constraints, QGauss<dim>(dof_handler.get_fe().degree + 1), data);
  }

std::pair<std::vector<unsigned int>,std::vector<std::pair<int,double>>> compute_cell_categories_and_dt(unsigned int degree, double cr)
{

  auto& triangulation = *tria__;
  auto& dof_handler = *dof__;

    std::vector<double> diameters(triangulation.n_active_cells());

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      diameters[cell->active_cell_index()] = cell->diameter();
    }
    std::sort(diameters.begin(), diameters.end());
    const double h_min = Utilities::MPI::min(diameters.front(), dof_handler.get_communicator());
    const double dt_min = cr * HelperFunctions::compute_dt_cfl(h_min, degree, 1.0);

  std::vector<unsigned int> cell_categories(triangulation.n_active_cells());
  std::set<unsigned int> unique_categories;
  for(const auto & cell : triangulation.active_cell_iterators())
    {
      if(cell->is_locally_owned())
        {
          cell_categories[cell->active_cell_index()] = static_cast<unsigned int>((cell->diameter()/h_min)+1e-3);
          unique_categories.insert(cell_categories[cell->active_cell_index()]);
        }
    }


  // create pairs of categories and dt. has to be sorted from large dt to small dt
  std::vector<std::pair<int,double>> categories_dt;
  for (auto c: unique_categories)
  {
      categories_dt.emplace_back(c, dt_min/c);
  }
  std::sort(categories_dt.begin(), categories_dt.end(),[](auto a, auto b){return b.second<a.second;});

  return {cell_categories, categories_dt};
}

  void run( MatrixFree<dim, Number> &matrix_free,
  Triangulation<dim>&tria,
  DoFHandler<dim>&dh,
           const double                   cr,
           const double                   end_time,
           const Function<dim>           &initial_condition)
  {
    matrix_free__ = &matrix_free;
    tria__ = &tria;
  dof__=&dh;
    const auto &dof_handler = matrix_free.get_dof_handler();
    const auto  degree      = dof_handler.get_fe().degree;


    VectorType solution;
    matrix_free.initialize_dof_vector(solution);
    std::vector<VectorType> st_predictor(degree+1);
    for (auto &v : st_predictor)
      matrix_free.initialize_dof_vector(v);

    HelperFunctions::set_initial_condition(matrix_free, initial_condition,solution);
    ADEROperator<dim, Number> ader(matrix_free);

    double h_local_max = std::numeric_limits<double>::lowest();
    for (const auto &cell : dof_handler.active_cell_iterators())
      h_local_max = std::max(h_local_max,
                             (cell->vertex(1) - cell->vertex(0)).norm_square());
    h_local_max = std::sqrt(h_local_max);
    const double h_max =
      Utilities::MPI::min(h_local_max, dof_handler.get_communicator());

    const double dt = cr * HelperFunctions::compute_dt_cfl(h_max, degree, 1.0);

  VectorType avg;
  matrix_free.initialize_dof_vector(avg);

    double       time     = 0;
    unsigned int timestep = 0;
    while (time < end_time)
      {
        if (timestep%10==0)
          {
        HelperFunctions::write_vtu(solution,
                                   matrix_free.get_dof_handler(),
                                   *matrix_free.get_mapping_info().mapping,
                                   degree,
                                   "step_89-ader" + std::to_string(timestep));
          }
        timestep++;

        // 1. The Picard Predictor (Cell-local, 0 Communication)
        //TODO: p refinement?
        //TODO: not every timestep!
        adapt_resolution(solution, st_predictor,avg);
        auto const&[cell_categories, categories_dt] = compute_cell_categories_and_dt(degree,cr);
        reinit_mf(cell_categories);

        avg=0.0;
        std::for_each(st_predictor.begin(),st_predictor.end(),[&](auto&st){st=0.0;});

        auto index = 0;

        ader.compute_predictor(st_predictor, solution, categories_dt[index].first, categories_dt[index].second);
        double sub_time =0.0;
        do_sub_stepping(ader,solution, st_predictor,avg, categories_dt, index+1, sub_time);
        ader.correct(solution, avg, categories_dt[index].first, categories_dt[index].second);

        time += dt;
      }
  }
};


template <int dim, typename Number>
void run( MatrixFree<dim, Number> &matrix_free,
  Triangulation<dim>&tria,
  DoFHandler<dim>&dh,
         const double                   end_time,
         const Function<dim>           &initial_condition)
{
  ADERTimeStepper<dim, Number> ader;
  ader.run(matrix_free,
    tria,
    dh,
           /*Cr*/ 0.1,
           end_time,
           initial_condition);
}

int main(int argc, char *argv[])
{
  constexpr int dim = 2;
  using Number      = double;

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  std::cout.precision(5);
  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  const unsigned int refinements = 3;
  const unsigned int degree      = 3;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  pcout << "Create grid..." << std::endl;

  GridGenerator::subdivided_hyper_rectangle(tria,
                                            {1, 1},
                                            {0.0, 0.0},
                                            {1.0, 1.0});
  tria.refine_global(refinements);

  // // Refine mesh adaptively once again in region around sphere and in the
  // // immediate wake
  // for(auto const & cell : tria.active_cell_iterators())
  //   if(cell->is_locally_owned())
  //     {
  //       if(cell->center()[0] > 0.5)
  //         {
  //             cell->set_refine_flag();
  //         }
  //     }
  // tria.execute_coarsening_and_refinement();
  //
  //
  // for(auto const & cell : tria.active_cell_iterators())
  //   if(cell->is_locally_owned())
  //     {
  //        if (cell->center()[0]>0.7&&cell->center()[0]<0.9&&cell->center()[1]>0.6&&cell->center()[1]<0.8)
  //         {
  //           cell->set_refine_flag();
  //         }
  //     }
  // tria.execute_coarsening_and_refinement();
  //
  //
  // for(auto const & cell : tria.active_cell_iterators())
  //   if(cell->is_locally_owned())
  //     {
  //       if (cell->center()[0]>0.8&&cell->center()[0]<1.0&&cell->center()[1]>0.4&&cell->center()[1]<0.7)
  //         {
  //           cell->set_refine_flag();
  //         }
  //     }
  // tria.execute_coarsening_and_refinement();



  pcout << " - Refinement level: " << refinements << std::endl;
  pcout << " - Number of cells: " << tria.n_cells() << std::endl;

  // Set up MatrixFree.

  pcout << "Create DoFHandler..." << std::endl;
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_DGQ<dim>(degree) ^ (dim + 1)));
  pcout << " - Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  AffineConstraints<Number> constraints;
  constraints.close();

  pcout << "Set up MatrixFree..." << std::endl;
  typename MatrixFree<dim, Number>::AdditionalData data;
  data.mapping_update_flags =
    update_gradients | update_values | update_quadrature_points;
  data.mapping_update_flags_inner_faces = update_values;
  data.mapping_update_flags_boundary_faces =
    update_quadrature_points | update_values;

      auto &       cell_categories = data.cell_vectorization_category;
      cell_categories.resize(tria.n_active_cells());
      for(const auto & cell : tria.active_cell_iterators())
        {
          if(cell->is_locally_owned())
            {
              // double d = cell->diameter();
              // if (d>0.17)
              //   {
                  cell_categories[cell->active_cell_index()] = 0;
             //    }
             // else if (d>0.08)
             //    {
             //      cell_categories[cell->active_cell_index()] = 1;
             //    }
             // else if (d>0.04)
             //   {
             //     cell_categories[cell->active_cell_index()] = 2;
             //   }
             // else if (d>0.02)
             //   {
             //     cell_categories[cell->active_cell_index()] = 3;
             //   }
             //
            }
        }
      data.cell_vectorization_categories_strict = true;
  data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::none;

  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(
    MappingQ1<dim>(), dof_handler, constraints, QGauss<dim>(degree + 1), data);

  const auto initial_solution_membrane =
    InitialConditionVibratingMembrane<dim>(2);

  pcout << "Run Simulation..." << std::endl;

  run(matrix_free, tria,dof_handler,1.0, initial_solution_membrane);

  return 0;
}
