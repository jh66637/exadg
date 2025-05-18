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

// deal.II
#include <deal.II/numerics/vector_tools.h>

// ExaDG
#include <exadg/acoustic_conservation_equations/spatial_discretization/spatial_operator.h>
#include <exadg/grid/mapping_dof_vector.h>
#include <exadg/operators/finite_element.h>
#include <exadg/operators/grid_related_time_step_restrictions.h>
#include <exadg/operators/quadrature.h>
#include <exadg/utilities/exceptions.h>

namespace ExaDG
{
namespace Acoustics
{
template<int dim, typename Number>
SpatialOperator<dim, Number>::SpatialOperator(
  std::shared_ptr<Grid<dim> const>               grid_in,
  std::shared_ptr<dealii::Mapping<dim> const>    mapping_in,
  std::shared_ptr<BoundaryDescriptor<dim> const> boundary_descriptor_in,
  std::shared_ptr<FieldFunctions<dim> const>     field_functions_in,
  Parameters const &                             parameters_in,
  std::string const &                            field_in,
  MPI_Comm const &                               mpi_comm_in)
  : Interface::SpatialOperator<Number>(),
    grid(grid_in),
    mapping(mapping_in),
    boundary_descriptor(boundary_descriptor_in),
    field_functions(field_functions_in),
    param(parameters_in),
    field(field_in),
    dof_handler_p(*grid_in->triangulation),
    dof_handler_u(*grid_in->triangulation),
    aero_acoustic_source_term(nullptr),
    mpi_comm(mpi_comm_in),
    pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_comm_in) == 0)
{
  pcout << std::endl
        << "Construct acoustic conservation equations operator ..." << std::endl
        << std::flush;

  // TODO:
  //  if(param.has_pml)
  //  {
  //    // TODO: quick fix. currently categorize_pml_cells is called twice to ensure n_pml_cells is
  //    // always correct
  //    std::vector<unsigned int> temp;
  //    n_pml_cells = PML::Utilities::categorize_pml_cells(dof_handler_p.get_triangulation(), temp);
  //    AssertThrow(n_pml_cells != 0, dealii::ExcMessage("Could not find a PML"));
  //  }

  initialize_dof_handler_and_constraints();

  pcout << std::endl << "... done!" << std::endl << std::flush;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::fill_matrix_free_data(
  MatrixFreeData<dim, Number> & matrix_free_data) const
{
  // append mapping flags
  matrix_free_data.append_mapping_flags(Operators::Kernel<dim, Number>::get_mapping_flags());

  // TODO:
  //  if(param.has_pml)
  //    matrix_free_data.append_mapping_flags(Operators::PMLKernel<dim,
  //    Number>::get_mapping_flags());

  if(param.right_hand_side)
    matrix_free_data.append_mapping_flags(
      ExaDG::Operators::RHSKernel<dim, Number>::get_mapping_flags());

  // mapping flags required for CFL condition
  if(param.calculation_of_time_step_size == TimeStepCalculation::CFL)
  {
    MappingFlags flags_cfl;
    flags_cfl.cells = dealii::update_quadrature_points;
    matrix_free_data.append_mapping_flags(flags_cfl);
  }

  // dof handler
  matrix_free_data.insert_dof_handler(&dof_handler_p, field + dof_index_p);
  matrix_free_data.insert_dof_handler(&dof_handler_u, field + dof_index_u);

  // constraint
  matrix_free_data.insert_constraint(&constraint_p, field + dof_index_p);
  matrix_free_data.insert_constraint(&constraint_u, field + dof_index_u);

  // quadrature for pressure
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_p =
    create_quadrature<dim>(param.grid.element_type, param.degree_p + 1);
  matrix_free_data.insert_quadrature(*quadrature_p, field + quad_index_p);

  // quadrature for velocity and PML auxiliary variable
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_u =
    create_quadrature<dim>(param.grid.element_type, param.degree_u + 1);
  matrix_free_data.insert_quadrature(*quadrature_u, field + quad_index_u);

  // quadrature that works for pressure and velocity
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_p_u =
    create_quadrature<dim>(param.grid.element_type, std::max(param.degree_p, param.degree_u) + 1);
  matrix_free_data.insert_quadrature(*quadrature_p_u, field + quad_index_p_u);

  // TODO:
  //  if(param.has_pml)
  //  {
  //    // divide into pml cells and pure acoustic cells to be able to evaluate
  //    // pml only in a subset of cells
  //    PML::Utilities::categorize_pml_cells(dof_handler_p.get_triangulation(),
  //                                         matrix_free_data.data.cell_vectorization_category);
  //
  //    // TODO: probably we don't need the strict categorization because for mixed batches
  //    // matrix free returns the maximum category which is always the pml category.
  //    // Not using strict categories will make us compute the pml equation in a few unnecessary
  //    cells
  //    // but this will not have any implications on the result. Not using strict categories
  //    enables
  //    // matrix free to run faster, so this is probably the way to go.
  //    // Note: For local time stepping this means that smaller cells need larger numbers since we
  //    // are allowed to perform smaller timesteps on larger cells.
  //    // @Kraxi: can you test if you get the same results with true and false and if you see
  //    // differences in the runtime?
  //    matrix_free_data.data.cell_vectorization_categories_strict = true;
  //  }


  // if(lts)
  //  {
  auto const & tria            = dof_handler_p.get_triangulation();
  auto &       cell_categories = matrix_free_data.data.cell_vectorization_category;
  cell_categories.resize(tria.n_active_cells());

  for(const auto & cell : tria.active_cell_iterators())
  {
    if(cell->is_locally_owned())
    {
      AssertIndexRange(cell->active_cell_index(), tria.n_active_cells());

      auto const p = cell->center();
      if(p[0] > 0.375 && p[0] < 0.625 && p[1] > 0.375 && p[1] < 0.625)
      {
        cell_categories[cell->active_cell_index()] = 1;
      }
      else if(p[0] > 0.25 && p[0] < 0.75 && p[1] > 0.25 && p[1] < 0.75)
      {
        cell_categories[cell->active_cell_index()] = 2;
      }
      else
      {
        cell_categories[cell->active_cell_index()] = 3;
      }
    }
  }

  matrix_free_data.data.cell_vectorization_categories_strict = true;
  // }
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::setup()
{
  // initialize MatrixFree and MatrixFreeData
  std::shared_ptr<dealii::MatrixFree<dim, Number>> mf =
    std::make_shared<dealii::MatrixFree<dim, Number>>();
  std::shared_ptr<MatrixFreeData<dim, Number>> mf_data =
    std::make_shared<MatrixFreeData<dim, Number>>();

  fill_matrix_free_data(*mf_data);

  mf->reinit(*get_mapping(),
             mf_data->get_dof_handler_vector(),
             mf_data->get_constraint_vector(),
             mf_data->get_quadrature_vector(),
             mf_data->data);

  // Subsequently, call the other setup function with MatrixFree/MatrixFreeData objects as
  // arguments.
  this->setup(mf, mf_data);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::setup(
  std::shared_ptr<dealii::MatrixFree<dim, Number> const> matrix_free_in,
  std::shared_ptr<MatrixFreeData<dim, Number> const>     matrix_free_data_in)
{
  pcout << std::endl
        << "Setup acoustic conservation equations operator ..." << std::endl
        << std::flush;

  // MatrixFree
  matrix_free      = matrix_free_in;
  matrix_free_data = matrix_free_data_in;

  initialize_operators();

  pcout << std::endl << "... done!" << std::endl << std::flush;
}

template<int dim, typename Number>
dealii::MatrixFree<dim, Number> const &
SpatialOperator<dim, Number>::get_matrix_free() const
{
  return *matrix_free;
}

template<int dim, typename Number>
std::string
SpatialOperator<dim, Number>::get_dof_name_pressure() const
{
  return field + dof_index_p;
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_dof_index_pressure() const
{
  return matrix_free_data->get_dof_index(get_dof_name_pressure());
}

template<int dim, typename Number>
std::string
SpatialOperator<dim, Number>::get_dof_name_velocity() const
{
  return field + dof_index_u;
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_dof_index_velocity() const
{
  return matrix_free_data->get_dof_index(get_dof_name_velocity());
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_pressure_velocity() const
{
  return matrix_free_data->get_quad_index(field + quad_index_p_u);
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_pressure() const
{
  return matrix_free_data->get_quad_index(field + quad_index_p);
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_velocity() const
{
  return matrix_free_data->get_quad_index(field + quad_index_u);
}

template<int dim, typename Number>
std::shared_ptr<dealii::Mapping<dim> const>
SpatialOperator<dim, Number>::get_mapping() const
{
  return mapping;
}

template<int dim, typename Number>
dealii::FiniteElement<dim> const &
SpatialOperator<dim, Number>::get_fe_p() const
{
  return *fe_p;
}

template<int dim, typename Number>
dealii::FiniteElement<dim> const &
SpatialOperator<dim, Number>::get_fe_u() const
{
  return *fe_u;
}

template<int dim, typename Number>
dealii::DoFHandler<dim> const &
SpatialOperator<dim, Number>::get_dof_handler_p() const
{
  return dof_handler_p;
}

template<int dim, typename Number>
dealii::DoFHandler<dim> const &
SpatialOperator<dim, Number>::get_dof_handler_u() const
{
  return dof_handler_u;
}

template<int dim, typename Number>
dealii::AffineConstraints<Number> const &
SpatialOperator<dim, Number>::get_constraint_p() const
{
  return constraint_p;
}

template<int dim, typename Number>
dealii::AffineConstraints<Number> const &
SpatialOperator<dim, Number>::get_constraint_u() const
{
  return constraint_u;
}

template<int dim, typename Number>
dealii::types::global_dof_index
SpatialOperator<dim, Number>::get_number_of_dofs() const
{
  return dof_handler_u.n_dofs() + dof_handler_p.n_dofs() + n_pml_cells * fe_u->n_dofs_per_cell();
}

/*
 * Initialization of vectors.
 */
template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_dof_vector(BlockVectorType & dst) const
{
  // TODO:
  // dst.reinit(param.has_pml ? 3 : 2);
  dst.reinit(2);

  matrix_free->initialize_dof_vector(dst.block(block_index_pressure), get_dof_index_pressure());
  matrix_free->initialize_dof_vector(dst.block(block_index_velocity), get_dof_index_velocity());
  // TODO:
  //  if(param.has_pml)
  //    matrix_free->initialize_dof_vector(dst.block(block_index_pml_aux),
  //    get_dof_index_velocity());

  dst.collect_sizes();
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_dof_vector_pressure(VectorType & dst) const
{
  matrix_free->initialize_dof_vector(dst, get_dof_index_pressure());
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::prescribe_initial_conditions(BlockVectorType & dst,
                                                           double const      time) const
{
  field_functions->initial_solution_pressure->set_time(time);
  field_functions->initial_solution_velocity->set_time(time);

  // This is necessary if Number == float
  using VectorTypeDouble = dealii::LinearAlgebra::distributed::Vector<double>;

  VectorTypeDouble pressure_double;
  VectorTypeDouble velocity_double;
  pressure_double = dst.block(block_index_pressure);
  velocity_double = dst.block(block_index_velocity);

  dealii::VectorTools::interpolate(*get_mapping(),
                                   dof_handler_p,
                                   *(field_functions->initial_solution_pressure),
                                   pressure_double);

  dealii::VectorTools::interpolate(*get_mapping(),
                                   dof_handler_u,
                                   *(field_functions->initial_solution_velocity),
                                   velocity_double);

  dst.block(block_index_pressure) = pressure_double;
  dst.block(block_index_velocity) = velocity_double;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::set_aero_acoustic_source_term(
  VectorType const & aero_acoustic_source_term_in)
{
  aero_acoustic_source_term = &aero_acoustic_source_term_in;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::evaluate(BlockVectorType &          dst,
                                       BlockVectorType const &    src,
                                       double const               time,
                                       dealii::types::material_id cell_category) const
{
  evaluate_acoustic_operator(dst, src, time, cell_category);

  // TODO:
  //  if(param.has_pml)
  //  {
  //    // add contributions to mass and momentum equation, and reset pml equation
  //    dst.block(block_index_pml_aux) = 0.0;
  //    pml_operator.evaluate_add(dst, src);
  //  }

  // shift to the right-hand side of the equation
  // TODO: we have do this here and remove -1.0 from apply_scaled_inverse_mass_operator
  // dst *= -1.0;

  // TODO:
  //  if(param.right_hand_side)
  //    rhs_operator.evaluate_add(dst.block(block_index_pressure), time);

  // TODO:
  //  if(param.aero_acoustic_source_term)
  //  {
  //    AssertThrow(aero_acoustic_source_term,
  //                dealii::ExcMessage("Aero-acoustic source term not valid."));
  //    dst.block(block_index_pressure) += *aero_acoustic_source_term;
  //  }

  apply_scaled_inverse_mass_operator(dst, dst, cell_category);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::add_dofs_of_cell_category(
  BlockVectorType &          dst,
  Number                     factor,
  BlockVectorType const &    src,
  dealii::types::material_id cell_category) const
{
  acoustic_operator.add_vectors(dst, factor, src, cell_category);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::copy_dofs_of_cell_category(
  BlockVectorType &          dst,
  BlockVectorType const &    src,
  dealii::types::material_id cell_category) const
{
  acoustic_operator.copy_dofs_of_cell_category(dst, src, cell_category);
}


template<int dim, typename Number>
void
SpatialOperator<dim, Number>::evaluate_acoustic_operator(
  BlockVectorType &          dst,
  BlockVectorType const &    src,
  double const               time,
  dealii::types::material_id cell_category) const
{
  acoustic_operator.evaluate(dst, src, time, cell_category);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::apply_scaled_inverse_mass_operator(
  BlockVectorType &          dst,
  BlockVectorType const &    src,
  dealii::types::material_id cell_category) const
{
  inverse_mass_pressure.apply_scale(dst.block(block_index_pressure),
                                    -1.0 * param.speed_of_sound * param.speed_of_sound,
                                    src.block(block_index_pressure),
                                    cell_category);
  inverse_mass_velocity.apply_scale(dst.block(block_index_velocity),
                                    -1.0,
                                    src.block(block_index_velocity),
                                    cell_category);


  // inverse_mass_pressure.apply_scale(dst.block(block_index_pressure),
  //                                   param.speed_of_sound * param.speed_of_sound,
  //                                   src.block(block_index_pressure),
  //                                   cell_category);
  // inverse_mass_velocity.apply(dst.block(block_index_velocity),
  //                             src.block(block_index_velocity),
  //                             cell_category);

  // TODO:
  //  if(param.has_pml)
  //    inverse_mass_velocity.apply(dst.block(block_index_pml_aux), src.block(block_index_pml_aux),
  //    numbers::pml_material_id);
}

template<int dim, typename Number>
std::vector<
  std::
    tuple<double, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>>
SpatialOperator<dim, Number>::calculate_time_step_lts() const
{
  // TODO: this function is very basic and has to be extended!
  // TODO: WE HAVE TO USE CELL CATEGORIES AND ENSURE CELL BATCHES WITH THE SAME CATEGORY ARE
  // STRICTLY SEPARATED!!!!
  // TODO: ENSURE WE COMMUNICATE IN CASE OF MULTIPLE PROCESSORS

  using scalar = dealii::VectorizedArray<Number>;
  using vector = dealii::Tensor<1, dim, scalar>;

  auto const small_dt = calculate_time_step_cfl();
  auto const large_dt = 2.0 * small_dt;

  std::cerr << "FK: small_dt " << small_dt << std::endl;

  auto const &                     mf = get_matrix_free();
  CellIntegrator<dim, dim, Number> fe_eval(mf);

  std::vector<unsigned int> large_cell_batches;
  std::vector<unsigned int> small_cell_batches;
  // std::vector<unsigned int> small_cell_adjacent_batches;

  for(unsigned int cell = 0; cell < mf.n_cell_batches(); ++cell)
  {
    fe_eval.reinit(cell);
    auto p = fe_eval.quadrature_point(0);

    std::cerr << "x " << p[0] << std::endl;
    std::cerr << "y " << p[1] << std::endl << std::endl;

    auto x_min = std::min_element(p[0].begin(), p[0].end());
    auto y_min = std::min_element(p[1].begin(), p[1].end());
    auto x_max = std::max_element(p[0].begin(), p[0].end());
    auto y_max = std::max_element(p[1].begin(), p[1].end());
    if(*x_min > 0.375 && *x_max < 0.625 && *y_min > 0.375 && *y_max < 0.625)
    {
      std::cerr << "small bactch " << cell << std::endl;
      small_cell_batches.push_back(cell);
    }
    else if(*x_min > 0.25 && *x_max < 0.75 && *y_min > 0.25 && *y_max < 0.75)
    {
      std::cerr << "small and large bactch " << cell << std::endl;
      // small_cell_adjacent_batches.push_back(cell);
      small_cell_batches.push_back(cell);
      large_cell_batches.push_back(cell);
    }
    else
    {
      std::cerr << "large bactch " << cell << std::endl;
      large_cell_batches.push_back(cell);
    }
  }

  FaceIntegrator<dim, dim, Number> fe_face_eval(mf);
  std::vector<unsigned int>        large_face_batches;
  std::vector<unsigned int>        small_face_batches;
  for(unsigned int face = 0; face < mf.n_inner_face_batches(); ++face)
  {
    fe_face_eval.reinit(face);
    auto p = fe_face_eval.quadrature_point(0);

    std::cerr << "face x " << p[0] << std::endl;
    std::cerr << "face y " << p[1] << std::endl << std::endl;

    auto x_min = std::min_element(p[0].begin(), p[0].end());
    auto y_min = std::min_element(p[1].begin(), p[1].end());
    auto x_max = std::max_element(p[0].begin(), p[0].end());
    auto y_max = std::max_element(p[1].begin(), p[1].end());

    // eps to ensure we pick up every face batch, in reality we only need faces that are originated
    // in the same cell category and faces that touch the same cell category and its neighbors
    auto const eps = 0.1;
    if(*x_min > 0.375 - eps && *x_max < 0.625 + eps && *y_min > 0.375 - eps && *y_max < 0.625 + eps)
    {
      std::cerr << "small face bactch " << face << std::endl;
      small_face_batches.push_back(face);
    }
    else
    {
      std::cerr << "large face bactch " << face << std::endl;
      large_face_batches.push_back(face);
    }
  }

  std::vector<unsigned int> small_face_bnd_batches;
  std::vector<unsigned int> large_face_bnd_batches;
  large_face_bnd_batches.resize(mf.n_boundary_face_batches());
  std::iota(large_face_bnd_batches.begin(),
            large_face_bnd_batches.end(),
            mf.n_inner_face_batches());

  // loop over cells of processor
  /*
    // optimization opportuniy. dont store cell_id
    std::vector<std::pair<unsigned int, double>> cell_id_size(mf.n_cell_batches());

    for(unsigned int cell = 0; cell < mf.n_cell_batches(); ++cell)
    {
      scalar dt_min = dealii::make_vectorized_array<Number>(std::numeric_limits<Number>::max());
      fe_eval.reinit(cell);
      for(unsigned int q = 0; q < fe_eval.n_q_points; ++q)
      {
        // TODO: assume spped of sound is 1 for now, we are only interested in cell sizes
        vector c;
        c             = 1.0;
        auto   invJ   = fe_eval.inverse_jacobian(q);
        std::cerr<< invJ << std::endl;
        invJ = transpose(invJ);
        std::cerr<< invJ << std::endl;
        scalar factor = 1.0 / (invJ * c).norm();
        dt_min        = std::min(dt_min, factor);

      }
      auto const dt_min_batch = *std::min_element(dt_min.begin(), dt_min.end());
      cell_id_size[cell]      = std::make_pair(cell, dt_min_batch);

    }

    auto max_el =
      *std::max_element(cell_id_size.begin(),
                        cell_id_size.end(),
                        [](auto const & a, auto const & b) { return a.second < b.second; });

    std::cerr << "Numer of cells in batch " << scalar::size() << std::endl;
    for(auto [i, s] : cell_id_size)
    {
      (s < 0.9 * max_el.second) ? small_cell_batches.push_back(i) : large_cell_batches.push_back(i);
      std::cerr << ((s < 0.9 * max_el.second) ? "small" : "large") << std::endl;
    }
  */


  // TODO: remove this once it works with once cell category

  large_face_batches.resize(mf.n_inner_face_batches());
  std::iota(large_face_batches.begin(), large_face_batches.end(), 0);
  large_cell_batches.resize(mf.n_cell_batches());
  std::iota(large_cell_batches.begin(), large_cell_batches.end(), 0);
  return {{small_dt, {}, {}, {}},
          {large_dt, large_cell_batches, large_face_batches, large_face_bnd_batches}};



  return {{small_dt, small_cell_batches, small_face_batches, small_face_bnd_batches},
          {large_dt, large_cell_batches, large_face_batches, large_face_bnd_batches}};
}


template<int dim, typename Number>
double
SpatialOperator<dim, Number>::calculate_time_step_cfl() const
{
  // In case of mixed-orders use the maximum fe_degree and the corresponding
  // quadrature rule.

  // The time-step size is not adapted every time-step. Thus, we are using
  // a constant function to pass in the speed of sound, even though it is
  // possible to optimize calculate_time_step_cfl_local() for this case.


  auto small_dt =
    this->param.cfl *
    calculate_time_step_cfl_local<dim, Number>(
      get_matrix_free(),
      get_dof_index_velocity(),
      get_quad_index_pressure_velocity(),
      std::make_shared<dealii::Functions::ConstantFunction<dim>>(param.speed_of_sound, dim),
      param.start_time /* will not be used (ConstantFunction) */,
      std::max(param.degree_p, param.degree_u),
      param.cfl_exponent_fe_degree,
      CFLConditionType::VelocityNorm,
      mpi_comm);

  std::cerr << "SMALLLL DTR" << small_dt << std::endl;

  // // testing purposes
  // calculate_time_step_lts<dim, Number>(small_dt, get_matrix_free());


  return small_dt;
}


template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_dof_handler_and_constraints()
{
  fe_p = create_finite_element<dim>(param.grid.element_type, true, 1, param.degree_p);
  fe_u = create_finite_element<dim>(param.grid.element_type, true, dim, param.degree_u);

  // enumerate degrees of freedom
  dof_handler_p.distribute_dofs(*fe_p);
  dof_handler_u.distribute_dofs(*fe_u);

  // close constraints
  constraint_u.close();
  constraint_p.close();

  // Output DoF information
  pcout << "Pressure:" << std::endl;
  print_parameter(pcout, "degree of 1D polynomials", param.degree_p);
  print_parameter(pcout, "number of dofs per cell", fe_p->n_dofs_per_cell());
  print_parameter(pcout, "number of dofs (total)", dof_handler_p.n_dofs());

  pcout << "Velocity:" << std::endl;
  print_parameter(pcout, "degree of 1D polynomials", param.degree_u);
  print_parameter(pcout, "number of dofs per cell", fe_u->n_dofs_per_cell());
  print_parameter(pcout, "number of dofs (total)", dof_handler_u.n_dofs());

  // TODO:
  //  if(param.has_pml)
  //  {
  //    pcout << "PML auxiliary:" << std::endl;
  //    print_parameter(pcout, "degree of 1D polynomials", param.degree_u);
  //    print_parameter(pcout, "number of dofs per cell", fe_u->n_dofs_per_cell());
  //    print_parameter(pcout, "number of dofs (total)", n_pml_cells * fe_u->n_dofs_per_cell());
  //  }

  pcout << "Total:" << std::endl;
  print_parameter(pcout,
                  "number of dofs per cell",
                  fe_p->n_dofs_per_cell() + fe_u->n_dofs_per_cell());
  // TODO:
  //  if(param.has_pml)
  //  {
  //    print_parameter(pcout,
  //                    "number of dofs per PML cell",
  //                    fe_p->n_dofs_per_cell() + 2 * fe_u->n_dofs_per_cell());
  //  }
  print_parameter(pcout, "number of dofs (total)", get_number_of_dofs());

  pcout << std::flush;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_operators()
{
  // inverse mass operator pressure
  {
    InverseMassOperatorData data;
    data.dof_index  = get_dof_index_pressure();
    data.quad_index = get_quad_index_pressure();
    inverse_mass_pressure.initialize(*matrix_free, data);
  }

  // inverse mass operator velocity
  {
    InverseMassOperatorData data;
    data.dof_index  = get_dof_index_velocity();
    data.quad_index = get_quad_index_velocity();
    inverse_mass_velocity.initialize(*matrix_free, data);
  }

  // acoustic operator
  {
    OperatorData<dim> data;
    data.dof_index_pressure   = get_dof_index_pressure();
    data.dof_index_velocity   = get_dof_index_velocity();
    data.quad_index           = get_quad_index_pressure_velocity();
    data.block_index_pressure = block_index_pressure;
    data.block_index_velocity = block_index_velocity;
    data.speed_of_sound       = param.speed_of_sound;
    data.formulation          = param.formulation;
    data.bc                   = boundary_descriptor;
    acoustic_operator.initialize(*matrix_free, data);
  }

  // pml operator
  // TODO:
  // if(param.has_pml)
  // {
  //   PMLOperatorData<dim> data;
  //   data.dof_index_pressure    = get_dof_index_pressure();
  //   data.dof_index_velocity    = get_dof_index_velocity();
  //   data.quad_index            = get_quad_index_pressure_velocity();
  //   data.block_index_pressure  = block_index_pressure;
  //   data.block_index_velocity  = block_index_velocity;
  //   data.block_index_auxiliary = block_index_pml_aux;
  //   data.pml_damping           = field_functions->pml_damping;
  //
  //   pml_operator.initialize(*matrix_free, data);
  // }

  // rhs operator
  if(param.right_hand_side)
  {
    RHSOperatorData<dim> data;
    data.dof_index  = get_dof_index_pressure();
    data.quad_index = get_quad_index_pressure();
    // no source terms are allowed inside a PML, so we have to skip them during the cell loop.
    data.has_pml = false;
    // TODO:
    // data.has_pml       = param.has_pml;
    data.kernel_data.f = field_functions->right_hand_side;
    rhs_operator.initialize(*matrix_free, data);
  }
}

template class SpatialOperator<2, float>;
template class SpatialOperator<3, float>;

template class SpatialOperator<2, double>;
template class SpatialOperator<3, double>;

} // namespace Acoustics
} // namespace ExaDG
