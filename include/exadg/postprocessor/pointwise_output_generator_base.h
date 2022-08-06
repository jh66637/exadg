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

#ifndef INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_BASE_H_
#define INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_BASE_H_

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/point.h>

#ifdef DEAL_II_WITH_HDF5
#  include <deal.II/base/hdf5.h>
#endif

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>

namespace ExaDG
{
template<int dim>
struct PointwiseOutputDataBase
{
  using point_value_type = typename dealii::Point<dim>::value_type;

  PointwiseOutputDataBase();

  bool        write_output;
  std::string directory;
  std::string filename;

  bool update_points_before_evaluation;

  double start_time;
  double end_time;
  double interval_time;

  std::vector<dealii::Point<dim>> evaluation_points;

  void
  print(dealii::ConditionalOStream & pcout) const;
};

template<int dim, typename Number>
class PointwiseOutputGeneratorBase
{
public:
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<Number>;
  using point_value_type = typename PointwiseOutputDataBase<dim>::point_value_type;

  void
  evaluate(VectorType const & solution, double const & time, int const & time_step_number);

protected:
  PointwiseOutputGeneratorBase(MPI_Comm const & comm);

  virtual ~PointwiseOutputGeneratorBase() = default;

  void
  setup_base(dealii::DoFHandler<dim> const &      dof_handler_in,
             dealii::Mapping<dim> const &         mapping_in,
             PointwiseOutputDataBase<dim> const & pointwise_output_data_in);

  void
  add_scalar_quantity_dataset(std::string const & name);

  void
  add_vector_quantity_datasets(std::string const & name);

  void
  write_scalar_quantity(std::string const &         name,
                        std::vector<Number> const & values,
                        unsigned int const          first_component = 0)
  {
    AssertThrow(first_component == 0,
                dealii::ExcMessage("This function can only be called with first_component==0."));

    componentwise_result = std::move(dealii::Vector<Number>(values.begin(), values.end()));
    write_component(name, componentwise_result);
  }

  template<int n_components>
  void
  write_scalar_quantity(std::string const &                                          name,
                        std::vector<dealii::Tensor<1, n_components, Number>> const & values,
                        unsigned int const first_component = 0)
  {
    this->extract_component(componentwise_result, values, first_component);
    write_component(name, componentwise_result);
  }

  template<int n_components>
  void
  write_vectorial_quantity(std::string const &                                          name,
                           std::vector<dealii::Tensor<1, n_components, Number>> const & values,
                           unsigned int const first_component = 0)
  {
    for(unsigned int d = 0; d < dim; ++d)
    {
      this->extract_component(componentwise_result, values, d + first_component);
      write_component(name + std::to_string(d), componentwise_result);
    }
  }

  template<int n_components>
  [[nodiscard]] std::vector<
    typename dealii::FEPointEvaluation<n_components, dim, dim, Number>::value_type>
  compute_point_values(VectorType const & solution) const
  {
    return dealii::VectorTools::point_values<n_components>(*remote_evaluator,
                                                           *dof_handler,
                                                           solution);
  }

private:
  virtual void
  do_evaluate(VectorType const & solution) = 0;

  void
  setup_remote_evaluator();

  void
  reinit_remote_evaluator();

  void
  create_hdf5_file();

#ifdef DEAL_II_WITH_HDF5
  void
  add_evaluation_points_dataset(dealii::HDF5::Group & group, std::string const & name);
#endif

  void
  fill_evaluation_points(std::string const & name);

#ifdef DEAL_II_WITH_HDF5
  void
  add_time_dataset(dealii::HDF5::Group & group, std::string const & name);
#endif

  void
  write_time(double time);

  void
  write_component(std::string const & name, dealii::Vector<Number> const & componentwise_result);

  template<int n_components>
  void
  extract_component(dealii::Vector<Number> &                                     dst,
                    std::vector<dealii::Tensor<1, n_components, Number>> const & values,
                    unsigned int const                                           comp) const
  {
    AssertIndexRange(comp, n_components);

    auto iter = dst.begin();
    for(auto const & val : values)
    {
      *iter = val[comp];
      ++iter;
    }
  }

  MPI_Comm const mpi_comm;
  unsigned int   counter;
  bool           reset_counter;

  dealii::SmartPointer<dealii::DoFHandler<dim> const>                 dof_handler;
  dealii::SmartPointer<dealii::Mapping<dim> const>                    mapping;
  PointwiseOutputDataBase<dim>                                        pointwise_output_data;
  dealii::Vector<Number>                                              componentwise_result;
  unsigned int                                                        n_out_samples;
  std::shared_ptr<dealii::Utilities::MPI::RemotePointEvaluation<dim>> remote_evaluator;
#ifdef DEAL_II_WITH_HDF5
  std::unique_ptr<dealii::HDF5::File> hdf5_file;
#endif
};

} // namespace ExaDG

#endif /* INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_BASE_H_*/
