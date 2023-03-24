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


// deal.II
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/point.h>

#ifdef DEAL_II_WITH_HDF5
#  include <deal.II/base/hdf5.h>
#endif

#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/numerics/vector_tools.h>


// ExaDG
#include <exadg/postprocessor/time_control.h>
#include <exadg/utilities/create_directories.h>
#include <exadg/utilities/extract_component_from_tensors.h>

namespace ExaDG
{
template<int dim>
struct PointwiseOutputDataBase
{
  using point_value_type = typename dealii::Point<dim>::value_type;

  PointwiseOutputDataBase();

  TimeControlData time_control_data;

  std::string directory;
  std::string filename;

  bool update_points_before_evaluation;

  std::vector<dealii::Point<dim>> evaluation_points;

  void
  print(dealii::ConditionalOStream & pcout) const;
};

template<int dim, typename VectorType>
class PointwiseOutputGeneratorBase
{
public:
  using Number = typename VectorType::value_type;

  using point_value_type = typename PointwiseOutputDataBase<dim>::point_value_type;

  void
  evaluate(VectorType const & solution, double const time, bool const unsteady)
  {
    AssertThrow(unsteady, dealii::ExcMessage("Only implemented for the unsteady case."));

    if(first_evaluation)
    {
      first_evaluation = false;
      AssertThrow(time_control.get_counter() == 0,
                  dealii::ExcMessage(
                    "Only implemented in the case that the simulation is not restarted"));
    }

    if(pointwise_output_data.update_points_before_evaluation)
      reinit_remote_evaluator();

    write_time(time);
    do_evaluate(solution);
  }

  TimeControl time_control;

protected:
  PointwiseOutputGeneratorBase(MPI_Comm const & comm) : mpi_comm(comm), first_evaluation(true)
  {
  }

  virtual ~PointwiseOutputGeneratorBase() = default;

  void
  setup_base(dealii::Triangulation<dim> const &   triangulation_in,
             dealii::Mapping<dim> const &         mapping_in,
             PointwiseOutputDataBase<dim> const & pointwise_output_data_in)
  {
#ifdef DEAL_II_WITH_HDF5
    triangulation         = &triangulation_in;
    pointwise_output_data = pointwise_output_data_in;


    AssertThrow(
      (get_unsteady_evaluation_type(pointwise_output_data_in.time_control_data) ==
       TimeControlData::UnsteadyEvalType::Interval) ||
        (get_unsteady_evaluation_type(pointwise_output_data_in.time_control_data) ==
         TimeControlData::UnsteadyEvalType::None),
      dealii::ExcMessage(
        "This module can currently only be used with time TimeControlData::UnsteadyEvalType::Interval"));

    time_control.setup(pointwise_output_data_in.time_control_data);

    if(pointwise_output_data.time_control_data.is_active &&
       pointwise_output_data.evaluation_points.size() > 0)
    {
      mapping = &mapping_in;

      // allocate memory for hyperslab
      componentwise_result.reinit(pointwise_output_data.evaluation_points.size());

      // number of samples to write into file
      n_out_samples = 1 + static_cast<unsigned int>(
                            std::ceil((pointwise_output_data.time_control_data.end_time -
                                       pointwise_output_data.time_control_data.start_time) /
                                      pointwise_output_data.time_control_data.trigger_interval));

      setup_remote_evaluator();

      create_hdf5_file();

      {
        auto group = hdf5_file->create_group("GeneralInformation");
        add_evaluation_points_dataset(group, "EvaluationPoints");
        write_evaluation_points("GeneralInformation/EvaluationPoints");
      }

      {
        auto group = hdf5_file->create_group("PhysicalInformation");
        add_time_dataset(group, "Time");
      }
    }
#else
    (void)triangulation_in;
    (void)mapping_in;
    (void)pointwise_output_data_in;
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }


  void
  add_quantity(std::string const & name, unsigned int const n_components)
  {
    AssertThrow(n_components > 0, dealii::ExcMessage("n_components has to be > 0."));

    auto const & [it, success] = name_to_components.try_emplace(name, n_components);
    AssertThrow(success, dealii::ExcMessage("Name already given to quantity dataset."));

#ifdef DEAL_II_WITH_HDF5
    auto group = hdf5_file->open_group("PhysicalInformation");
    for(unsigned int comp = 0; comp < n_components; ++comp)
    {
      group.template create_dataset<Number>(
        (n_components == 1) ? name : name + std::to_string(comp),
        std::vector<hsize_t>{pointwise_output_data.evaluation_points.size(), 1 * n_out_samples});
    }
#else
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }


  template<int n_components>
  void
  write_quantity(std::string const &                                          name,
                 std::vector<dealii::Tensor<1, n_components, Number>> const & values,
                 unsigned int const                                           first_component)
  {
    unsigned int const components = name_to_components.at(name);
    for(unsigned int comp = 0; comp < components; ++comp)
    {
      extract_component_from_tensors(componentwise_result, values, first_component + comp);
      write_component((components == 1) ? name : name + std::to_string(comp), componentwise_result);
    }
  }

  void
  write_quantity(std::string const & name, std::vector<Number> const & values)
  {
    write_component(name, values);
  }

  template<int n_components>
  [[nodiscard]] std::vector<
    typename dealii::FEPointEvaluation<n_components, dim, dim, Number>::value_type>
  compute_point_values(dealii::LinearAlgebra::distributed::Vector<Number> const & solution,
                       dealii::DoFHandler<dim> const &                            dof_handler) const
  {
    return dealii::VectorTools::point_values<n_components>(*remote_evaluator,
                                                           dof_handler,
                                                           solution);
  }

private:
  virtual void
  do_evaluate(VectorType const & solution) = 0;

  void
  setup_remote_evaluator()
  {
    remote_evaluator =
      std::make_shared<dealii::Utilities::MPI::RemotePointEvaluation<dim>>(1e-6, false, 0);
    reinit_remote_evaluator();
  }

  void
  reinit_remote_evaluator()
  {
    remote_evaluator->reinit(pointwise_output_data.evaluation_points, *triangulation, *mapping);
    AssertThrow(remote_evaluator->all_points_found(),
                dealii::ExcMessage("Not all remote points found."));
  }

  void
  create_hdf5_file()
  {
#ifdef DEAL_II_WITH_HDF5
    ExaDG::create_directories(pointwise_output_data.directory, mpi_comm);
    hdf5_file = std::make_unique<dealii::HDF5::File>(pointwise_output_data.directory +
                                                       pointwise_output_data.filename + ".h5",
                                                     dealii::HDF5::File::FileAccessMode::create,
                                                     mpi_comm);
#else
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }

  void
  write_evaluation_points(std::string const & name)
  {
#ifdef DEAL_II_WITH_HDF5
    auto dataset = hdf5_file->open_dataset(name);

    if(dealii::Utilities::MPI::this_mpi_process(mpi_comm) == 0)
    {
      dealii::FullMatrix<point_value_type> evaluation_points(
        dataset.get_dimensions()[0],
        dataset.get_dimensions()[1],
        &pointwise_output_data.evaluation_points[0][0]);

      std::vector<hsize_t> hyperslab_offset = {0, 0};

      dataset.write_hyperslab(evaluation_points, hyperslab_offset, dataset.get_dimensions());
    }
    else
    {
      dataset.template write_none<point_value_type>();
    }
#else
    (void)name;
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }

#ifdef DEAL_II_WITH_HDF5
  void
  add_evaluation_points_dataset(dealii::HDF5::Group & group, std::string const & name)
  {
    std::vector<hsize_t> hdf5_point_dims{pointwise_output_data.evaluation_points.size(), dim};

    group.template create_dataset<point_value_type>(name, hdf5_point_dims);
  }

  void
  add_time_dataset(dealii::HDF5::Group & group, std::string const & name)
  {
    group.template create_dataset<Number>(name, std::vector<hsize_t>{1, 1 * n_out_samples});
  }
#endif

  void
  write_time(double time)
  {
#ifdef DEAL_II_WITH_HDF5
    auto dataset = hdf5_file->open_dataset("PhysicalInformation/Time");

    std::vector<hsize_t> hyperslab_offset = {0, time_control.get_counter()};

    if(dealii::Utilities::MPI::this_mpi_process(mpi_comm) == 0)
    {
      dataset.write_hyperslab(std::vector<Number>{static_cast<Number>(time)},
                              hyperslab_offset,
                              std::vector<hsize_t>{1, 1});
    }
    else
    {
      dataset.template write_none<Number>();
    }
#else
    (void)time;
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }

  template<typename ComponentwiseContainerType>
  void
  write_component(std::string const & name, ComponentwiseContainerType const & componentwise_result)
  {
#ifdef DEAL_II_WITH_HDF5
    auto dataset = hdf5_file->open_dataset("PhysicalInformation/" + name);

    std::vector<hsize_t> hyperslab_offset = {0, time_control.get_counter()};
    std::vector<hsize_t> hyperslab_dim    = {pointwise_output_data.evaluation_points.size(), 1};

    if(dealii::Utilities::MPI::this_mpi_process(mpi_comm) == 0)
      dataset.write_hyperslab(componentwise_result, hyperslab_offset, hyperslab_dim);
    else
      dataset.template write_none<Number>();
#else
    (void)name;
    (void)componentwise_result;
    AssertThrow(false, dealii::ExcMessage("deal.II is not compiled with HDF5!"));
#endif
  }

  MPI_Comm const mpi_comm;

  dealii::SmartPointer<dealii::Triangulation<dim> const>              triangulation;
  dealii::SmartPointer<dealii::Mapping<dim> const>                    mapping;
  PointwiseOutputDataBase<dim>                                        pointwise_output_data;
  dealii::Vector<Number>                                              componentwise_result;
  unsigned int                                                        n_out_samples;
  std::shared_ptr<dealii::Utilities::MPI::RemotePointEvaluation<dim>> remote_evaluator;
  bool                                                                first_evaluation;

  std::map<std::string, unsigned int> name_to_components;

#ifdef DEAL_II_WITH_HDF5
  std::unique_ptr<dealii::HDF5::File> hdf5_file;
#endif
};

} // namespace ExaDG

#endif /* INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_POINTWISE_OUTPUT_GENERATOR_BASE_H_*/
