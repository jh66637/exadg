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

#ifndef INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_ERROR_CALCULATION_H_
#define INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_ERROR_CALCULATION_H_

// deal.II
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/la_parallel_vector.h>

// ExaDG
#include <exadg/utilities/print_functions.h>

namespace ExaDG
{
template<int dim>
struct ErrorCalculationData
{
  ErrorCalculationData()
    : analytical_solution_available(false),
      calculate_relative_errors(true),
      calculate_H1_seminorm_error(false),
      error_calc_start_time(std::numeric_limits<double>::max()),
      error_calc_interval_time(std::numeric_limits<double>::max()),
      calculate_every_time_steps(std::numeric_limits<unsigned int>::max()),
      write_errors_to_file(false),
      directory("output/"),
      name("all fields")
  {
  }

  void
  print(dealii::ConditionalOStream & pcout, bool unsteady)
  {
    print_parameter(pcout, "Calculate error", analytical_solution_available);
    if(analytical_solution_available == true && unsteady == true)
    {
      print_parameter(pcout, "Calculate relative errors", calculate_relative_errors);
      print_parameter(pcout, "Calculate H1-seminorm error", calculate_H1_seminorm_error);
      print_parameter(pcout, "Error calculation start time", error_calc_start_time);
      print_parameter(pcout, "Error calculation interval time", error_calc_interval_time);
      print_parameter(pcout, "Calculate error every time steps", calculate_every_time_steps);
      print_parameter(pcout, "Write errors to file", write_errors_to_file);
      if(write_errors_to_file)
        print_parameter(pcout, "Directory", directory);
      print_parameter(pcout, "Name", name);
    }
  }

  // to calculate the error an analytical solution to the problem has to be available
  bool analytical_solution_available;

  std::shared_ptr<dealii::Function<dim>> analytical_solution;

  // relative or absolute errors?
  // If calculate_relative_errors == false, this implies that absolute errors are calculated
  bool calculate_relative_errors;

  // by default, only the L2-error is computed. Other norms have to be explicitly specified by the
  // user.
  bool calculate_H1_seminorm_error;

  // before then no error calculation will be performed
  double error_calc_start_time;

  // specifies the time interval in which error calculation is performed
  double error_calc_interval_time;

  // calculate error every time steps
  unsigned int calculate_every_time_steps;

  // write errors to file?
  bool write_errors_to_file;

  // directory and name (used as filename and as identifier for screen output)
  std::string directory;
  std::string name;
};

template<int dim, typename Number>
class ErrorCalculator
{
public:
  typedef dealii::LinearAlgebra::distributed::Vector<Number> VectorType;

  ErrorCalculator(MPI_Comm const & comm);

  void
  setup(dealii::DoFHandler<dim> const &   dof_handler,
        dealii::Mapping<dim> const &      mapping,
        ErrorCalculationData<dim> const & error_data);

  void
  evaluate(VectorType const & solution, double const & time, int const & time_step_number);

private:
  void
  do_evaluate(VectorType const & solution_vector, double const time  ,    const std::set<dealii::types::material_id>& material_ids={});

  MPI_Comm const mpi_comm;

  unsigned int error_counter;
  bool         reset_counter;

  bool clear_files_L2, clear_files_H1_seminorm;

  dealii::SmartPointer<dealii::DoFHandler<dim> const> dof_handler;
  dealii::SmartPointer<dealii::Mapping<dim> const>    mapping;

  ErrorCalculationData<dim> error_data;
};

} // namespace ExaDG

#endif /* INCLUDE_COMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_ERROR_CALCULATION_H_ */
