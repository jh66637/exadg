/*
 * driver.h
 *
 *  Created on: 25.03.2020
 *      Author: fehn
 */

#ifndef INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_DRIVER_H_
#define INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_DRIVER_H_

#include <exadg/compressible_navier_stokes/spatial_discretization/dg_operator.h>
#include <exadg/compressible_navier_stokes/time_integration/time_int_explicit_runge_kutta.h>
#include <exadg/compressible_navier_stokes/user_interface/analytical_solution.h>
#include <exadg/compressible_navier_stokes/user_interface/application_base.h>
#include <exadg/compressible_navier_stokes/user_interface/boundary_descriptor.h>
#include <exadg/compressible_navier_stokes/user_interface/field_functions.h>
#include <exadg/compressible_navier_stokes/user_interface/input_parameters.h>
#include <exadg/functions_and_boundary_conditions/verify_boundary_conditions.h>
#include <exadg/grid/mapping_degree.h>
#include <exadg/grid/mesh.h>
#include <exadg/matrix_free/matrix_free_wrapper.h>
#include <exadg/utilities/print_general_infos.h>

namespace ExaDG
{
namespace CompNS
{
using namespace dealii;

// Select the operator to be applied
enum class Operator
{
  ConvectiveTerm,
  ViscousTerm,
  ViscousAndConvectiveTerms,
  InverseMassMatrix,
  InverseMassMatrixDstDst,
  VectorUpdate,
  EvaluateOperatorExplicit
};

inline std::string
enum_to_string(Operator const enum_type)
{
  std::string string_type;

  switch(enum_type)
  {
    // clang-format off
    case Operator::ConvectiveTerm:            string_type = "ConvectiveTerm";           break;
    case Operator::ViscousTerm:               string_type = "ViscousTerm";              break;
    case Operator::ViscousAndConvectiveTerms: string_type = "ViscousAndConvectiveTerms";break;
    case Operator::InverseMassMatrix:         string_type = "InverseMassMatrix";        break;
    case Operator::InverseMassMatrixDstDst:   string_type = "InverseMassMatrixDstDst";  break;
    case Operator::VectorUpdate:              string_type = "VectorUpdate";             break;
    case Operator::EvaluateOperatorExplicit:  string_type = "EvaluateOperatorExplicit"; break;

    default:AssertThrow(false, ExcMessage("Not implemented.")); break;
      // clang-format on
  }

  return string_type;
}

inline void
string_to_enum(Operator & enum_type, std::string const string_type)
{
  // clang-format off
  if     (string_type == "ConvectiveTerm")            enum_type = Operator::ConvectiveTerm;
  else if(string_type == "ViscousTerm")               enum_type = Operator::ViscousTerm;
  else if(string_type == "ViscousAndConvectiveTerms") enum_type = Operator::ViscousAndConvectiveTerms;
  else if(string_type == "InverseMassMatrix")         enum_type = Operator::InverseMassMatrix;
  else if(string_type == "InverseMassMatrixDstDst")   enum_type = Operator::InverseMassMatrixDstDst;
  else if(string_type == "VectorUpdate")              enum_type = Operator::VectorUpdate;
  else if(string_type == "EvaluateOperatorExplicit")  enum_type = Operator::EvaluateOperatorExplicit;
  else AssertThrow(false, ExcMessage("Unknown operator type. Not implemented."));
  // clang-format on
}

inline unsigned int
get_dofs_per_element(std::string const & input_file,
                     unsigned int const  dim,
                     unsigned int const  degree)
{
  (void)input_file;

  unsigned int const dofs_per_element = (dim + 2) * std::pow(degree + 1, dim);

  return dofs_per_element;
}

template<int dim, typename Number = double>
class Driver
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  Driver(MPI_Comm const & comm);

  void
  setup(std::shared_ptr<ApplicationBase<dim, Number>> application,
        unsigned int const                            degree,
        unsigned int const                            refine_space,
        unsigned int const                            refine_time,
        bool const                                    is_test,
        bool const                                    is_throughput_study);

  void
  solve();

  void
  print_performance_results(double const total_time, bool const is_test) const;

  std::tuple<unsigned int, types::global_dof_index, double>
  apply_operator(unsigned int const  degree,
                 std::string const & operator_type,
                 unsigned int const  n_repetitions_inner,
                 unsigned int const  n_repetitions_outer,
                 bool const          is_test) const;

private:
  // MPI communicator
  MPI_Comm const & mpi_comm;

  // output to std::cout
  ConditionalOStream pcout;

  // application
  std::shared_ptr<ApplicationBase<dim, Number>> application;

  InputParameters param;

  // triangulation
  std::shared_ptr<parallel::TriangulationBase<dim>> triangulation;

  // mapping
  std::shared_ptr<Mesh<dim>> mesh;

  // periodic boundaries
  std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;

  std::shared_ptr<FieldFunctions<dim>>           field_functions;
  std::shared_ptr<BoundaryDescriptor<dim>>       boundary_descriptor_density;
  std::shared_ptr<BoundaryDescriptor<dim>>       boundary_descriptor_velocity;
  std::shared_ptr<BoundaryDescriptor<dim>>       boundary_descriptor_pressure;
  std::shared_ptr<BoundaryDescriptorEnergy<dim>> boundary_descriptor_energy;

  std::shared_ptr<MatrixFreeData<dim, Number>> matrix_free_data;
  std::shared_ptr<MatrixFree<dim, Number>>     matrix_free;

  std::shared_ptr<DGOperator<dim, Number>> comp_navier_stokes_operator;

  std::shared_ptr<PostProcessorBase<dim, Number>> postprocessor;

  std::shared_ptr<TimeIntExplRK<Number>> time_integrator;

  // Computation time (wall clock time)
  mutable TimerTree timer_tree;
};

} // namespace CompNS
} // namespace ExaDG

#endif /* INCLUDE_EXADG_COMPRESSIBLE_NAVIER_STOKES_DRIVER_H_ */
