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

#ifndef INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_AB_LTS_BASE_H_
#define INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_AB_LTS_BASE_H_


#include <exadg/time_integration/ab_constants.h>
#include <exadg/time_integration/am_constants.h>
#include <exadg/time_integration/push_back_vectors.h>
#include <exadg/time_integration/time_int_multistep_base.h>
#include <exadg/utilities/print_solver_results.h>

namespace ExaDG
{
/**
 * This class implements the purely explicit Adams--Bashforth--Moulton predictor corrector method.
 */
template<typename Operator, typename VectorType>
class TimeIntAdamsBashforthMoultonBase : public TimeIntMultistepBase
{
  using Number = typename VectorType::value_type;

public:
  TimeIntAdamsBashforthMoultonBase(std::shared_ptr<Operator> pde_operator_in,
                                   double const              start_time_,
                                   double const              end_time_,
                                   unsigned int const        max_number_of_time_steps_,
                                   unsigned int const        order_,
                                   bool const                start_with_low_order_,
                                   bool const                adaptive_time_stepping_,
                                   bool const                local_time_stepping_,
                                   RestartData const &       restart_data_,
                                   MPI_Comm const &          mpi_comm_,
                                   bool const                is_test_)
    : TimeIntMultistepBase(start_time_,
                           end_time_,
                           max_number_of_time_steps_,
                           order_,
                           start_with_low_order_,
                           adaptive_time_stepping_,
                           restart_data_,
                           mpi_comm_,
                           is_test_),
      pde_operator(pde_operator_in),
      // order of predictor can be chosen one below order of corrector
      ab(order_ - 1, start_with_low_order_),
      am(order_, start_with_low_order_),
      vec_evaluated_operators(order_ - 1),
      local_time_stepping(local_time_stepping_)
  {
    AssertThrow(order_ >= 1,
                dealii::ExcMessage("Oder of ABM time integrator has to be at least 1."));
  }

  void
  print_iterations() const
  {
    // explicit time integration -> no iterations
    print_list_of_iterations(pcout, {"Adams-Bashforth-Moulton"}, {0});
  }

  void
  ale_update()
  {
    AssertThrow(false, dealii::ExcMessage("not yet implemented"));
  }

  VectorType const &
  get_solution() const
  {
    return solution;
  }

protected:
  Operator const &
  get_underlying_operator() const
  {
    return *pde_operator;
  }

private:
  void
  update_time_integrator_constants() final
  {
    ab.update(time_step_number, adaptive_time_stepping, time_steps);
    am.update(time_step_number, adaptive_time_stepping, time_steps);
  }

  void
  allocate_vectors() final
  {
    pde_operator->initialize_dof_vector(solution);
    pde_operator->initialize_dof_vector(prediction);

    pde_operator->initialize_dof_vector(evaluated_operator_np);
    for(auto & evaluated_operator : vec_evaluated_operators)
      pde_operator->initialize_dof_vector(evaluated_operator);
  }

  void
  initialize_current_solution() final
  {
    pde_operator->prescribe_initial_conditions(solution, get_time());
  }

  void
  initialize_former_multistep_dof_vectors() final
  {
    if(start_with_low_order)
    {
      if(vec_evaluated_operators.size() > 0)
        pde_operator->evaluate(vec_evaluated_operators[0], solution, get_time());
    }
    else // start with high order
    {
      if(!local_time_stepping)
      {
        // fill evaluated operators
        VectorType temp_sol;
        pde_operator->initialize_dof_vector(temp_sol);
        for(unsigned int i = 0; i < vec_evaluated_operators.size(); ++i)
        {
          pde_operator->prescribe_initial_conditions(temp_sol, get_previous_time(i));
          pde_operator->evaluate(vec_evaluated_operators[i], temp_sol, get_previous_time(i));
        }
      }
      else
      {
        std::vector<CellBatchInfo> infos;
        infos.push_back({get_time_step_size() / 4.0, {1}, 2});
        infos.push_back({get_time_step_size() / 2.0, {2, 3}, 4});
        infos.push_back({get_time_step_size() / 1.0, {4, 5}, dealii::numbers::invalid_material_id});
        std::sort(infos.begin(),
                  infos.end(),
                  [](auto const & a, auto const & b) { return a.dt < b.dt; });

        // fill evaluated operators
        VectorType temp_sol;
        pde_operator->initialize_dof_vector(temp_sol);
        for(unsigned int i = 0; i < vec_evaluated_operators.size(); ++i)
        {
          for(auto info : infos)
          {
            double const previous_time = get_time() - i * info.dt;
            pde_operator->prescribe_initial_conditions(temp_sol, previous_time);
            for(auto cell_category : info.cell_categories)
            {
              pde_operator->evaluate(vec_evaluated_operators[i],
                                     temp_sol,
                                     get_previous_time(i),
                                     cell_category);
            }
          }
        }
      }
    }
  }

  void
  setup_derived() final
  {
  }

  void
  do_timestep_predict()
  {
    dealii::Timer timer;
    timer.restart();

    predrict_solution(prediction, solution, vec_evaluated_operators);

    timer_tree->insert({"Timeloop", "Adams-Bashforth-Moulton"}, timer.wall_time());
  }

  void
  do_timestep_correct()
  {
    dealii::Timer timer;
    timer.restart();

    // evaluate operator given the predicted solution
    pde_operator->evaluate(evaluated_operator_np, prediction, get_next_time());
    // correct solution
    correct_solution(solution, evaluated_operator_np, vec_evaluated_operators);
    // correct operator by evaluating operator with correct solution
    pde_operator->evaluate(evaluated_operator_np, solution, get_next_time());

    // write output
    if(this->print_solver_info() and not(this->is_test))
    {
      pcout << std::endl << "Adams-Bashforth-Moulton:";
      print_wall_time(pcout, timer.wall_time());
    }

    timer_tree->insert({"Timeloop", "Adams-Bashforth-Moulton"}, timer.wall_time());
  }

  void
  do_lts_timestep_predict(VectorType &                      p,
                          VectorType const &                s,
                          double const                      dt,
                          dealii::types::material_id        cell_category,
                          ABTimeIntegratorConstants const & ab_in)
  {
    predrict_solution_lts(p, s, vec_evaluated_operators, dt, cell_category, ab_in);
  }

  void
  do_lts_timestep_correct(VectorType &               s,
                          VectorType const &         p,
                          double                     dt,
                          double                     next_time,
                          dealii::types::material_id cell_category)
  {
    // evaluate operator given the predicted solution
    pde_operator->evaluate(evaluated_operator_np, p, next_time, cell_category);
    // correct solution
    correct_solution_lts(s, evaluated_operator_np, vec_evaluated_operators, dt, cell_category);
    // correct operator by evaluating operator with correct solution
    pde_operator->evaluate(evaluated_operator_np, s, next_time, cell_category);
  }

  struct CellBatchInfo
  {
    double                                  dt;
    std::vector<dealii::types::material_id> cell_categories;
    dealii::types::material_id              attached_cell_category;
  };

  void
  do_lts_sub_timestep_solve(VectorType &                                    dst,
                            VectorType const &                              src,
                            VectorType const &                              src_old,
                            double                                          time,
                            double                                          dt,
                            std::vector<dealii::types::material_id> const & categories,
                            dealii::types::material_id const                attached,
                            ABTimeIntegratorConstants                       ab_attached)
  {
    if(attached != dealii::numbers::invalid_material_id)
    {
      do_lts_timestep_predict(dst, src_old, dt, attached, ab_attached);
    }
    for(auto const category : categories)
    {
      do_lts_timestep_predict(dst, src, dt, category, ab);
    }
    for(auto const category : categories)
    {
      pde_operator->evaluate(evaluated_operator_np, dst, time + dt, category);
    }
    for(auto const category : categories)
    {
      prepare_vectors_for_next_timestep_lts(category);
    }
  }

  void
  do_lts_sub_timesteps_recursively(VectorType &                       dst,
                                   VectorType const &                 src,
                                   double                             t,
                                   double const                       dt,
                                   unsigned int                       index,
                                   std::vector<CellBatchInfo> const & infos,
                                   ABTimeIntegratorConstants          ab_intermediate)
  {
    if(index > 0)
    {
      auto intermediate = src;

      do_lts_sub_timesteps_recursively(
        intermediate, src, t, 0.5 * dt, index - 1, infos, ab_intermediate);

      do_lts_sub_timestep_solve(intermediate,
                                src,
                                src,
                                t,
                                dt,
                                infos[index].cell_categories,
                                infos[index].attached_cell_category,
                                ab_intermediate);

      do_lts_sub_timesteps_recursively(
        dst, intermediate, t + dt, 0.5 * dt, index - 1, infos, ab_intermediate);

      do_lts_sub_timestep_solve(dst,
                                intermediate,
                                src,
                                t + dt,
                                dt,
                                infos[index].cell_categories,
                                infos[index].attached_cell_category,
                                ab);
    }
    else
    {
      auto intermediate = src;
      do_lts_sub_timestep_solve(intermediate,
                                src,
                                src,
                                t,
                                dt,
                                infos[index].cell_categories,
                                infos[index].attached_cell_category,
                                ab_intermediate);

      do_lts_sub_timestep_solve(dst,
                                intermediate,
                                src,
                                t + dt,
                                dt,
                                infos[index].cell_categories,
                                infos[index].attached_cell_category,
                                ab);
    }
  }

  void
  do_timestep_solve() final
  {
    if(!local_time_stepping)
    {
      do_timestep_predict();
      do_timestep_correct();
    }
    else
    {
      double const t  = get_time();
      double const dt = get_time_step_size();

      VectorType solution_old = solution;

      std::vector<CellBatchInfo> infos;
      infos.push_back({get_time_step_size() / 4.0, {1}, 2});
      infos.push_back({get_time_step_size() / 2.0, {2, 3}, 4});
      infos.push_back({get_time_step_size() / 1.0, {4, 5}, dealii::numbers::invalid_material_id});
      std::sort(infos.begin(),
                infos.end(),
                [](auto const & a, auto const & b) { return a.dt < b.dt; });

      auto ts = get_time_step_vector();
      ts[0] *= 0.5;
      auto ab_intermediate = ab;
      ab_intermediate.update(get_current_order(), true, ts);

      do_lts_sub_timesteps_recursively(
        solution, solution_old, t, 0.5 * dt, infos.size() - 2, infos, ab_intermediate);

      do_lts_sub_timestep_solve(solution,
                                solution_old,
                                solution_old,
                                t,
                                dt,
                                infos.back().cell_categories,
                                dealii::numbers::invalid_material_id,
                                ab);
    }
  }

  void
  prepare_vectors_for_next_timestep_lts(dealii::types::material_id cell_category)
  {
    if(vec_evaluated_operators.size() > 0)
    {
      for(std::size_t i = vec_evaluated_operators.size() - 1; i > 0; --i)
      {
        pde_operator->copy_dofs_of_cell_category(vec_evaluated_operators[i],
                                                 vec_evaluated_operators[i - 1],
                                                 cell_category);
      }
      pde_operator->copy_dofs_of_cell_category(vec_evaluated_operators[0],
                                               evaluated_operator_np,
                                               cell_category);
    }
  }

  void
  correct_solution(VectorType &                    dst,
                   VectorType const &              op_np,
                   std::vector<VectorType> const & ops) const
  {
    dst.add(static_cast<Number>(get_time_step_size() * am.get_gamma0()), op_np);
    for(unsigned int i = 0; i < this->am.get_order() - 1; ++i)
      dst.add(static_cast<Number>(get_time_step_size() * am.get_alpha(i)), ops[i]);
  }

  void
  predrict_solution(VectorType &                    dst,
                    VectorType const &              src,
                    std::vector<VectorType> const & ops) const
  {
    dst = src;
    for(unsigned int i = 0; i < this->ab.get_order(); ++i)
      dst.add(static_cast<Number>(get_time_step_size() * ab.get_alpha(i)), ops[i]);
  }

  void
  correct_solution_lts(VectorType &                    dst,
                       VectorType const &              op_np,
                       std::vector<VectorType> const & ops,
                       double                          dt,
                       dealii::types::material_id      cell_category) const
  {
    pde_operator->add_dofs_of_cell_category(dst,
                                            static_cast<Number>(dt * am.get_gamma0()),
                                            op_np,
                                            cell_category);

    for(unsigned int i = 0; i < this->am.get_order() - 1; ++i)
    {
      pde_operator->add_dofs_of_cell_category(dst,
                                              static_cast<Number>(dt * am.get_alpha(i)),
                                              ops[i],
                                              cell_category);
    }
  }

  void
  predrict_solution_lts(VectorType &                      dst,
                        VectorType const &                src,
                        std::vector<VectorType> const &   ops,
                        double                            dt,
                        dealii::types::material_id        cell_category,
                        ABTimeIntegratorConstants const & ab_in) const
  {
    pde_operator->copy_dofs_of_cell_category(dst, src, cell_category);
    for(unsigned int i = 0; i < ab_in.get_order(); ++i)
    {
      pde_operator->add_dofs_of_cell_category(dst,
                                              static_cast<Number>(dt * ab_in.get_alpha(i)),
                                              ops[i],
                                              cell_category);
    }
  }

  void
  prepare_vectors_for_next_timestep() final
  {
    // for local timestepping this is already done during do_timestep_solve()
    if(vec_evaluated_operators.size() > 0 && !local_time_stepping)
    {
      push_back(vec_evaluated_operators);
      std::swap(vec_evaluated_operators[0], evaluated_operator_np);
    }
  }

  void
  read_restart_vectors(boost::archive::binary_iarchive & ia) final
  {
    ia >> solution;
    ia >> prediction;
  }

  void
  write_restart_vectors(boost::archive::binary_oarchive & oa) const final
  {
    oa << solution;
    oa << prediction;
  }

  void
  solve_steady_problem() final
  {
    AssertThrow(false, dealii::ExcMessage("Steady not implemented."));
  }

private:
  // spatial pde operator
  std::shared_ptr<Operator> pde_operator;

  // Time integration constants
  ABTimeIntegratorConstants ab;
  AMTimeIntegratorConstants am;

  // solution vector
  VectorType solution;
  // temporary vector to store prediction
  VectorType prediction;

  // store evaluated operators from previous time steps
  VectorType              evaluated_operator_np;
  std::vector<VectorType> vec_evaluated_operators;

  bool const local_time_stepping;
};

} // namespace ExaDG

#endif /* INCLUDE_EXADG_TIME_INTEGRATION_TIME_INT_AB_LTS_BASE_H_*/
