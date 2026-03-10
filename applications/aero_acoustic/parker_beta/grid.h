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

#ifndef APPLICATIONS_GRID_TOOLS_MESH_TESTCASE_PAUL_H_
#define APPLICATIONS_GRID_TOOLS_MESH_TESTCASE_PAUL_H_

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

namespace ExaDG
{
namespace Parker
{
double const       obstacle_length      = 193.6e-3;
double const       r                    = 0.5 * 12.1e-3;
double const       plate_length         = obstacle_length - r;
unsigned int const elements_along_plate = 8;
double const       h_ref                = plate_length / double(elements_along_plate);

double const coarsening_levels       = 1;
double const h_ref_outer             = std::pow(2, coarsening_levels) * h_ref;
double const front                   = -5.0 * h_ref_outer - h_ref;
double const end                     = plate_length + 15.0 * h_ref_outer + h_ref;
double const original_channel_height = 0.244;
double const channel_height = std::floor(original_channel_height / h_ref_outer) * h_ref_outer;
double const channel_width  = 4.0 * Parker::h_ref;

template<int dim>
std::vector<unsigned int>
balance_refinements(dealii::Point<dim> const & p1, dealii::Point<dim> const & p2, double const h)
{
  std::vector<unsigned int> repetitions(dim);
  for(uint d = 0; d < dim; ++d)
  {
    repetitions[d] = (uint)std::floor((std::abs(p2[d] - p1[d]) + 1e-6) / h);
    if(repetitions[d] == 0)
      repetitions[d] = 1;
  }

  return repetitions;
}

template<int dim>
void
create_obstacle(dealii::Triangulation<dim> & obstacle, double const channel_width)
{
  // cylindrical front
  dealii::Triangulation<dim> front;
  {
    dealii::Triangulation<dim> tria;
    dealii::GridGenerator::hyper_cube_with_cylindrical_hole(
      tria, r, h_ref, channel_width, std::floor(channel_width / h_ref));

    using CellIter = typename dealii::Triangulation<dim>::active_cell_iterator;
    std::set<CellIter> cells_to_remove;
    for(auto const & cell : tria.active_cell_iterators())
      if(cell->center()[0] > 0.0)
        cells_to_remove.insert(cell);

    dealii::GridGenerator::create_triangulation_with_removed_cells(tria, cells_to_remove, front);

    front.set_all_manifold_ids(-1);
    for(auto const & face : front.active_face_iterators())
      if(face->at_boundary())
      {
        for(unsigned int l = 0; l < face->n_lines(); ++l)
        {
          dealii::Point<2> center_yz{face->line(l)->center()[0], face->line(l)->center()[1]};
          if(center_yz.norm() < 1.1 * r)
            face->line(l)->set_manifold_id(10);
        }
      }
    front.set_manifold(10, dealii::CylindricalManifold<dim>(2));
  }

  // transition top
  dealii::Triangulation<dim> top;
  {
    dealii::Point<dim> P1{0.0, r, 0};
    dealii::Point<dim> P2{end, h_ref, channel_width};
    auto const         refinements = balance_refinements(P1, P2, h_ref);
    dealii::GridGenerator::subdivided_hyper_rectangle(top, refinements, P1, P2);
  }
  // transition bottom
  dealii::Triangulation<dim> bottom;
  {
    dealii::Point<dim> P1{0.0, -r, 0};
    dealii::Point<dim> P2{end, -h_ref, channel_width};
    auto const         refinements = balance_refinements(P1, P2, h_ref);
    dealii::GridGenerator::subdivided_hyper_rectangle(bottom, refinements, P1, P2);
  }

  // wake
  dealii::Triangulation<dim> wake;
  {
    dealii::Point<dim> P1{plate_length, -r, 0};
    dealii::Point<dim> P2{end, r, channel_width};
    auto const         refinements = balance_refinements(P1, P2, h_ref);
    dealii::GridGenerator::subdivided_hyper_rectangle(wake, refinements, P1, P2);
  }

  dealii::GridGenerator::merge_triangulations({&front, &top, &bottom, &wake}, obstacle, 1e-6, true);

  // set boundary ids: wall=0, inlet=1, outlet=2,, symmetry=3, non-matching>90
  for(const auto & face : obstacle.active_face_iterators())
    if(face->at_boundary())
    {
      if(face->center()[0] > end - 1e-6)
        face->set_boundary_id(2);
      else if(face->center()[2] < 1e-6 || face->center()[2] > channel_width - 1e-6)
        face->set_boundary_id(3);
      else if(std::abs(face->center()[1]) < r + 1e-6)
        face->set_boundary_id(0);
      else
        face->set_boundary_id(99);
    }
}

template<int dim>
void
create_triangulation_surrounding(dealii::Triangulation<dim> & surrounding,
                                 double const                 channel_width)
{
  dealii::Triangulation<dim> tria;
  dealii::Point<dim>         p1{front, -0.5 * channel_height, 0};
  dealii::Point<dim>         p2{end, 0.5 * channel_height, channel_width};

  dealii::GridGenerator::subdivided_hyper_rectangle(tria,
                                                    balance_refinements(p1, p2, h_ref_outer),
                                                    p1,
                                                    p2);

  // remove cells
  using CellIter = typename dealii::Triangulation<dim>::active_cell_iterator;
  std::set<CellIter> cells_to_remove;
  for(unsigned int l = 0; l < coarsening_levels; ++l)
    for(auto const & cell : tria.active_cell_iterators())
      if(cell->is_locally_owned())
        if(cell->center()[0] > -h_ref && std::abs(cell->center()[1]) < h_ref)
          cells_to_remove.insert(cell);

  dealii::GridGenerator::create_triangulation_with_removed_cells(tria,
                                                                 cells_to_remove,
                                                                 surrounding);

  // set boundary ids: wall=0, inlet=1, outlet=2,, symmetry=3, non-matching>90
  for(const auto & face : surrounding.active_face_iterators())
    if(face->at_boundary())
    {
      if(face->center()[0] > end - 1e-6)
        face->set_boundary_id(2);
      else if(face->center()[0] < front + 1e-6)
        face->set_boundary_id(1);
      else if(std::abs(face->center()[1]) < h_ref + 1e-6 && face->center()[0] > -h_ref - 1e-6)
        face->set_boundary_id(98);
      else
        face->set_boundary_id(3);
    }
}



template<int dim>
void
create_triangulation_acoustic(dealii::Triangulation<dim> & tria_final,
                              double const                 channel_width,
                              unsigned int const           n_refinements_wake = 0)
{
  // outer mesh
  dealii::Triangulation<dim> surrounding;
  create_triangulation_surrounding(surrounding, channel_width);
  // wake
  dealii::Triangulation<dim> obstacle;
  create_obstacle(obstacle, channel_width);

  // final
  dealii::GridGenerator::merge_triangulations(
    {&obstacle, &surrounding}, tria_final, 0.0, true, true);
  tria_final.set_manifold(10, dealii::CylindricalManifold<dim>(2));

  // additional refinements in wake
  for(unsigned int l = 0; l < n_refinements_wake; ++l)
  {
    for(auto const & cell : tria_final.active_cell_iterators())
      if(cell->is_locally_owned())
        if(cell->center()[0] > plate_length - h_ref &&
           cell->center()[0] < plate_length + 2.0 * h_ref &&
           std::abs(cell->center()[1]) < 2.0 * h_ref)
          cell->set_refine_flag();

    tria_final.execute_coarsening_and_refinement();
  }
}

template<int dim>
void
create_triangulation_fluid(dealii::Triangulation<dim> & tria_final,
                           double const                 channel_width,
                           unsigned int const           n_refinements_global,
                           unsigned int const           n_refinements_boundary_layer = 0)
{
  // outer mesh
  dealii::Triangulation<dim> surrounding;
  dealii::Triangulation<dim> temp;
  create_triangulation_surrounding(temp, channel_width);
  temp.refine_global(coarsening_levels);
  dealii::GridGenerator::flatten_triangulation(temp, surrounding);

  // wake
  dealii::Triangulation<dim> obstacle;
  create_obstacle(obstacle, channel_width);

  // final
  dealii::GridGenerator::merge_triangulations({&obstacle, &surrounding}, tria_final, 1e-6, true);
  tria_final.set_manifold(10, dealii::CylindricalManifold<dim>(2));

  // set boundary ids: wall=0, inlet=1, outlet=2,, symmetry=3, non-matching>90
  for(const auto & face : tria_final.active_face_iterators())
    if(face->at_boundary())
    {
      if(face->center()[0] < front + 1e-6)
        face->set_boundary_id(1);
      else if(face->center()[0] > end - 1e-6)
        face->set_boundary_id(2);
      else if(std::abs(face->center()[1]) < r + 1e-6 && face->center()[0] < plate_length + 1e-6)
        face->set_boundary_id(0);
      else
        face->set_boundary_id(3);
    }

  // refine global
  tria_final.refine_global(n_refinements_global);

  // refine boundary layer
  for(uint i = 0; i < n_refinements_boundary_layer; ++i)
  {
    for(const auto & cell : tria_final.active_cell_iterators())
      for(const auto f : cell->face_indices())
        if(cell->face(f)->boundary_id() == 0)
          cell->set_refine_flag();

    tria_final.execute_coarsening_and_refinement();
  }
}


} // namespace Parker
} // namespace ExaDG
#endif /*APPLICATIONS_GRID_TOOLS_MESH_TESTCASE_PAUL_H_*/
