/*
 * BoundaryDescriptorLaplace.h
 *
 *  Created on: Nov 8, 2016
 *      Author: fehn
 */

#ifndef INCLUDE_BOUNDARYDESCRIPTORLAPLACE_H_
#define INCLUDE_BOUNDARYDESCRIPTORLAPLACE_H_

using namespace dealii;

#include <deal.II/base/function.h>
#include <deal.II/base/types.h>

template<int dim>
struct BoundaryDescriptorLaplace
{
  std::map<types::boundary_id,std_cxx11::shared_ptr<Function<dim> > > dirichlet;
  std::map<types::boundary_id,std_cxx11::shared_ptr<Function<dim> > > neumann;
};


#endif /* INCLUDE_BOUNDARYDESCRIPTORLAPLACE_H_ */