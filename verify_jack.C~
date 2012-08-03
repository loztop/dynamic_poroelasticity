   // C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "transient_system.h"
#include "perf_log.h"
#include "boundary_info.h"
#include "utility.h"
#include "elem.h"
#include "mesh_data.h"
#include "gmsh_io.h"

// Some (older) compilers do not offer full stream 
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"


#include "assemble.h"
#include "solid_system.h"
#include <petsc_linear_solver.h>
using namespace std;

void verify_jack(EquationSystems& es) {
TransientLinearImplicitSystem & newton_update =
   es.get_system<TransientLinearImplicitSystem> ("Newton-update");

unsigned int N=50;
double h_jac=0.02;
double rhs_orig;
double rhs_new;
double anal_jac;
double num_jac;
es.get_system("Newton-update").assemble();
 MeshBase::const_node_iterator nd_end =
      es.get_mesh().nodes_end();
  for (MeshBase::const_node_iterator nd = es.get_mesh().nodes_begin();
      nd != nd_end; ++nd) {
    Node *node = *nd;
for (unsigned int d = 0; d < 3; ++d) {
      unsigned int source_dof = node->dof_number(1, d, 0);
rhs_orig=newton_update.rhs->el(source_dof);
//std::cout << "rhs_orig "<< rhs_orig <<std::endl;
(*node)(d)=(*node)(d)+h_jac;
es.get_system("Newton-update").assemble();
rhs_new=newton_update.rhs->el(source_dof);
//std::cout << "rhs_new "<< rhs_new <<std::endl;

//Restore the node to its original position
(*node)(d)=(*node)(d)-h_jac;
anal_jac=newton_update.matrix->operator() (source_dof,source_dof);
std::cout << "anal_jac: "<<anal_jac<<std::endl;
num_jac=(rhs_new-rhs_orig)/h_jac;
std::cout << "num_jac: "<<num_jac<<std::endl;
   }
  }

         
/*
//inspect a particular element
const MeshBase& inspect_mesh = es.get_mesh();
 Elem* inspect_el=inspect_mesh.elem(5);
inspect_el->print_info(); */

}
