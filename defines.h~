#ifndef DEFINES_H_
#define DEFINES_H_

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <petsc_linear_solver.h>

// Basic include file needed for the libmesh functionality.
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
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "perf_log.h"
#include "o_string_stream.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "quadrature.h"
#include "elem.h"
#include "math.h"
#include "gmsh_io.h"
#include "mesh_data.h"
#include <tecplot_io.h>
//New header files
#include "solid_system.h"
#include "functions.h"

#include <time.h>

// Bring in everything from the libMesh namespace
using namespace libMesh;
using namespace std;

#define HERM 0


///------------------------------//



#define PORO 1

#define VERIFY_TEST_FINAL 0
#define CHAP 1
#define COMPRESSIBLE 0
#define INCOMPRESSIBLE 1
#define INCOMPRESSIBLE_STRESS 0 //I think the linearization is not quite correct especially the fourth order tensor.
#define INCOMPRESSIBLE_CHEAT 0

#define SOLVE_SOLID 1

#define FLUID 0
#define FLUID_VEL 1

#define ASSEMBLE_PRESSURE 0
#define ASSEMBLE_PRESSURE_GRAD 0
#define ASSEMBLE_RESULTS 1


#define LOG_PERFORMANCE 0
#define VERIFY_JACK 0
#define PETSC_MUMPS 1

#define PRINT_NON_LINEAR 0

#define CUBE 1
#define BEAM 0
#define CUSTOM_MESH 0

#define STATIC 0
#define DYNAMIC 1
#define SOLID_VELOCITY 1


#define FIXED_MESH 0
#define MOVING_MESH 1

#define TEST_OUTPUT 0


//Main
#define ADAPTIVE_TIME 0
#define WRITE_TEC 0
#define WRITE_MESH 0

#define UN_MINUS_ONE 0

//Most stable seems to be vn = U_n+1^{k} - un / dt

//Assemble
#define GRAVITY 0
#define LOG_ASSEMBLE_PERFORMANCE 0
#define ELASTICITY 1

//Boundary conditions
#define TRACTION_BC 0
#define DIRICHLET_CLASSIC 1
#define PENALTY 0
#define MOVING_DIRICHLET_BCS 0

#define CHAP_SWELL 1


//Fluid BCs
#define NEUMANN_PRESSURE	1
//Different Neumann pressure bcs
#define MOMENTUM_NEUMANN_PRESSURE 1
#define MASS_NEUMANN_PRESSURE 0

#define DIRICHLET_PRESSURE	1
#define DIRICHLET_VELOCITY	1
#define DISK_FLOW	0	
#define SEALED_CUBE	0	


#endif /* DEFINES_H_ */

