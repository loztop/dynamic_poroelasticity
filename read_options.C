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
//#include "defines.h"
#include <petsc_linear_solver.h>
using namespace std;


void read_options(unsigned int &  n_timesteps, std::string& result_file_name,const char*& plot_out_file_name, int& argc, char**& argv){




n_timesteps= atoi( argv[1] );
result_file_name =  argv[2] ;
plot_out_file_name =  argv[3] ;


    std::cout<<"n_timesteps "<< n_timesteps <<" \n";
    std::cout<<"result_file_name "<< result_file_name <<" \n";
    std::cout<<"plot_out_file_name "<< plot_out_file_name <<" \n";


}
