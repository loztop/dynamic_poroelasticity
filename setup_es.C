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
#include "defines.h"

//#include "solid_system.h"
#include <petsc_linear_solver.h>
using namespace std;


#define ELEMENT_TYPE LAGRANGE 
#define ORDER_HIGH SECOND
#define ORDER_LOW FIRST 


// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_solid (EquationSystems& es,
                      const std::string& system_name);
void assemble_fluid (EquationSystems& es,
                      const std::string& system_name);
#if FLUID_VEL
void assemble_fluid_vel (EquationSystems& es,
                      const std::string& system_name);
#endif

#if ASSEMBLE_PRESSURE
void assemble_pressure (EquationSystems& es,
                      const std::string& system_name);
#endif
#if ASSEMBLE_RESULTS
void assemble_pressure_lol (EquationSystems& es,
                      const std::string& system_name);
#endif

#if ASSEMBLE_PRESSURE_GRAD
void assemble_pressure_grad (EquationSystems& es,
                      const std::string& system_name);
#endif

void setup_es(EquationSystems& equation_systems) {
 TransientLinearImplicitSystem & ref_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Reference-Configuration");
TransientLinearImplicitSystem & system = equation_systems.add_system<TransientLinearImplicitSystem> ("Newton-update");

#if FLUID
 TransientLinearImplicitSystem & fluid_system = equation_systems.add_system<TransientLinearImplicitSystem> ("fluid-system");
 fluid_system.add_variable ("fluid_m", ORDER_HIGH,ELEMENT_TYPE);
 fluid_system.add_variable ("fluid_aux1",  ORDER_HIGH,ELEMENT_TYPE);
 fluid_system.add_variable ("fluid_aux2",  ORDER_HIGH,ELEMENT_TYPE);
 fluid_system.add_variable ("fluid_aux3",  ORDER_LOW,ELEMENT_TYPE);
 fluid_system.attach_assemble_function (assemble_fluid);
 TransientLinearImplicitSystem&  fluid =
        equation_systems.get_system<TransientLinearImplicitSystem>("fluid-system");
#endif

  #if FLUID_VEL
 TransientLinearImplicitSystem & fluid_system_vel = equation_systems.add_system<TransientLinearImplicitSystem> ("fluid-system-vel");
 fluid_system_vel.add_variable ("fluid_U_vel", ORDER_HIGH,ELEMENT_TYPE);
 fluid_system_vel.add_variable ("fluid_V_vel", ORDER_HIGH,ELEMENT_TYPE);
 fluid_system_vel.add_variable ("fluid_W_vel", ORDER_HIGH,ELEMENT_TYPE);
 fluid_system_vel.add_variable ("fluid_P", ORDER_LOW,ELEMENT_TYPE);
#if FLUID_P_CONST
fluid_system_vel.add_variable ("fluid_M", ORDER_LOW,ELEMENT_TYPE);
#endif
 fluid_system_vel.attach_assemble_function (assemble_fluid_vel);
#endif


#if ASSEMBLE_PRESSURE
TransientLinearImplicitSystem & pressure_system = equation_systems.add_system<TransientLinearImplicitSystem> ("pressure-system");
pressure_system.add_variable ("fluid_pressure", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("jacobian", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("p_aux1", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("p_aux2", ORDER_LOW,ELEMENT_TYPE);
 pressure_system.attach_assemble_function (assemble_pressure);
#endif

#if ASSEMBLE_RESULTS
TransientLinearImplicitSystem & results_systems = equation_systems.add_system<TransientLinearImplicitSystem> ("pressure-system-lol");
results_systems.add_variable ("fluid_m", ORDER_HIGH,ELEMENT_TYPE);
 results_systems.add_variable ("Jacobian", ORDER_HIGH,ELEMENT_TYPE);
 results_systems.add_variable ("res_aux1",ORDER_HIGH,ELEMENT_TYPE); 
  results_systems.add_variable ("res_aux2", ORDER_LOW,ELEMENT_TYPE);
 results_systems.attach_assemble_function (assemble_pressure_lol);
#endif

#if ASSEMBLE_PRESSURE_GRAD
TransientLinearImplicitSystem & pressure_gard_system = equation_systems.add_system<TransientLinearImplicitSystem> ("pressure-grad-system");
pressure_gard_system.add_variable ("p_grad_u", ORDER_HIGH,ELEMENT_TYPE);
 pressure_gard_system.add_variable ("p_grad_v", ORDER_HIGH,ELEMENT_TYPE);
 pressure_gard_system.add_variable ("p_grad_w", ORDER_HIGH,ELEMENT_TYPE);
 pressure_gard_system.add_variable ("p_grad_aux", ORDER_LOW,ELEMENT_TYPE);
 pressure_gard_system.attach_assemble_function (assemble_pressure_grad);
 TransientLinearImplicitSystem&  pressure_grad =
        equation_systems.get_system<TransientLinearImplicitSystem>("pressure-grad-system");
#endif

#if ANALNEO
TransientLinearImplicitSystem & anal_system = equation_systems.add_system<TransientLinearImplicitSystem> ("anal-system");
anal_system.add_variable ("anal_u", ORDER_HIGH,ELEMENT_TYPE);
 anal_system.add_variable ("anal_v", ORDER_HIGH,ELEMENT_TYPE);
 anal_system.add_variable ("anal_w", ORDER_HIGH,ELEMENT_TYPE);
 anal_system.add_variable ("anal_aux", ORDER_LOW,ELEMENT_TYPE);
// pressure_gard_system.attach_assemble_function (assemble_pressure_grad);
 TransientLinearImplicitSystem&  anal =
        equation_systems.get_system<TransientLinearImplicitSystem>("anal-system");
#endif

TransientLinearImplicitSystem & last_non_linear_soln_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");



#if UN_MINUS_ONE
 TransientLinearImplicitSystem & unm1_system = equation_systems.add_system<TransientLinearImplicitSystem> ("unm1-system");
#endif

system.add_variable ("u_nu", ORDER_HIGH,ELEMENT_TYPE);
system.add_variable ("v_nu", ORDER_HIGH,ELEMENT_TYPE);
system.add_variable ("w_nu", ORDER_HIGH,ELEMENT_TYPE);
system.attach_assemble_function (assemble_solid);

last_non_linear_soln_system.add_variable ("u", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("v", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("w", ORDER_HIGH,ELEMENT_TYPE);

ref_system.add_variable ("u_ref", ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("v_ref", ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("w_ref", ORDER_HIGH,ELEMENT_TYPE);

#if SOLID_VELOCITY
  velocity_system.add_variable ("u_vel", ORDER_HIGH,ELEMENT_TYPE);
  velocity_system.add_variable ("v_vel", ORDER_HIGH,ELEMENT_TYPE);
  velocity_system.add_variable ("w_vel", ORDER_HIGH,ELEMENT_TYPE);
  #if INCOMPRESSIBLE  
  velocity_system.add_variable ("craP",ORDER_LOW,ELEMENT_TYPE);
#endif
#endif

#if UN_MINUS_ONE
  unm1_system.add_variable ("u_nm1", ORDER_HIGH,ELEMENT_TYPE);
  unm1_system.add_variable ("v_nm1", ORDER_HIGH,ELEMENT_TYPE);
  unm1_system.add_variable ("w_nm1", ORDER_HIGH,ELEMENT_TYPE);
  #if INCOMPRESSIBLE
  unm1_system.add_variable ("craP", ORDER_LOW,ELEMENT_TYPE);
  #endif
#endif

#if INCOMPRESSIBLE
  system.add_variable ("p_nu", ORDER_LOW,ELEMENT_TYPE);
  last_non_linear_soln_system.add_variable ("p", ORDER_LOW,ELEMENT_TYPE);
  ref_system.add_variable ("p_ref",ORDER_LOW,ELEMENT_TYPE);
#endif  





}
