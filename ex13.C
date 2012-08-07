//export PETSC_ARCH=linux-gnu   (debug)
//export METHOD=dbg (debug libmesh)

#include "defines.h"
#include "assemble.h"
//#include "solid_system.h"
#include <iostream>
#include <fstream>
//#define ELEMENT_TYPE LAGRANGE 
//#define ORDER_HIGH SECOND
//#define ORDER_LOW FIRST 

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

// The main program.
int main (int argc, char** argv)
{

#if CUSTOM_MESH
  std::string mesh_file_name ("/ecslab/lorenzb/gmsh-2.5.0-Linux/lozo_annulus_coarse.msh");
//    std::string mesh_file_name ("pipethick.msh");
#endif

#if WRITE_MESH
std::string mesh_out_file_name ("parallel.msh");
#endif

#if CHAP
ofstream outFile;         // Step #2 - Declare outstream file.
outFile.open("results.txt");        // Step #3 - Open outFile.
#endif

//Some solver Options
const unsigned int n_nonlinear_steps = 25;
const Real nonlinear_tolerance       = 1.e-1;
const Real initial_linear_solver_tol = 1.e-8;

  // Initialize libMesh.
LibMeshInit init (argc, argv);

Mesh mesh(3);
MeshData mesh_data(mesh);

#if CUSTOM_MESH
 GmshIO(mesh).read(mesh_file_name);
#endif

#if CUBE
 MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1.5,	
                                       0., 1.5,
				                               0., 1.5,
                                       HEX27);
#endif
test(1);
mesh.all_second_order();
mesh.prepare_for_use();
// Print information about the mesh to the screen.
mesh.print_info();
EquationSystems equation_systems (mesh);  

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
 fluid_system_vel.attach_assemble_function (assemble_fluid_vel);
 TransientLinearImplicitSystem&  fluid_vel =
        equation_systems.get_system<TransientLinearImplicitSystem>("fluid-system-vel");
#endif

#if ASSEMBLE_PRESSURE
TransientLinearImplicitSystem & pressure_system = equation_systems.add_system<TransientLinearImplicitSystem> ("pressure-system");
pressure_system.add_variable ("fluid_pressure", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("jacobian", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("p_aux1", ORDER_HIGH,ELEMENT_TYPE);
 pressure_system.add_variable ("p_aux2", ORDER_LOW,ELEMENT_TYPE);
 pressure_system.attach_assemble_function (assemble_pressure);
 TransientLinearImplicitSystem&  pressure =
        equation_systems.get_system<TransientLinearImplicitSystem>("pressure-system");
#endif

#if ASSEMBLE_RESULTS
TransientLinearImplicitSystem & results_systems = equation_systems.add_system<TransientLinearImplicitSystem> ("pressure-system-lol");
results_systems.add_variable ("fluid_m", ORDER_HIGH,ELEMENT_TYPE);
 results_systems.add_variable ("Jacobian", ORDER_HIGH,ELEMENT_TYPE);
 results_systems.add_variable ("res_aux1",ORDER_HIGH,ELEMENT_TYPE);
 results_systems.add_variable ("res_aux2", ORDER_LOW,ELEMENT_TYPE);
 results_systems.attach_assemble_function (assemble_pressure_lol);
 TransientLinearImplicitSystem&  results_sys =
        equation_systems.get_system<TransientLinearImplicitSystem>("pressure-system-lol");
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

TransientLinearImplicitSystem & last_non_linear_soln_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");


#if SOLID_VELOCITY
  TransientLinearImplicitSystem & velocity_system = equation_systems.add_system<TransientLinearImplicitSystem> ("velocity-system");
#endif


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

#if LOG_PERFORMANCE
  // Create a performance-logging object for this example
  PerfLog perf_log("Solid Solver");
#endif   
// Now we begin the timestep loop to compute the time-accurate
// solution of the equations.
Real dt = 1;
Real time     = 0;
unsigned int n_timesteps = 1;

#if DYNAMIC
n_timesteps = 50;
dt = 0.01;
ExodusII_IO exo= ExodusII_IO(equation_systems.get_mesh());
#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif
Real dt_n_minus_1 = dt;
#endif 



#if PETSC_MUMPS
  PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
  PC pc = petsc_linear_solver->pc();
  int ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

equation_systems.parameters.set<Real> ("dt")   = dt;
equation_systems.init ();
equation_systems.print_info();
test(2);

// Get a reference to the systems to use later.
TransientLinearImplicitSystem&  newton_update =
        equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
	
TransientLinearImplicitSystem&  last_non_linear_soln =
        equation_systems.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");

TransientLinearImplicitSystem&  reference =
        equation_systems.get_system<TransientLinearImplicitSystem>("Reference-Configuration");

#if SOLID_VELOCITY
TransientLinearImplicitSystem&  velocity =
        equation_systems.get_system<TransientLinearImplicitSystem>("velocity-system");
#endif

#if UN_MINUS_ONE
TransientLinearImplicitSystem&  unm1 =
        equation_systems.get_system<TransientLinearImplicitSystem>("unm1-system");
#endif

// Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
const MeshBase::const_node_iterator nd_end =
      equation_systems.get_mesh().local_nodes_end();
  for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
    const Node *node = *nd; 
  test(4);
    for (unsigned int d = 0; d < 3; ++d) {

      unsigned int dest_dof = node->dof_number(reference.number(), d, 0);
      Real value = (*node)(d);
      reference.current_local_solution->set(dest_dof, value);
      reference.solution->set(dest_dof, value);
    }
}
test(5);
reference.solution->close();
reference.current_local_solution->close();
reference.update();   

AutoPtr<NumericVector<Number> > change_in_newton_update (newton_update.solution->clone());
AutoPtr<NumericVector<Number> > newton_solution (newton_update.solution->clone());
newton_solution->zero();

#if FLUID
AutoPtr<NumericVector<Number> > change_in_fluid (fluid_system.solution->clone());
#endif
#if FLUID_VEL
AutoPtr<NumericVector<Number> > change_in_fluid (fluid_system_vel.solution->clone());
#endif

// Load in the reference mesh as an initial guess (a neater way to do this ?)
last_non_linear_soln.solution->zero(); 
last_non_linear_soln.solution->add(1.,(*reference.solution));
last_non_linear_soln.solution->close();
last_non_linear_soln.current_local_solution->zero(); 
last_non_linear_soln.current_local_solution->add(1.,(*reference.current_local_solution));
last_non_linear_soln.current_local_solution->close();
last_non_linear_soln.old_local_solution->zero(); 
last_non_linear_soln.old_local_solution->add(1.,(*reference.current_local_solution));
last_non_linear_soln.old_local_solution->close();
last_non_linear_soln.update(); 

#if SOLID_VELOCITY
//Initial velocity is zero
  velocity.solution->zero(); 
  velocity.solution->close();
  velocity.current_local_solution->zero(); 
  velocity.current_local_solution->close();
  velocity.update(); 
#endif

#if UN_MINUS_ONE
unm1.solution->zero(); 
unm1.solution->add(1.,(*reference.solution));
unm1.solution->close();
unm1.current_local_solution->zero(); 
unm1.current_local_solution->add(1.,(*reference.current_local_solution));
unm1.current_local_solution->close();
unm1.old_local_solution->zero(); 
unm1.old_local_solution->add(1.,(*reference.current_local_solution));
unm1.old_local_solution->close();
unm1.update(); 
#endif

#if FLUID
fluid_system.solution->zero(); 
fluid_system.solution->close();
#endif

#if FLUID_VEL
fluid_system_vel.solution->zero(); 
fluid_system_vel.solution->close();
#endif


  for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
    {
      equation_systems.parameters.set<Real> ("time") = time;
      double progress = (t_step+0.00000000001) / (n_timesteps+0.00000000001);
      equation_systems.parameters.set<Real>("progress") = progress;
      equation_systems.parameters.set<unsigned int>("step") = t_step; 

#if ADAPTIVE_TIME
if (progress>0.09){ equation_systems.parameters.set<Real> ("dt") = 0.05; }
if (progress>0.2){  equation_systems.parameters.set<Real> ("dt") = 0.1;  }
if (progress>0.5){  equation_systems.parameters.set<Real> ("dt") = 0.1;  }
#endif

      time += dt;

      std::cout << "\n\n*** Solving time step " << t_step << ", time = " << time <<  ", progress = " << progress << " ***" << std::endl;

// Now we need to update the solution vector from the
// previous time step.  This is done directly through
// the reference to the system.
#if DYNAMIC
     *last_non_linear_soln.old_local_solution = *last_non_linear_soln.current_local_solution;
#endif

#if DYNAMIC && FLUID
     *fluid_system.old_local_solution = *fluid_system.current_local_solution;
#endif

#if ASSEMBLE_PRESSURE
     *pressure_system.old_local_solution = *pressure_system.current_local_solution;
#endif
#if ASSEMBLE_PRESSURE_GRAD
     *pressure_gard_system.old_local_solution = *pressure_gard_system.current_local_solution;
#endif
    
 #if UN_MINUS_ONE 
//*unm1.old_local_solution == U_(n-1)
    *unm1.old_local_solution = *unm1.current_local_solution; 
    *unm1.current_local_solution=*last_non_linear_soln.old_local_solution;
#endif   

// Now we begin the nonlinear loop
for (unsigned int l=0; l<n_nonlinear_steps; ++l)
        {
      equation_systems.parameters.set<Real> ("non_lin_step") = l;

      	  std::cout<<"\n Non-linear iteration " << l << std::endl;

//Makes sure linear solver tolearance hasn't changed
//  equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

change_in_newton_update->zero();
change_in_newton_update->add(*newton_update.solution);


#if VERIFY_JACK

//verify_jack(equation_systems);

#endif

//Prepare the newton update system for it's linear solve
newton_update.solution->zero();  
newton_update.update();

// Assemble & solve the linear system.
#if SOLVE_SOLID

#if !PETSC_MUMPS
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 2000;
#endif
#if PETSC_MUMPS
 petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
   pc = petsc_linear_solver->pc();
   ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

clock_t begin_solid_solve=clock();
equation_systems.get_system("Newton-update").solve();
clock_t end_solid_solve=clock();
equation_systems.reinit();


#endif 
	//update the final solution
	// xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
	//Apply a full Newton-step
Real K=1; //Newton step size
	last_non_linear_soln.solution->add(-K*1,*newton_update.solution);
	last_non_linear_soln.solution->close();
 test(2);
  	last_non_linear_soln.current_local_solution->add(-K*1,*newton_update.current_local_solution);
  	last_non_linear_soln.current_local_solution->close();
  	last_non_linear_soln.update();
test(3);
//std::cout<<" last_non_linear_soln " <<  (*last_non_linear_soln.current_local_solution)(123)<<std::endl;

#if SOLID_VELOCITY //This is to compute Vn=(Un+1^{k}-Un)/dtn-1 !! is this strange ?
//Calculate the new velocity using a classic finite difference 

//Should be dt_n_minus_1 instead of dt
  velocity.solution->zero(); 
  velocity.solution->add(1.0/dt_n_minus_1,*last_non_linear_soln.current_local_solution);
  velocity.solution->add(-1.0/dt_n_minus_1,*last_non_linear_soln.old_local_solution);
  velocity.solution->close();
  velocity.current_local_solution->zero(); 
  velocity.current_local_solution->add(1.0/dt_n_minus_1,*last_non_linear_soln.current_local_solution);
  velocity.current_local_solution->add(-1.0/dt_n_minus_1,*last_non_linear_soln.old_local_solution);
  velocity.current_local_solution->close();
  velocity.update();

  std::cout<<"dt " <<  dt <<std::endl;


  std::cout<<" last_non_linear_soln current " <<  (*last_non_linear_soln.current_local_solution)(123)<<std::endl;
  std::cout<<" last_non_linear_soln old " <<  (*last_non_linear_soln.old_local_solution)(123)<<std::endl;

  std::cout<<" velocity " <<  (*velocity.current_local_solution)(123)<<std::endl;

  dt_n_minus_1=dt;
#endif
test(4);

#if MOVING_MESH
//Don't do mesh updating anymore
//Update the mesh - maybe can do this more elegantly by using built in functions or writing my own
//This culd be reimplemented using the buildt inFEMSystem class and usind set_mesh_x_var etc 
 //Have to use an element iterator to make sure we deal with ghost nodes correctly.
MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 
  for ( ; el != end_el; ++el)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
  for (unsigned int n=0; n<elem->n_nodes(); n++){
       Node *node = elem->get_node(n);
  for (unsigned int d = 0; d < 3; ++d) {
        unsigned int source_dof = node->dof_number(1, d, 0);
        Real value = last_non_linear_soln.current_local_solution->el(source_dof);
  (*node)(d)=value;
        }
  }
}
#endif

//Find some interesting dofs to inspect
Node node = mesh.node(10);
//std::cout<<" node a " <<  node  <<std::endl;

unsigned int u_dof_a = node.dof_number(1, 0, 0);
unsigned int v_dof_a = node.dof_number(1, 1, 0);
unsigned int w_dof_a = node.dof_number(1, 2, 0);
unsigned int p_dof_a = node.dof_number(1, 3, 0);

node = mesh.node(100);
//std::cout<<" node b " <<  node  <<std::endl;

unsigned int u_dof_b = node.dof_number(1, 0, 0);
unsigned int v_dof_b = node.dof_number(1, 1, 0);
unsigned int w_dof_b = node.dof_number(1, 2, 0);
unsigned int p_dof_b = node.dof_number(1, 3, 0);

node = mesh.node(150);
//std::cout<<" node c " <<  node  <<std::endl;

unsigned int u_dof_c = node.dof_number(1, 0, 0);
unsigned int v_dof_c = node.dof_number(1, 1, 0);
unsigned int w_dof_c = node.dof_number(1, 2, 0);
unsigned int p_dof_c = node.dof_number(1, 3, 0);


//Assemble the pressure of the porouse medium 'and solve'
#if ASSEMBLE_PRESSURE
equation_systems.get_system("pressure-system").solve();
std::cout<<" pressure a " <<  (*pressure.current_local_solution)(u_dof_a)<<std::endl;
#endif



#if ASSEMBLE_PRESSURE_GRAD
equation_systems.get_system("pressure-grad-system").solve();
std::cout<<" pressure_grad 123 " <<  (*pressure_grad.current_local_solution)(123)<<std::endl;
#endif



change_in_newton_update->add (-1., *newton_update.solution);
change_in_newton_update->close();
Real norm_delta = change_in_newton_update->l2_norm();

Real norm_delta_fluid = 999.999;

 if ((norm_delta < 10.1) ){

  std::cout<<" Solving fluid " <<  std::endl;

//----------------- Solve the Fluid -------------------------------///
#if FLUID
//Save previous solution so we can compare against teh new solution
change_in_fluid->zero();
change_in_fluid->add(*fluid_system.solution);

fluid.solution->zero();  
fluid.update();

#if PETSC_MUMPS
petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
   pc = petsc_linear_solver->pc();
   ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

clock_t begin_fluid_solve=clock();
equation_systems.get_system("fluid-system").solve();
clock_t end_fluid_solve=clock();


equation_systems.reinit();
#endif


#if FLUID_VEL
//Save previous solution so we can compare against teh new solution
change_in_fluid->zero();
change_in_fluid->add(*fluid_system_vel.solution);

fluid_vel.solution->zero();  
fluid_vel.update();

#if PETSC_MUMPS
petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
   pc = petsc_linear_solver->pc();
   ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

clock_t begin_fluid_solve=clock();
equation_systems.get_system("fluid-system-vel").solve();
clock_t end_fluid_solve=clock();

equation_systems.reinit();
#endif
// End of Fluid solve ------------------------

#if ASSEMBLE_RESULTS
equation_systems.get_system("pressure-system-lol").solve();
#endif

#if PRINT_NON_LINEAR
 last_non_linear_soln.update();      
 newton_update.update(); 
 print_resluts(equation_systems,time,l);
#endif

  #if FLUID || FLUID_VEL
          #if FLUID
          change_in_fluid->add (-1., *fluid_system.solution);
          #endif
           #if FLUID_VEL
          change_in_fluid->add (-1., *fluid_system_vel.solution);
          #endif
          change_in_fluid->close();
           norm_delta_fluid = change_in_fluid->l2_norm();
          #endif 

}else{l=l-1;}

          /*****///Convergence and Norm Computing Stuff///***********/
          change_in_newton_update->add (-1., *newton_update.solution);
          change_in_newton_update->close();
           norm_delta = change_in_newton_update->l2_norm();

        

          // How many iterations were required to solve the linear system?
          const unsigned int n_linear_iterations = newton_update.n_linear_iterations();
          
          // What was the final residual of the linear system?
          const Real final_linear_residual = newton_update.final_linear_residual();
          	  
          #if SOLVE_SOLID
          std::cout << "Solid: Linear conv at step: "
                    << n_linear_iterations
                    << ", resid: "
                    << final_linear_residual 
                    << ", time: " << double(diffclock(end_solid_solve,begin_solid_solve)) << " ms"
                    <<std::endl;
          #endif

          #if FLUID || FLUID_VEL
          #if FLUID
          const unsigned int n_linear_iterations_fluid = fluid_system.n_linear_iterations();
          const Real final_linear_residual_fluid = fluid_system.final_linear_residual();  
          #endif
          #if FLUID_VEL
          const unsigned int n_linear_iterations_fluid = fluid_system_vel.n_linear_iterations();
          const Real final_linear_residual_fluid = fluid_system_vel.final_linear_residual();  
          #endif
          std::cout << "Fluid: Linear conv at step: "
                    << n_linear_iterations_fluid
                    << ", resid: "
                    << final_linear_residual_fluid 
                   // << ", time: " << double(diffclock(end_fluid_solve,begin_fluid_solve)) << " ms"
                    <<std::endl;
          #endif
                    
          std::cout   << "Solid Nonlinear convergence: ||u - u_old|| = "
                    << norm_delta
                    << std::endl;
          #if FLUID || FLUID_VEL

          std::cout   << "Fluid convergence: ||w - w_old|| = "
                    << norm_delta_fluid
                    << std::endl;
          #endif



#if (FLUID || FLUID_VEL) && SOLVE_SOLID
 if ((norm_delta < nonlinear_tolerance) && (norm_delta_fluid < nonlinear_tolerance)  ){
    std::cout << " Nonlinear solver converged at step "<<std::endl; 
    break;
  } 
#endif

#if ! (FLUID||FLUID_VEL) &&  SOLVE_SOLID
 if ((norm_delta < nonlinear_tolerance) ){
    std::cout << "Nonlinear solver converged at step "<<std::endl; break;} 
#endif


        } // end nonlinear loop

#if CHAP && ASSEMBLE_RESULTS
Point p1(0,0,0);
Point p2(0.75,0.75,0.75);
Point p3(1.5,1.5,1.5);

outFile << fluid_system_vel.point_value(2,p1) << " ";
outFile<<fluid_system_vel.point_value(2,p2)<< " ";
outFile<<fluid_system_vel.point_value(2,p3)<< " ";

outFile << results_systems.point_value(1,p1) << " ";
outFile<<results_systems.point_value(1,p2)<< " ";
outFile<<results_systems.point_value(1,p3)<< " ";

outFile << results_systems.point_value(2,p1) << " ";
outFile<<results_systems.point_value(2,p2)<< " ";
outFile<<results_systems.point_value(2,p3)<< " ";

outFile<< endl;   
#endif

 last_non_linear_soln.update();      
 newton_update.update(); 

 // Write out every nth timestep to file.
 const unsigned int write_interval = 1;
      
#ifdef LIBMESH_HAVE_EXODUS_API
      if ((t_step+1)%write_interval == 0)
        { 
  #if STATIC  
        OStringStream file_name;
        file_name << "/ecslab/lorenzb/Dropbox/libresults/cube_with_stokes";
        OSSRealzeroright(file_name,3,0, t_step + 1);
        file_name << ".e";	  
	ExodusII_IO(mesh).write_equation_systems (file_name.str(),equation_systems);

  #if WRITE_MESH
	GmshIO(mesh).write(mesh_out_file_name);
  #endif
  #endif  

  #if DYNAMIC
 
  std::stringstream file_name;
  #if ! DAHOAM 
  file_name << "/ecslab/lorenzb/Dropbox/libresults/";
  #endif
  #if DAHOAM 
  file_name << "/home/lorenz/Dropbox/libresults/";
  #endif
  file_name << "poro_newtest_";
  file_name << std::setw(2) << std::setfill('0') << t_step;
  file_name << ".e-s.";
  file_name << std::setw(3) << std::setfill('0') << t_step+1;
  exo.write_timestep(file_name.str(), equation_systems,t_step+1,time);

#if WRITE_TEC
  std::stringstream file_name_tec;
  file_name_tec << "cube_darcy"<< t_step<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
#endif
#endif
     }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    } // end timestep loop.
    // All done.  

#if CHAP
outFile.close();      // Step #5 - Close file.
#endif

  return 0;
}
