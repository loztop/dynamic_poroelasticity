//export PETSC_ARCH=linux-gnu   (debug)
//export METHOD=dbg (debug libmesh)

#include "defines.h"
#include "assemble.h"
#include <iostream>
#include <fstream>

#define WRITE_RESULTS 1

// The main program.
int main (int argc, char** argv)
{

LibMeshInit init (argc, argv);

#include "setup_mesh.txt"

EquationSystems equation_systems (mesh);  
setup_es(equation_systems);
Real time     = 0;
Real end_time     = 1.5;

const unsigned int n_nonlinear_steps = 25;
const Real nonlinear_tolerance       = 1.e-2;
const Real initial_linear_solver_tol = 1.e-18;
unsigned int n_timesteps = 100;
Real dt = end_time/n_timesteps;


ExodusII_IO exo= ExodusII_IO(equation_systems.get_mesh());
#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif
Real dt_n_minus_1 = dt;

std::cout<<"dt = "<<dt<<std::endl;

const char* plot_out_file_name ="plot.txt";
std::string result_file_name ("SquashSponge2");

if (argc >2){
read_options(n_timesteps,result_file_name,plot_out_file_name,argc, argv);
}

dt = end_time/n_timesteps;

//dt=1;

ofstream outFile; // Step #2 - Declare outstream file.
outFile.open(plot_out_file_name, std::ios_base::app);        // Step #3 - Open outFile.
//outFile <<"\n" << result_file_name  << " " << n_timesteps << "\n";


//1000*(iHat*fluid_U_vel + jHat*fluid_V_vel+ kHat*fluid_W_vel)
//iHat*fluid_P + jHat*fluid_P+ kHat*fluid_P

//10000*(iHat*fluid_U_vel + jHat*fluid_V_vel+ kHat*fluid_W_vel)+1*(iHat*a + jHat*b+ kHat*c)

#if LOG_PERFORMANCE
  PerfLog perf_log("Solid Solver");
#endif   

equation_systems.parameters.set<Real> ("dt")   = dt;
equation_systems.init ();
//equation_systems.print_info();

// Get a reference to the systems to use later.
#include "get_references_to_systems.txt"

#if PETSC_MUMPS
  PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
  PC pc = petsc_linear_solver->pc();
  int ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

//Find nodes for result output later;
    const Node *node_0 ; 
    const Node *node_1 ; 
    const Node *node_2 ; 

// Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
const MeshBase::const_node_iterator nd_end =
      equation_systems.get_mesh().local_nodes_end();
  for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
    const Node *node = *nd; 


  if ( abs((*node)(0)-0) + abs((*node)(1)-0) + abs((*node)(2)-0) < 0.001){
     node_0 = *nd; 
    std::cout<<"node "<< (*node_0)(0)<<std::endl;
  }
  if ( abs((*node)(0)-1) + abs((*node)(1)-0) + abs((*node)(2)-0) < 0.001){
      node_1 = *nd; 
    std::cout<<"node1 "<< (*node_1)(0)<<std::endl;
  }
  if ( abs((*node)(0)-0) + abs((*node)(1)-1) + abs((*node)(2)-1) < 0.001){
      node_2 = *nd; 
    std::cout<<"node2 "<< (*node_2)(0)<<std::endl;
  }


//Copy initial mesh into reference.solution
    for (unsigned int d = 0; d < 3; ++d) {

      unsigned int dest_dof = node->dof_number(reference.number(), d, 0);
      Real value = (*node)(d);
      reference.current_local_solution->set(dest_dof, value);
      reference.solution->set(dest_dof, value);
    }


}

reference.solution->close();
reference.current_local_solution->close();
reference.update();   

AutoPtr<NumericVector<Number> > change_in_newton_update (newton_update.solution->clone());
AutoPtr<NumericVector<Number> > newton_solution (newton_update.solution->clone());
newton_solution->zero();

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

#if FLUID_VEL
fluid_system_vel.solution->zero(); 
fluid_system_vel.solution->close();
#endif

  for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
    {
      time += dt;
      equation_systems.parameters.set<Real> ("time") = time;
      double progress = (t_step+0.000000001) / (n_timesteps+0.000000001);
      equation_systems.parameters.set<Real>("progress") = progress;
      equation_systems.parameters.set<unsigned int>("step") = t_step; 

//#include "adaptive_time_stepping.txt"

      std::cout << "\n\n*** Solving time step " << t_step << ", time = " << time <<  ", progress = " << progress << " ***" << std::endl;

// Now we need to update the solution vector from the
// previous time step.  This is done directly through
// the reference to the system.
#if DYNAMIC
     *last_non_linear_soln.old_local_solution = *last_non_linear_soln.current_local_solution;
#endif

#if DYNAMIC && FLUID_VEL
     *fluid_system_vel.old_local_solution = *fluid_system_vel.current_local_solution;
#endif

#if ASSEMBLE_PRESSURE
     *pressure_system.old_local_solution = *pressure_system.current_local_solution;
#endif
#if ASSEMBLE_PRESSURE_GRAD
     *pressure_gard_system.old_local_solution = *pressure_gard_system.current_local_solution;
#endif
    
 #if UN_MINUS_ONE 
//*unm1.old_local_solution == U_(n-1)
    *unm1.old_local_solution = *unm1.current_local_solution;  *unm1.current_local_solution=*last_non_linear_soln.old_local_solution;
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
#if !PETSC_MUMPS
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 2000;
equation_systems.parameters.set<Real>("linear solver minimum tolerance") = 1e-18;
#endif

#if PETSC_MUMPS
  PetscLinearSolver<Number>* petsc_linear_solver_solid =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
  PC pc_solid = petsc_linear_solver_solid->pc();
  int ierr_solid = PCSetType(pc_solid, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr_solid);
  ierr_solid = PCFactorSetMatSolverPackage(pc_solid,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr_solid);
#endif

clock_t begin_solid_solve=clock();
equation_systems.get_system("Newton-update").solve();
clock_t end_solid_solve=clock();
equation_systems.reinit();
test(1);  
//update the final solution
	// xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
	//Apply a full Newton-step
Real K=1; //Newton step size
	last_non_linear_soln.solution->add(-1*K,*newton_update.solution);
	last_non_linear_soln.solution->close();
 test(2);
  	last_non_linear_soln.current_local_solution->add(-1,*newton_update.current_local_solution);
  	last_non_linear_soln.current_local_solution->close();
  	last_non_linear_soln.update();
test(3);
//std::cout<<" last_non_linear_soln " <<  (*last_non_linear_soln.current_local_solution)(123)<<std::endl;

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

//Assemble the pressure of the porouse medium 'and solve'
#if ASSEMBLE_PRESSURE
equation_systems.get_system("pressure-system").solve();
#endif
#if ASSEMBLE_PRESSURE_GRAD
equation_systems.get_system("pressure-grad-system").solve();
#endif

change_in_newton_update->add (-1., *newton_update.solution);
change_in_newton_update->close();
Real norm_delta = change_in_newton_update->l2_norm();

Real norm_delta_fluid = 999.999;

 if ((norm_delta < 120.000001) ){

 // std::cout<<" Solving fluid " <<  std::endl;

//------ Solve the Fluid--------------///
#if FLUID
//Save previous solution so we can compare against teh new solution
change_in_fluid->zero();
change_in_fluid->add(*fluid_system.solution);

fluid.solution->zero();  
fluid.update();

#if !PETSC_MUMPS
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 2000;
equation_systems.parameters.set<Real>("linear solver tolerance") = 1e-20;
#endif

#if PETSC_MUMPS
PetscLinearSolver<Number>* petsc_linear_solver_fluid =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
  PC pc_solid = petsc_linear_solver_solid->pc();
  int ierr_solid = PCSetType(pc_solid, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr_solid);
  ierr_solid = PCFactorSetMatSolverPackage(pc_solid,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr_solid);
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
fluid_system_vel.solution->zero();  
fluid_system_vel.update();



#if PETSC_MUMPS
PetscLinearSolver<Number>* petsc_linear_solver_fluid =dynamic_cast<PetscLinearSolver<Number>*>(fluid_system_vel.get_linear_solver());
  PC pc_fluid = petsc_linear_solver_fluid->pc();
  int ierr_fluid = PCSetType(pc_fluid, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr_fluid);
  ierr_fluid = PCFactorSetMatSolverPackage(pc_fluid,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr_fluid);
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

          #if  FLUID_VEL
          const unsigned int n_linear_iterations_fluid = fluid_system_vel.n_linear_iterations();
          const Real final_linear_residual_fluid = fluid_system_vel.final_linear_residual();  
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


#if ( FLUID_VEL) && SOLVE_SOLID
 if ((norm_delta < nonlinear_tolerance) && (norm_delta_fluid < nonlinear_tolerance)  ){
 //if ((norm_delta < nonlinear_tolerance) ){
    std::cout << "Nonlinear solver converged at step "<<std::endl;
    break;
  } 
#endif

#if ! (FLUID_VEL) &&  SOLVE_SOLID
 if ((norm_delta < nonlinear_tolerance) ){
    std::cout << "Nonlinear solver converged at step "<<std::endl; break;} 
#endif
        } // end nonlinear loop


//#include "write_results.txt" 

 last_non_linear_soln.update();      
 newton_update.update(); 
 // Write out every nth timestep to file.
 const unsigned int write_interval = 1;
#include "write_variable_results_and_mesh.txt"
    } // end timestep loop.
    // All done.  
  return 0;
}
