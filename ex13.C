//export PETSC_ARCH=linux-gnu   (debug)
//export METHOD=dbg (debug libmesh)

#include "defines.h"
#include "assemble.h"
#include <iostream>
#include <fstream>

// The main program.
int main (int argc, char** argv)
{




  // Initialize libMesh.
LibMeshInit init (argc, argv);

//setup mesh and equations system  
#if CUSTOM_MESH
  std::string mesh_file_name ("/ecslab/lorenzb/gmsh-2.5.0-Linux/lozo_annulus_coarse.msh");
//    std::string mesh_file_name ("pipethick.msh");
#endif
#if WRITE_MESH
std::string mesh_out_file_name ("parallel.msh");
#endif

Mesh mesh(3);
MeshData mesh_data(mesh);

#if CUSTOM_MESH
 GmshIO(mesh).read(mesh_file_name);
#endif

#if CUBE
 /*
 MeshTools::Generation::build_cube (mesh,
                                       10, 5, 5,
                                       0.0, 4.0,	
                                       -0.5, 0.5,
				                               0., 1,
                                       HEX27);
*/

 MeshTools::Generation::build_cube (mesh,
                                       3, 2,2,
                                       0.0, 1.5,  
                                       0.0, 1.5,
                                       0.0, 1.5,
                                       HEX27);

#endif
mesh.all_second_order();
mesh.prepare_for_use();
mesh.print_info();
EquationSystems equation_systems (mesh);  
setup_es(equation_systems);
Real time     = 0;
Real end_time     = 25;


const unsigned int n_nonlinear_steps = 20;
const Real nonlinear_tolerance       = 1.e-2;
const Real initial_linear_solver_tol = 1.e-8;
unsigned int n_timesteps = 25;
Real dt = end_time/n_timesteps;



ExodusII_IO exo= ExodusII_IO(equation_systems.get_mesh());
#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif
Real dt_n_minus_1 = dt;

const char* plot_out_file_name ="plot.txt";
std::string result_file_name ("Chap_transfer_no_slip");

if (argc >2){
read_options(n_timesteps,result_file_name,plot_out_file_name,argc, argv);
}

dt = end_time/n_timesteps;

#if CHAP
ofstream outFile;         // Step #2 - Declare outstream file.
outFile.open(plot_out_file_name, std::ios_base::app);        // Step #3 - Open outFile.
//outFile <<"\n" << result_file_name  << " " << n_timesteps << "\n";

#endif


#if LOG_PERFORMANCE
  PerfLog perf_log("Solid Solver");
#endif   

equation_systems.parameters.set<Real> ("dt")   = dt;
equation_systems.init ();
//equation_systems.print_info();

// Get a reference to the systems to use later.
TransientLinearImplicitSystem&  newton_update =   equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
TransientLinearImplicitSystem&  last_non_linear_soln = equation_systems.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
TransientLinearImplicitSystem&  reference =   equation_systems.get_system<TransientLinearImplicitSystem>("Reference-Configuration");
#if FLUID_VEL
TransientLinearImplicitSystem&  fluid_system_vel = equation_systems.get_system<TransientLinearImplicitSystem>("fluid-system-vel");
#endif
#if ASSEMBLE_RESULTS
TransientLinearImplicitSystem&  results_system =    equation_systems.get_system<TransientLinearImplicitSystem>("pressure-system-lol");
#endif
#if ANALNEO
TransientLinearImplicitSystem&  anal =    equation_systems.get_system<TransientLinearImplicitSystem>("anal-system");
#endif

#if PETSC_MUMPS
  PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
  PC pc = petsc_linear_solver->pc();
  int ierr = PCSetType(pc, PCLU);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,"mumps");
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

#if UN_MINUS_ONE
TransientLinearImplicitSystem&  unm1 =
        equation_systems.get_system<TransientLinearImplicitSystem>("unm1-system");
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
  if ( abs((*node)(0)-0.5) + abs((*node)(1)-0.5) + abs((*node)(2)-0.5) < 0.001){
      node_1 = *nd; 
    std::cout<<"node1 "<< (*node_1)(0)<<std::endl;
  }
  if ( abs((*node)(0)-1.5) + abs((*node)(1)-1.5) + abs((*node)(2)-1.5) < 0.001){
      node_2 = *nd; 
    std::cout<<"node2 "<< (*node_2)(0)<<std::endl;
  }

    for (unsigned int d = 0; d < 3; ++d) {

      unsigned int dest_dof = node->dof_number(reference.number(), d, 0);
      Real value = (*node)(d);
      reference.current_local_solution->set(dest_dof, value);
      reference.solution->set(dest_dof, value);
    }


#if ANALNEO
   Real a= 0.1;
Real b= 0.1;
Real c1= 0.1;
      unsigned int dest_dof_X = node->dof_number(anal.number(), 0, 0);
      Real value_X = (*node)(0);
      Real x_val=value_X + a*(value_X*value_X)/2.0;
      anal.current_local_solution->set(dest_dof_X, x_val);
      anal.solution->set(dest_dof_X, x_val);
    
     unsigned int dest_dof_Y = node->dof_number(anal.number(), 1, 0);
      Real value_Y = (*node)(1);
      Real y_val=value_Y + b*(value_Y*value_Y)/2.0;
      anal.current_local_solution->set(dest_dof_Y, y_val);
      anal.solution->set(dest_dof_Y, y_val);

      unsigned int dest_dof_Z = node->dof_number(anal.number(), 2, 0);
      Real value_Z = (*node)(2);
      Real z_val=value_Z*( 1.0/( 1+ a*value_X )  )*( 1.0/( 1+ b*value_Y )  );
      anal.current_local_solution->set(dest_dof_Z, z_val);
      anal.solution->set(dest_dof_Z, z_val);
#endif
/*
      unsigned int dest_dof = node->dof_number(reference.number(), 3, 0);
      Real value = 0;

      if (dest_dof <12345678){ 
      reference.current_local_solution->set(dest_dof, value);
      reference.solution->set(dest_dof, value);
      }
*/
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
      time += dt;
      equation_systems.parameters.set<Real> ("time") = time;
      double progress = (t_step+0.00000000001) / (n_timesteps+0.00000000001);
      equation_systems.parameters.set<Real>("progress") = progress;
      equation_systems.parameters.set<unsigned int>("step") = t_step; 


#if ANALNEO

 for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
    const Node *node = *nd; 

   Real a= 0.1*progress;
Real b= 0.1*progress;
Real c1= 0.1;

      unsigned int dest_dof_X = node->dof_number(anal.number(), 0, 0);
      Real value_X = reference.current_local_solution->el(dest_dof_X);
      Real x_val=value_X + a*(value_X*value_X)/2.0;
      anal.current_local_solution->set(dest_dof_X, x_val);
      anal.solution->set(dest_dof_X, x_val);
    
     unsigned int dest_dof_Y = node->dof_number(anal.number(), 1, 0);
      Real value_Y = reference.current_local_solution->el(dest_dof_Y);
      Real y_val=value_Y + b*(value_Y*value_Y)/2.0;
      anal.current_local_solution->set(dest_dof_Y, y_val);
      anal.solution->set(dest_dof_Y, y_val);

      unsigned int dest_dof_Z = node->dof_number(anal.number(), 2, 0);
      Real value_Z = reference.current_local_solution->el(dest_dof_Z);
      Real z_val=value_Z*( 1.0/( 1+ a*value_X )  )*( 1.0/( 1+ b*value_Y )  );
      anal.current_local_solution->set(dest_dof_Z, z_val);
      anal.solution->set(dest_dof_Z, z_val);

/*

 unsigned int dest_dof_X = node->dof_number(anal.number(), 0, 0);
      Real value_X = reference.current_local_solution->el(dest_dof_X);
      Real x_val=value_X ;
      anal.current_local_solution->set(dest_dof_X, x_val);
      anal.solution->set(dest_dof_X, x_val);
    
     unsigned int dest_dof_Y = node->dof_number(anal.number(), 1, 0);
      Real value_Y = reference.current_local_solution->el(dest_dof_Y);
      Real y_val=value_Y;
      anal.current_local_solution->set(dest_dof_Y, y_val);
      anal.solution->set(dest_dof_Y, y_val);

      unsigned int dest_dof_Z = node->dof_number(anal.number(), 2, 0);
      Real value_Z = reference.current_local_solution->el(dest_dof_Z);
      Real z_val=0.5*value_Z*value_X;
      anal.current_local_solution->set(dest_dof_Z, z_val);
      anal.solution->set(dest_dof_Z, z_val);
*/
}

#endif









#if ADAPTIVE_TIME
if (progress>0.09){ equation_systems.parameters.set<Real> ("dt") = 0.05; }
if (progress>0.2){  equation_systems.parameters.set<Real> ("dt") = 0.1;  }
if (progress>0.5){  equation_systems.parameters.set<Real> ("dt") = 0.1;  }
#endif
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


#endif 
	//update the final solution
	// xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
	//Apply a full Newton-step
Real K=1; //Newton step size
	last_non_linear_soln.solution->add(-1,*newton_update.solution);
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
/*
//Find some interesting dofs to inspect
Node node = mesh.node(10);
//std::cout<<" node a " <<  node  <<std::endl;

unsigned int u_dof_a = node.dof_number(1, 0, 0);
unsigned int v_dof_a = node.dof_number(1, 1, 0);
unsigned int w_dof_a = node.dof_number(1, 2, 0);
unsigned int p_dof_a = node.dof_number(1, 3, 0);
*/

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

 if ((norm_delta < 1230.1) ){

  //std::cout<<" Solving fluid " <<  std::endl;

//----------------- Solve the Fluid -------------------------------///
#if FLUID
//Save previous solution so we can compare against teh new solution
change_in_fluid->zero();
change_in_fluid->add(*fluid_system.solution);

fluid.solution->zero();  
fluid.update();

#if !PETSC_MUMPS
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 2000;
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

/*
#if CHAP && ASSEMBLE_RESULTS
Point p1(0,0,0);
Point p2(0.5,0.5,0.5); //These co ords are not in the reference pos
Point p3(1.5,1.5,1.5);

//outFile << " time: "<< time << "\n";

//outFile << results_system.point_value(0,p1) << " ";
//outFile<<results_system.point_value(0,p2)<< " ";
//outFile<<results_system.point_value(0,p3)<< " ";

outFile << results_system.solution->el(node_0->dof_number(reference.number(), 1, 0)) << " ";
outFile << results_system.solution->el(node_1->dof_number(reference.number(), 1, 0)) << " ";
outFile << results_system.solution->el(node_2->dof_number(reference.number(), 1, 0)) << " ";

outFile << fluid_system_vel.solution->el(node_0->dof_number(reference.number(), 3, 0)) << " ";
outFile << fluid_system_vel.solution->el(node_1->dof_number(reference.number(), 3, 0)) << " ";
outFile << fluid_system_vel.solution->el(node_2->dof_number(reference.number(), 3, 0)) << " ";


outFile << results_system.point_value(1,p1) << " ";
outFile<<results_system.point_value(1,p2)<< " ";
outFile<<results_system.point_value(1,p3)<< " ";

outFile << fluid_system_vel.point_value(3,p1) << " ";
outFile<<fluid_system_vel.point_value(3,p2)<< " ";
outFile<<fluid_system_vel.point_value(3,p3)<< " ";

outFile<< endl;   
#endif
*/

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
//  file_name << "/ecslab/lorenzb/Dropbox/";
      file_name << "SIM_";

  file_name << result_file_name;
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
