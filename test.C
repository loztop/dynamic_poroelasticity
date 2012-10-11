#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"
#include "poro_elastic_cc.h"

using namespace std;

void test(int a) {
<<<<<<< HEAD
//std::cout << "ex "<< a << std::endl;
=======
//zstd::cout << "ex "<< a << std::endl;
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
}


void test2(int a) {
std::cout << "test2 "<< a << std::endl;
}

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
} 


#if PORO
void compute_fluid_input(EquationSystems& es, unsigned int qp, const Elem* elem, Real m_old, Real& p_fluid, Real& J){


J=12;

p_fluid=23;

   TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");
      System& aux_system = es.get_system<System>("Reference-Configuration");

 const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  FEType fe_vel_type = last_non_linear_soln.variable_type(0);
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  std::vector<unsigned int> dof_indices_u;
const DofMap & dof_map = last_non_linear_soln .get_dof_map();
dof_map.dof_indices (elem, dof_indices_u, 0);
const unsigned int n_u_dofs = dof_indices_u.size(); 


const unsigned int p_var = last_non_linear_soln .variable_number ("p");
FEType fe_pres_type = last_non_linear_soln .variable_type(p_var);
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  std::vector<unsigned int> dof_indices_p;
dof_map.dof_indices (elem, dof_indices_p, p_var);
const unsigned int n_p_dofs = dof_indices_p.size(); 


  PoroelasticConfig material(dphi,psi);

            	            	std::cout << "lol " <<  std::endl;
            	            	std::cout << "elem " <<elem<<  std::endl;

//fe_vel->reinit  (elem);
fe_pres->reinit (elem);
//material


//Compute required variables from the solid
    VectorValue<Gradient> grad_u_mat;  
          std::vector<unsigned int> undefo_index;
      
grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> u_undefo;
      aux_system.get_dof_map().dof_indices(elem, undefo_index,d);
      aux_system.current_local_solution->get(undefo_index, u_undefo);
      for (unsigned int l = 0; l != n_u_dofs; l++){
    //    grad_u_mat(d).add_scaled(dphi[l][qp], u_undefo[l]); 
      }
    }


Real   p_solid = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
            {
            	std::cout << "psi " << psi.size() << std::endl;

            	std::cout << "psi[l][qp] " << psi[l][qp] << std::endl;
            	std::cout << "last_non_linear_soln.current_local_solution->el(dof_indices_p[l]) " << last_non_linear_soln.current_local_solution->el(dof_indices_p[l]) << std::endl;

             //p_solid += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);

            }

material.init_for_qp(grad_u_mat, p_solid, 0,0,0);

std::cout << "material.J " << material.J << std::endl;

std::cout << "material.p_solid " << material.p_solid << std::endl;




}

#endif

void print_resluts(EquationSystems& equation_systems, Real time, Real l){
 const unsigned int write_interval = 1;
      unsigned int t_step    = equation_systems.parameters.get<unsigned int>("step");
ExodusII_IO exo= ExodusII_IO(equation_systems.get_mesh());
 //  Real time    = equation_systems.parameters.get<unsigned int>("time");

#ifdef LIBMESH_HAVE_EXODUS_API
      if ((t_step+1)%write_interval == 0)
        { 
  #if STATIC  
        OStringStream file_name;
        file_name << "/ecslab/lorenzb/Dropbox/libresults/cube_with_stokes";
<<<<<<< HEAD

=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
        OSSRealzeroright(file_name,3,0, l + 1);
        file_name << ".e";    
  ExodusII_IO(mesh).write_equation_systems (file_name.str(),equation_systems);

  #if WRITE_MESH
  GmshIO(mesh).write(mesh_out_file_name);
  #endif
  #endif  

  #if DYNAMIC
 
  std::stringstream file_name;
  #if ! DAHOAM 
<<<<<<< HEAD
//  file_name << "/ecslab/lorenzb/Dropbox/libresults/";
                  file_name << "lolz";

=======
  file_name << "/ecslab/lorenzb/Dropbox/libresults/";
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
  #endif
  #if DAHOAM 
  file_name << "/home/lorenz/Dropbox/libresults/";
  #endif
<<<<<<< HEAD
  //file_name << "poro_fill_non_linear_t_step_hex"<< l<<"_";
  file_name << std::setw(2) << std::setfill('0') << l;
  file_name << ".e-s.";
  file_name << std::setw(3) << std::setfill('0') << l+1;
  std::cout<<" Printed "<< file_name.str() <<std::endl;

=======
  file_name << "poro_fill_non_linear_t_step_hex"<< l<<"_";
  file_name << std::setw(2) << std::setfill('0') << l;
  file_name << ".e-s.";
  file_name << std::setw(3) << std::setfill('0') << l+1;
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
  exo.write_timestep(file_name.str(), equation_systems,l+1,time);
std::cout<<" Printed "<< file_name.str() <<std::endl;

#if WRITE_TEC
  std::stringstream file_name_tec;
  file_name_tec << "cube_darcy"<< l<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
#endif
#endif
     }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
}
