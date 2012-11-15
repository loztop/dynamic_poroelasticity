#include "defines.h"
#include "assemble.h"
#include <cmath>


using namespace std;
#define PI 3.14159265



// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_fluid_bcs (EquationSystems& es)
{
    
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
    
 TransientLinearImplicitSystem &  system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system-vel");

   TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

    System& aux_system = es.get_system<System>("Reference-Configuration");


    const unsigned int u_var = system.variable_number ("fluid_U_vel");
    const unsigned int v_var = system.variable_number ("fluid_V_vel");
    const unsigned int w_var = system.variable_number ("fluid_W_vel");
    const unsigned int p_var = system.variable_number ("fluid_P");

#if FLUID_P_CONST
    const unsigned int m_var = system.variable_number ("fluid_M");
#endif
    
    FEType fe_vel_type = system.variable_type(u_var);
    
    FEType fe_pres_type = system.variable_type(p_var);
  
    AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
      
    AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  
    fe_vel->attach_quadrature_rule (&qrule);
    fe_pres->attach_quadrature_rule (&qrule);
    
    const std::vector<Real>& JxW = fe_vel->get_JxW();
    
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  
      const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();


    const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
    
    const DofMap & dof_map = system.get_dof_map();
  
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
  
    DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), Kum(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), Kvm(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), Kwm(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke), Kpm(Ke),
    Kmu(Ke), Kmv(Ke), Kmw(Ke), Kmp(Ke), Kmm(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe),
    Fp(Fe);
#if FLUID_P_CONST
   DenseSubVector<Number>
   Fm(Fe);
#endif

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;

#if FLUID_P_CONST
  std::vector<unsigned int> dof_indices_m;
#endif

#if DYNAMIC
VectorValue<Gradient> grad_u_mat_old;
const Real dt    = es.parameters.get<Real>("dt");
const Real progress    = es.parameters.get<Real>("progress");
unsigned int step    = es.parameters.get<unsigned int>("step");
const Real non_lin_step    = es.parameters.get<Real>("non_lin_step");
    const Real time    = es.parameters.get<Real>("time");
#endif



 #if NEUMANN_PRESSURE
AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));             
AutoPtr<QBase> qface(fe_vel_type.default_quadrature_rule(3-1));  
fe_face->attach_quadrature_rule (qface.get());


AutoPtr<FEBase> fe_face_pres (FEBase::build(3, fe_pres_type));             
AutoPtr<QBase> qface_pres(fe_pres_type.default_quadrature_rule(3-1));  
fe_face_pres->attach_quadrature_rule (qface_pres.get());	
#endif

  test(1);

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    test(2);

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

test(3);

std::vector< int > rows;
std::vector< int > pressure_rows;

    for ( ; el != end_el; ++el)
    {       
      const Elem* elem = *el;
test(67);

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_w, w_var);
#if INCOMPRESSIBLE
      dof_map.dof_indices (elem, dof_indices_p, p_var);
#endif
#if FLUID_P_CONST
      dof_map.dof_indices (elem, dof_indices_m, m_var);
#endif
      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size(); 
      const unsigned int n_w_dofs = dof_indices_w.size(); 
#if INCOMPRESSIBLE
      const unsigned int n_p_dofs = dof_indices_p.size();
#endif
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);
test(69);

      // Similarly, the \p DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      #if INCOMPRESSIBLE
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
      #endif
      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
      #if INCOMPRESSIBLE
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      #endif
      
      Kwu.reposition (w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
      #if INCOMPRESSIBLE
      Kwp.reposition (w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);
      
      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
      #endif
      
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      #if INCOMPRESSIBLE
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      #endif
test(1);
//Now start actually applying the BCS        
for (unsigned int s=0; s<elem->n_sides(); s++){
   if (elem->neighbor(s) == NULL)
   {		
   AutoPtr<Elem> side (elem->build_side(s));     

//#include "orig_setup_fluid_bcs.txt"
//#include "breathing_fluid_bcs.txt"

//#include "pipe_fluid_bcs.txt"
//#include "confined_compression_fluid_bcs.txt"

//#include "unconfined_compression_fluid_bcs.txt"

//#include "airway_blow_fluid_bcs.txt"

 }//end elem->neighbor(s) == NULL 
 }// end s=0; s<elem->n_sides(); s++

test(100);
 dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

//#if NEUMANN_PRESSURE
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
//#endif

} // end of element loop


#if DIRICHLET_VELOCITY
  system.matrix->close();
  system.matrix->zero_rows(rows, 1.0);
#endif

#if DIRICHLET_PRESSURE
   system.matrix->close();
   system.matrix->zero_rows(pressure_rows, 1.0);
#endif

test(5);
system.rhs->close();
test(6);

system.matrix->close();
 //    std::cout<<"Fluid rhs->l2_norm () "<<system.rhs->l2_norm ()<<std::endl;
//system.rhs->print(std::cout);
//system.matrix->print(std::cout);


  return;
}



