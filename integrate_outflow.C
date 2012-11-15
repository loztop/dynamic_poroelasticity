#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"
#include "poro_elastic_cc.h"
#include "mooney_cc.h"

// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void integrate_outflow (EquationSystems& es)
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


AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));             
AutoPtr<QBase> qface(fe_vel_type.default_quadrature_rule(3-1));  
fe_face->attach_quadrature_rule (qface.get());

    // AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_vel_type));
    // QGauss qface (dim-1, fe_vel_type.default_quadrature_order());
    // const std::vector<Real>& JxW_face = fe_face->get_JxW();
    // const std::vector<std::vector<Real> >& psi_face = fe_face->get_phi();	

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 


std::vector< int > rows;
std::vector< int > pressure_srows;

Number outflow = 0;

    for ( ; el != end_el; ++el)
    {       
      const Elem* elem = *el;
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
test(2);


  const std::vector<std::vector<Real> >&  phi_face =  fe_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Point>& face_normals = fe_face->get_normals();
  fe_face->reinit(elem,s);  

for (unsigned int ns=0; ns<side->n_nodes(); ns++)
    {
       for (unsigned int n=0; n<elem->n_nodes(); n++){
    Node *node = elem->get_node(n);
       Point p;
  for (unsigned int d = 0; d < 3; ++d) {
        unsigned int source_dof = node->dof_number(1, d, 0);
        Real value = aux_system.current_local_solution->el(source_dof);
        p(d)=value;
    }
    if ((elem->node(n) == side->node(ns)) && (p(0)<0.01)   )

    for (unsigned int qp=0; qp<qface->n_points(); qp++)
    {       
  if ((elem->node(n) == side->node(ns)) && (p(0)<0.01)   )
{




  for (unsigned int l=0; l<n_u_dofs; l++)
{ 
   outflow += system.current_local_solution->el(dof_indices_u[l])*phi_face[l][qp];
        } 

  } //end if 

 

  }    
}



for (unsigned int ns=0; ns<side->n_nodes(); ns++)
    {
     for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
	     Point p;
	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
      	Real value = aux_system.current_local_solution->el(source_dof);
      	p(d)=value;
   	}
    
 
} // end nodes on side loop

}//for loop
  
}

 }//end elem->neighbor(s) == NULL 
 }// end s=0; s<elem->n_sides(); s++

 
} // end of element loop

std::cout<< " outflow "<< outflow <<std::endl;

const char* plot_out_file_name ="outflow.txt";
 ofstream outFile;    
outFile.open(plot_out_file_name, std::ios_base::app);

outFile << outflow << " ";
outFile.close();
  return;
}




