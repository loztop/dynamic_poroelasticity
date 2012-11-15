#include "defines.h"
#include "assemble.h"
#include "anal_neo_cc.h"

#include <cmath>
using namespace std;
#define PI 3.14159265
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_bcs (EquationSystems& es)
{

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
TransientLinearImplicitSystem & newton_update =
    es.get_system<TransientLinearImplicitSystem> ("Newton-update");

TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

   const System & ref_sys = es.get_system("Reference-Configuration"); 
  
  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = last_non_linear_soln.variable_number ("u");
  const unsigned int v_var = last_non_linear_soln.variable_number ("v");
  const unsigned int w_var = last_non_linear_soln.variable_number ("w");
#if INCOMPRESSIBLE
  const unsigned int p_var = last_non_linear_soln.variable_number ("p");
#endif 
FEType fe_vel_type = last_non_linear_soln.variable_type(u_var);
FEType fe_vel_type_ref = ref_sys.variable_type(u_var);

  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_pres_type = last_non_linear_soln.variable_type(p_var);
  AutoPtr<FEBase> fe_vel  (FEBase::build(3, fe_vel_type));
  AutoPtr<FEBase> fe_pres (FEBase::build(3, fe_pres_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

std::vector< int > rows;
 
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");
    const Real time    = es.parameters.get<Real>("time");

 
 //Build face
 #if TRACTION_BC
AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));          
AutoPtr<QBase> qface(fe_vel_type.default_quadrature_rule(3-1));
fe_face->attach_quadrature_rule (qface.get());

AutoPtr<FEBase> fe_face_ref (FEBase::build(3, fe_vel_type_ref));          
AutoPtr<QBase> qface_ref(fe_vel_type_ref.default_quadrature_rule(3-1));
fe_face_ref->attach_quadrature_rule (qface_ref.get());


    // AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_vel_type));
    // QGauss qface (dim-1, fe_vel_type.default_quadrature_order());
    // const std::vector<Real>& JxW_face = fe_face->get_JxW();
    // const std::vector<std::vector<Real> >& psi_face = fe_face->get_phi();	
#endif
 
  const DofMap & dof_map = last_non_linear_soln.get_dof_map();

test(81);
  // K will be the jacobian
  // F will be the Residual
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), 
    Kvu(Ke), Kvv(Ke), Kvw(Ke), 
    Kwu(Ke), Kwv(Ke), Kww(Ke); 
    
#if INCOMPRESSIBLE
  DenseSubMatrix<Number>  Kup(Ke),Kvp(Ke),Kwp(Ke), Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);
 #endif;
    
  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);
#if INCOMPRESSIBLE
  DenseSubVector<Number>    Fp(Fe);
#endif
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  
#if INCOMPRESSIBLE
  std::vector<unsigned int> dof_indices_p;
#endif
 
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

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
      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size(); 
      const unsigned int n_w_dofs = dof_indices_w.size(); 
#if INCOMPRESSIBLE
      const unsigned int n_p_dofs = dof_indices_p.size();
#endif

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

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






//#include "breathing_cube_solid_bcs.txt"
//#include "airway_blow_solid_bcs.txt"

//#include "confined_compression_solid_bcs.txt"

//#include "unconfined_compression_solid_bcs.txt"

#include "pipe_solid_bcs.txt"


      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

#if TRACTION_BC 
      newton_update.matrix->add_matrix (Ke, dof_indices);
      newton_update.rhs->add_vector    (Fe, dof_indices);
#endif

} // end of element loop

#if DIRICHLET_CLASSIC
        newton_update.matrix->close();
	newton_update.matrix->zero_rows(rows, 1.0);
#endif

     newton_update.rhs->close();
     newton_update.matrix->close();
     std::cout<<"Solid rhs->l2_norm () "<<newton_update.rhs->l2_norm ()<<std::endl;
//newton_update.rhs->print(std::cout);
  return;
}



