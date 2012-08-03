

#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
#include "solid_system.h"
#include "poro_elastic_cc.h"


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_fluid (EquationSystems& es,
                      const std::string& system_name)
{


    libmesh_assert (system_name == "fluid-system");
    
    const MeshBase& mesh = es.get_mesh();
    
    const unsigned int dim = mesh.mesh_dimension();
    

  
 TransientLinearImplicitSystem &  system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system");

 TransientLinearImplicitSystem &  fluid_system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system");

    TransientLinearImplicitSystem & pressure_system =
    es.get_system<TransientLinearImplicitSystem> ("pressure-system");

const System & ref_sys = es.get_system("Reference-Configuration"); 

    const unsigned int u_var = system.variable_number ("fluid_m");
    const unsigned int v_var = system.variable_number ("fluid_aux1");
    const unsigned int w_var = system.variable_number ("fluid_aux2");
    const unsigned int p_var = system.variable_number ("fluid_aux3");
    
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
    

AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));
                
    AutoPtr<QBase> qface(fe_vel_type.default_quadrature_rule(3-1));
  
    fe_face->attach_quadrature_rule (qface.get());

    const DofMap & dof_map = system.get_dof_map();
  
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
  
    DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe),
    Fp(Fe);
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;
#if DYNAMIC
  VectorValue<Gradient> grad_u_mat_old;
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");
       unsigned int step    = es.parameters.get<unsigned int>("step");
           const Real non_lin_step    = es.parameters.get<Real>("non_lin_step");
std::cout<<" non_lin_step " << non_lin_step<< std::endl;

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
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size(); 
      const unsigned int n_w_dofs = dof_indices_w.size(); 
      const unsigned int n_p_dofs = dof_indices_p.size();
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
      
      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      
      Kwu.reposition (w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
      Kwp.reposition (w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);
      
      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      
      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Assemble the u-velocity row
          // uu coupling
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++){
             Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);        
           }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++){
              Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
            }
          for (unsigned int i=0; i<n_w_dofs; i++)
            for (unsigned int j=0; j<n_w_dofs; j++){
              Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
            }
              for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++){
              Kpp(i,j) += JxW[qp]*(psi[i][qp]*psi[j][qp]);
            }

std::vector<unsigned int> undefo_index; 
VectorValue<Real> grad_p_mat;  
std::vector<Number> press_defo;
pressure_system.get_dof_map().dof_indices(elem, undefo_index,0);
pressure_system.current_local_solution->get(undefo_index, press_defo);
for (unsigned int l = 0; l != n_u_dofs; l++){
     grad_p_mat.add_scaled(dphi[l][qp], press_defo[l] ); 
}

   // std::cout << "grad_p_mat inn " <<grad_p_mat<< std::endl;


VectorValue<Real> grad_j_mat;  
pressure_system.get_dof_map().dof_indices(elem, undefo_index,1);
      pressure_system.current_local_solution->get(undefo_index, press_defo);
for (unsigned int l = 0; l != n_v_dofs; l++){
     grad_j_mat.add_scaled(dphi[l][qp], press_defo[l] ); 
}



Real m_old=0;
for (unsigned int l=0; l<n_u_dofs; l++)
{
  m_old += phi[l][qp]*fluid_system.old_local_solution->el(dof_indices_u[l]);
}


Real u_ref=0;
for (unsigned int l=0; l<n_u_dofs; l++)
{
  u_ref += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
}

//std::cout<<"u_ref "<< u_ref<<std::endl;

Real q=0;
Point point = elem->centroid();
//if ( (point(0)<1.4 )   && (point(0)>0.5 )  &&  (point(1)<1.4 )  && (point(1)>0.1 ) &&   (point(2)<1.4 )   && (point(2)>0.1 )  ){
//if ( (point(0)>0.8)    ){
//q= point(0)*0.2;
q= u_ref*0.1;

//}

Real J=0;
for (unsigned int l=0; l<n_v_dofs; l++)
{
  J += phi[l][qp]*pressure_system.current_local_solution->el(dof_indices_v[l]);
}

Real   p_fluid = 0.;
for (unsigned int l=0; l<n_v_dofs; l++)
            {
         p_fluid += phi[l][qp]*pressure_system.current_local_solution->el(dof_indices_u[l]);
            }

            for (unsigned int j=0; j<n_u_dofs; j++){
            
              //Pressure contribution
            //  Fu(j) += -J*dt*JxW[qp]*dphi[j][qp]*(grad_p_mat);
//New DK version

              Fu(j) += -dt*JxW[qp]*(grad_p_mat)*(phi[j][qp]*(grad_j_mat) +J*dphi[j][qp]);

              //Source contribution
              Fu(j) += J*dt*JxW[qp]*phi[j][qp]*q;
              //Previous time step contribution
              Fu(j) += JxW[qp]*phi[j][qp]*m_old;

            }

             for (unsigned int j=0; j<n_v_dofs; j++){
             // Fv(j) += JxW[qp]*grad_p_mat(0)*phi[j][qp];
              Fv(j) += -JxW[qp]*dphi[j][qp](0)*p_fluid;
            }
             for (unsigned int j=0; j<n_w_dofs; j++){
              Fw(j) += -JxW[qp]*dphi[j][qp](0)*(grad_p_mat(0));
          }

            for (unsigned int j=0; j<n_p_dofs; j++){
            //  Fp(j) += -J*dt*JxW[qp]*dpsi[j][qp]*(grad_p_mat);
          }


if( (step ==157) && (non_lin_step <2)  ) {
  std::cout << "J " << J<< std::endl;
    std::cout << "grad_p_mat " <<grad_p_mat<< std::endl;
}
        } // end of the quadrature point qp-loop





for (unsigned int s=0; s<elem->n_sides(); s++){
   if (elem->neighbor(s) == NULL)
   {    
   AutoPtr<Elem> side (elem->build_side(s));
          const Real penalty = 1;

      const std::vector<std::vector<Real> >&  phi_face =   fe_face->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
        const std::vector<Real>& JxW_face = fe_face->get_JxW();
        const std::vector<Point>& qface_point = fe_face->get_xyz();
        const std::vector<Point>& face_normals =fe_face->get_normals();
        fe_face->reinit(elem,s);  

  for (unsigned int qp=0; qp<qface->n_points(); qp++)
                  {

Real   p_fluid = 0.;
for (unsigned int l=0; l<phi_face.size(); l++){
  p_fluid += phi_face[l][qp]*pressure_system.current_local_solution->el(dof_indices_u[l]);
}

Real J=0;
for (unsigned int l=0; l<n_v_dofs; l++)
{
  J += phi[l][qp]*pressure_system.current_local_solution->el(dof_indices_v[l]);
}

VectorValue<Real> grad_p_mat;  
for (unsigned int l = 0; l != n_u_dofs; l++){
     grad_p_mat.add_scaled(dphi_face[l][qp], pressure_system.current_local_solution->el(dof_indices_u[l])); 
}

    //std::cout << "grad_p_mat boud " <<grad_p_mat<< std::endl;

                    const Number  value = 0;
                                                         
                    for (unsigned int i=0; i<phi_face.size(); i++){

           Fu(i) += dt*JxW_face[qp]*J*phi_face[i][qp]*face_normals[qp]*grad_p_mat;

           Fv(i) += JxW_face[qp]*face_normals[qp](0)*phi_face[i][qp]*p_fluid;

           Fw(i) += JxW_face[qp]*face_normals[qp](0)*phi_face[i][qp]*grad_p_mat(0);



        }
    
    } //end qp
}}


        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      } // end of element loop
    
    return;
  }


