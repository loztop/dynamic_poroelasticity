#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"
#include "poro_elastic_cc.h"

#if ASSEMBLE_PRESSURE
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_pressure (EquationSystems& es,
                      const std::string& system_name)
{

    libmesh_assert (system_name == "pressure-system");
    
    const MeshBase& mesh = es.get_mesh();
    
    const unsigned int dim = mesh.mesh_dimension();
    
    LinearImplicitSystem & system =
      es.get_system<LinearImplicitSystem> ("pressure-system");
  
  #if FLUID
 TransientLinearImplicitSystem & fluid_system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system");
#endif






    const unsigned int u_var = system.variable_number ("fluid_pressure");
    const unsigned int v_var = system.variable_number ("jacobian");
    const unsigned int w_var = system.variable_number ("p_aux1");
    const unsigned int p_var = system.variable_number ("p_aux2");
    
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
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe),
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;
        const Real dt    = es.parameters.get<Real>("dt");
   unsigned int step    = es.parameters.get<unsigned int>("step");
    const Real non_lin_step    = es.parameters.get<Real>("non_lin_step");


  
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

   TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");
      System& aux_system = es.get_system<System>("Reference-Configuration");
Real   p_solid = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
            {
         p_solid += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
            }

Real   w_disp = 0.;
for (unsigned int l=0; l<n_w_dofs; l++)
       {
         w_disp += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_w[l]);
      }

Real   u_disp = 0.;
for (unsigned int l=0; l<n_w_dofs; l++)
       {
         u_disp += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_u[l]);
      }

    VectorValue<Gradient> grad_u_mat;  
          std::vector<unsigned int> undefo_index; 
  grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> u_undefo;
      aux_system.get_dof_map().dof_indices(elem, undefo_index,d);
      aux_system.current_local_solution->get(undefo_index, u_undefo);
      for (unsigned int l = 0; l != n_u_dofs; l++){
     grad_u_mat(d).add_scaled(dphi[l][qp], u_undefo[l]); 
      }
    }

Real J; 
Real m_old=0;
Real m_current=0;

#if FLUID
for (unsigned int l=0; l<n_u_dofs; l++)
{
 m_old += phi[l][qp]*fluid_system.old_local_solution->el(dof_indices_u[l]);
 m_current += phi[l][qp]*fluid_system.current_local_solution->el(dof_indices_u[l]);
}
#endif

#if PORO
PoroelasticConfig material(dphi,phi);
#endif

#if !PORO && INCOMPRESSIBLE
NonlinearNeoHookeCurrentConfig material(dphi,phi);
#endif

#if !PORO && COMPRESSIBLE
NonlinearNeoHookeCurrentConfig material(dphi);
material.init_for_qp(grad_u_mat, p_solid, 0);
#endif

#if INCOMPRESSIBLE
Real p_fluid=0;
material.init_for_qp(grad_u_mat, p_solid, 0, m_old,p_fluid);
#endif

#if PORO
material.calculate_fluid_pressure();
p_fluid=material.p_fluid;
#endif
if( (step ==70) && (non_lin_step <1)  ) {
  std::cout << "material.C " << material.C << std::endl;
  std::cout << "material.I_3 " << material.I_3 << std::endl;
  std::cout << "material.I_1 " << material.I_1 << std::endl;
  std::cout << "material.D " << material.D << std::endl;
  std::cout << "material.m " << material.m << std::endl;
  std::cout << "material.f " << material.f << std::endl;
  std::cout << "material.p_solid" << material.p_solid << std::endl;
 // std::cout << "material.p_fluid " << material.p_fluid << std::endl;
}

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);

       for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);

       for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++){
              Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
             //  Kww(i,j) += -JxW[qp]*(phi[i][qp]*dphi[j][qp](0));

            }
       for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++)
              Kpp(i,j) += JxW[qp]*(psi[i][qp]*psi[j][qp]);
       #if PORO || FLUID
          for (unsigned int j=0; j<n_u_dofs; j++){
        Fu(j) += JxW[qp]*phi[j][qp]*p_fluid;
           }
      #endif

          for (unsigned int j=0; j<n_v_dofs; j++){
             Fv(j) += JxW[qp]*phi[j][qp]*material.J;
           // Fv(j) += -JxW[qp](dphi[j][qp](0) + dphi[j][qp](1) +dphi[j][qp](2)   )*1000+78;
          }


           for (unsigned int j=0; j<n_v_dofs; j++){
             //Fw(j) += -JxW[qp](dphi[j][qp](0) + dphi[j][qp](1) +dphi[j][qp](2))*1000;
            // Fw(j) += -JxW[qp]*(0*dphi[j][qp](0) + 1*dphi[j][qp](1) +0*dphi[j][qp](2))*1;

           // Fw(j) += JxW[qp]*(phi[j][qp])*1;

#if FLUID
           for (unsigned int i=0; i<n_v_dofs; i++){
            Fw(j) += JxW[qp]*(phi[j][qp])*dphi[i][qp](1)*p_fluid;
           }
#endif
            }

         for (unsigned int j=0; j<n_p_dofs; j++){
             Fp(j) += JxW[qp]*psi[j][qp]*(p_solid);
            //  Fp(j) += JxW[qp]*psi[j][qp]*p_solid;

          }
        } // end of the quadrature point qp-loop

        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      } // end of element loop
    
    return;    
  }

#endif



