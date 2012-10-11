#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"
#include "poro_elastic_cc.h"

<<<<<<< HEAD

// The matrix assembly function to be called at each time step to
// prepare for the linear solve.

#if FLUID_VEL
=======
#if FLUID_VEL
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
void assemble_fluid_vel (EquationSystems& es,
                      const std::string& system_name)
{


    libmesh_assert (system_name == "fluid-system-vel");
    
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
<<<<<<< HEAD

#if FLUID_P_CONST
    const unsigned int m_var = system.variable_number ("fluid_M");
#endif
=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
    
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
<<<<<<< HEAD
      const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();

=======
    
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
    const DofMap & dof_map = system.get_dof_map();
  
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
  
    DenseSubMatrix<Number>
<<<<<<< HEAD
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), 
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), 
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), 
=======
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe),
    Fp(Fe);
<<<<<<< HEAD
#if FLUID_P_CONST
   DenseSubVector<Number>  Fm(Fe);
       DenseSubMatrix<Number>  Kwm(Ke), Kvm(Ke), Kum(Ke),Kmu(Ke), Kpm(Ke), Kmv(Ke), Kmw(Ke), Kmp(Ke), Kmm(Ke);
#endif
=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;
<<<<<<< HEAD

#if FLUID_P_CONST
  std::vector<unsigned int> dof_indices_m;
#endif

=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
#if DYNAMIC
VectorValue<Gradient> grad_u_mat_old;
const Real dt    = es.parameters.get<Real>("dt");
const Real progress    = es.parameters.get<Real>("progress");
unsigned int step    = es.parameters.get<unsigned int>("step");
const Real non_lin_step    = es.parameters.get<Real>("non_lin_step");
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
<<<<<<< HEAD
#if FLUID_P_CONST
      dof_map.dof_indices (elem, dof_indices_m, m_var);
#endif
=======

>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size(); 
      const unsigned int n_w_dofs = dof_indices_w.size(); 
      const unsigned int n_p_dofs = dof_indices_p.size();

<<<<<<< HEAD
#if FLUID_P_CONST
      const unsigned int n_m_dofs = dof_indices_m.size();
#endif
      
=======
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76

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

<<<<<<< HEAD




=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);

<<<<<<< HEAD
#if FLUID_P_CONST

      Kum.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs , n_u_dofs, n_m_dofs);
      Kvm.reposition (v_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs, n_v_dofs, n_m_dofs);
      Kwm.reposition (w_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs, n_w_dofs, n_m_dofs);
      Kpm.reposition (p_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs, n_p_dofs, n_m_dofs);

      Kmu.reposition (p_var*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_m_dofs, n_u_dofs);
      Kmv.reposition (p_var*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_m_dofs, n_v_dofs);
      Kmw.reposition (p_var*n_u_dofs + n_p_dofs, w_var*n_u_dofs, n_m_dofs, n_w_dofs);
      Kmp.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_m_dofs, n_p_dofs);
      Kmm.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs, n_m_dofs, n_m_dofs);

      Fm.reposition (p_var*n_u_dofs + n_p_dofs, n_m_dofs);
#endif


fe_vel->reinit  (elem);
fe_pres->reinit (elem);

=======
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

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


#if PORO
PoroelasticConfig material(dphi,phi);
#endif

#if INCOMPRESSIBLE
Real p_solid=0;
for (unsigned int l=0; l<n_p_dofs; l++)
{
  p_solid += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
}
#endif


#if INCOMPRESSIBLE
Real p_fluid=0;
Real m_old=0;
material.init_for_qp(grad_u_mat, p_solid, 0, m_old,p_fluid);
Real J=material.J;
Real fchap=material.fchap;
Real fchapd=material.fchapd;
Real M=material.M;
Real Kperm=material.Kperm;
#endif

  Number   div_vs = 0.;
  for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> u_undefo_current;
      std::vector<Number> u_undefo_old;
      last_non_linear_soln.get_dof_map().dof_indices(elem, undefo_index,d);
      last_non_linear_soln.current_local_solution->get(undefo_index, u_undefo_current);
      last_non_linear_soln.old_local_solution->get(undefo_index, u_undefo_old);
      for (unsigned int l = 0; l != n_u_dofs; l++)
        div_vs+=(dphi[l][qp](d)*u_undefo_current[l]-dphi[l][qp](d)*u_undefo_old[l])/dt; 
    }

<<<<<<< HEAD
    Number   vel_u_s = 0.;
    Number   vel_v_s = 0.;
    Number   vel_w_s = 0.;

for (unsigned int l=0; l<n_u_dofs; l++)
{
              vel_u_s += (phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_u[l])- phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_u[l]) )/dt;
              vel_v_s += (phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_v[l])- phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_v[l]) )/dt;
              vel_w_s += (phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_w[l])- phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_w[l]) )/dt;

}

//std::cout<< "div_vs "<< div_vs<<std::endl;

 Number   lambda_dot = 0.;
      std::vector<Number> u_undefo_current;
      std::vector<Number> u_undefo_old;
      last_non_linear_soln.get_dof_map().dof_indices(elem, undefo_index,3);
      last_non_linear_soln.current_local_solution->get(undefo_index, u_undefo_current);
      last_non_linear_soln.old_local_solution->get(undefo_index, u_undefo_old);
      for (unsigned int l = 0; l != n_p_dofs; l++)
        lambda_dot+=(phi[l][qp]*u_undefo_current[l]-phi[l][qp]*u_undefo_old[l])/dt; 
    

Number   lambda = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
{
              lambda += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
}

Number   lambda_old = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
{
              lambda_old += psi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_p[l]);
}


#if FLUID_P_CONST
Number   m_olD = 0.;
for (unsigned int l=0; l<n_m_dofs; l++)
{
              m_olD += psi[l][qp]*system.old_local_solution->el(dof_indices_m[l]);
}



Number   m_current = 0.;
for (unsigned int l=0; l<n_m_dofs; l++)
{
              m_current += psi[l][qp]*system.current_local_solution->el(dof_indices_m[l]);
}
#endif

//std::cout<<" dt " <<  dt<<std::endl;
=======


>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76

//Calculate the pressure from the previous time step at given qp.
Number   p_old = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
{
              p_old += psi[l][qp]*system.old_local_solution->el(dof_indices_p[l]);
}
#if DECOUPLE
div_vs=0.1;
p_old=0.1;
#endif
<<<<<<< HEAD

Real factor=1;


=======
Real factor=1;
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
             // Assemble the u-velocity row
            //Mass Matrix needed for darcy flow
          for (unsigned int i=0; i<n_u_dofs; i++){

            for (unsigned int j=0; j<n_u_dofs; j++){
              //w.v term (u)

              Kuu(i,j) +=  (1.0/Kperm)*JxW[qp]*(phi[i][qp]*phi[j][qp]);
              Kvv(i,j) +=  (1.0/Kperm)*JxW[qp]*(phi[i][qp]*phi[j][qp]);
              Kww(i,j) +=  (1.0/Kperm)*JxW[qp]*(phi[i][qp]*phi[j][qp]);

              // Laplacian for stokes flow /Brinkman
<<<<<<< HEAD
              Real brink_mu=  0.01;
=======
              Real brink_mu=0.001;
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
              Kuu(i,j) += factor*brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              Kvv(i,j) += factor*brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              Kww(i,j) += factor*brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

            }


          for (unsigned int j=0; j<n_p_dofs; j++){

        // Weak form of grad p term
              Kup(i,j) += -factor*JxW[qp]*psi[j][qp]*dphi[i][qp](0);
              Kvp(i,j) += -factor*JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              Kwp(i,j) += -factor*JxW[qp]*psi[j][qp]*dphi[i][qp](2);
<<<<<<< HEAD
}
}


#if ! FLUID_P_CONST && ! SOLID_P_CONST  

          for (unsigned int i=0; i<n_u_dofs; i++){
          for (unsigned int j=0; j<n_p_dofs; j++){
           //the velocity divergence term 
              Kpu(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
              Kpv(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              Kpw(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);   
          }
	}
=======


           //the velocity divergence term 
              Kpu(j,i) += dt/dt*JxW[qp]*psi[j][qp]*dphi[i][qp](0);
              Kpv(j,i) += dt/dt*JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              Kpw(j,i) += dt/dt*JxW[qp]*psi[j][qp]*dphi[i][qp](2);
            }
          }
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76


for (unsigned int i=0; i<n_p_dofs; i++){

        //Backward euler contribution to rhs
<<<<<<< HEAD
    //    Fp(i) +=  (p_old/(J*M*fchap))*JxW[qp]*psi[i][qp];
        //Fp(i) += -dt*div_vs*JxW[qp]*psi[i][qp];
        Fp(i) += -1*div_vs*JxW[qp]*psi[i][qp];

            if(lambda>0 && lambda_dot > 0){
      Fp(i) +=  1*(p_old/(2*J*lambda*dt))*JxW[qp]*psi[i][qp];
}

  for (unsigned int j=0; j<n_p_dofs; j++){
  //  Kpp(i,j) += (1/(J*M*fchap))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
  //  Kpp(i,j) += -dt*(fchapd/(M*pow(fchap,2.0)))*div_vs*JxW[qp]*psi[i][qp]*(psi[j][qp]);

    if(lambda>0 && lambda_dot > 0){
   Kpp(i,j) += 1*(1.0/(2*J*lambda_dot))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
Kpp(i,j) += 1*(1.0/(2*J*lambda*dt))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
}
    }    
  }

#endif


#if SOLID_P_CONST

        for (unsigned int i=0; i<n_u_dofs; i++){
          for (unsigned int j=0; j<n_p_dofs; j++){
           //the velocity divergence term 
              Kpu(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
              Kpv(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              Kpw(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);   
          }
  }


for (unsigned int i=0; i<n_p_dofs; i++){
        //Fp(i) += -1*div_vs*JxW[qp]*psi[i][qp];
        Fp(i) += -1*vel_u_s*JxW[qp]*dpsi[i][qp](0);
        Fp(i) += -1*vel_v_s*JxW[qp]*dpsi[i][qp](1);   
        Fp(i) += -1*vel_w_s*JxW[qp]*dpsi[i][qp](2);   
   
   
  }

#endif



#if FLUID_P_CONST

//new kCl ITERATION

Number Penalty=1;
          for (unsigned int i=0; i<n_p_dofs; i++){
          for (unsigned int j=0; j<n_p_dofs; j++){

              Kpp(j,i) += Penalty*JxW[qp]*psi[j][qp]*psi[i][qp];
          }
	}


for (unsigned int i=0; i<n_p_dofs; i++){
       Fp(i) += Penalty*((lambda+lambda)/2.0)*JxW[qp]*psi[i][qp];
}

#if !KCL
    for (unsigned int i=0; i<n_u_dofs; i++){
          for (unsigned int j=0; j<n_m_dofs; j++){
           //the velocity divergence term 

              Kmu(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
              Kmv(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              Kmw(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);   

//std::cout<<" Kmu(j,i) " << Kmu(j,i)<<std::endl;
          }
	}




  for (unsigned int i=0; i<n_m_dofs; i++){
  for (unsigned int j=0; j<n_m_dofs; j++){
   
   Kmm(i,j) += (1.0/(J*dt))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
    }    
  }


for (unsigned int i=0; i<n_m_dofs; i++){
Fm(i) += (1.0/(J*dt))*m_olD*JxW[qp]*psi[i][qp];
}
#endif


#if KCL&& FLUID_P_CONST
  for (unsigned int i=0; i<n_m_dofs; i++){
  for (unsigned int j=0; j<n_m_dofs; j++){
   Kmm(i,j) += (1.0/(J*dt))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
    }    
  }


for (unsigned int i=0; i<n_m_dofs; i++){
Fm(i) += (1.0/(J*dt))*m_olD*JxW[qp]*psi[i][qp];
Fm(i) += Kperm*lambda*JxW[qp]*dpsi[i][qp](0)*dpsi[i][qp](0);
Fm(i) += Kperm*lambda*JxW[qp]*dpsi[i][qp](1)*dpsi[i][qp](1);
Fm(i) += Kperm*lambda*JxW[qp]*dpsi[i][qp](2)*dpsi[i][qp](2);
}

#endif

#endif
=======
        Fp(i) +=  (p_old/(J*M*fchap))*JxW[qp]*psi[i][qp];
        Fp(i) += -dt*div_vs*JxW[qp]*psi[i][qp];

  for (unsigned int j=0; j<n_p_dofs; j++){
    Kpp(i,j) += (1/(J*M*fchap))*JxW[qp]*psi[i][qp]*(psi[j][qp]);
    Kpp(i,j) += -dt*(fchapd/(M*pow(fchap,2.0)))*div_vs*JxW[qp]*psi[i][qp]*(psi[j][qp]);
    }    
  }

Real q=0;
Point point = elem->centroid();

if ( (point(0)>0.8)    ){
q= 0.0;
}
            for (unsigned int j=0; j<n_u_dofs; j++){
             // Fu(j) += JxW[qp]*dphi[j][qp](0)*p_fluid;
            }
             for (unsigned int j=0; j<n_v_dofs; j++){
             // Fv(j) += JxW[qp]*dphi[j][qp](1)*p_fluid;
            }
             for (unsigned int j=0; j<n_w_dofs; j++){
             // Fw(j) += JxW[qp]*dphi[j][qp](2)*p_fluid;
          }
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76

        } // end of the quadrature point qp-loop


<<<<<<< HEAD
    dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      } // end of element loop
      std::cout<<" ok "<<std::endl;

     assemble_fluid_bcs(es);
        std::cout<<" nah "<<std::endl;
=======


        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      } // end of element loop

     assemble_fluid_bcs(es);
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76


    return;
  }

#endif FLUID_VEL



