#include "defines.h"
#include "assemble.h"
#include "neohooke_cc.h"
#include "poro_elastic_cc.h"
#include "mooney_cc.h"

//#include "solid_system.h"


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_solid (EquationSystems& es,
                      const std::string& system_name)
{

//es.print_info();

#if LOG_ASSEMBLE_PERFORMANCE
  PerfLog perf_log("Assemble");
  perf_log.push("assemble stiffness");
#endif
    // Get a reference to the auxiliary system
  //TransientExplicitSystem& aux_system = es.get_system<TransientExplicitSystem>("Newton-update");

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert (system_name == "Newton-update");
  
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  test(61);
  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & newton_update =
   es.get_system<TransientLinearImplicitSystem> ("Newton-update");


  // New
   TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

#if FLUID && PORO
 TransientLinearImplicitSystem & fluid_system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system");
#endif

#if FLUID_VEL && PORO
 TransientLinearImplicitSystem & fluid_system_vel =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system-vel");
#endif


#if VELOCITY
TransientLinearImplicitSystem&  velocity = es.get_system<TransientLinearImplicitSystem>("velocity-system");
#endif

#if UN_MINUS_ONE
TransientLinearImplicitSystem & unm1 =
    es.get_system<TransientLinearImplicitSystem> ("unm1-system");
#endif
test(62);
const System & ref_sys = es.get_system("Reference-Configuration"); 
test(63);
  
  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = last_non_linear_soln .variable_number ("u");
  const unsigned int v_var = last_non_linear_soln .variable_number ("v");
  const unsigned int w_var = last_non_linear_soln .variable_number ("w");
#if INCOMPRESSIBLE
  const unsigned int p_var = last_non_linear_soln .variable_number ("p");
#endif 
  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = last_non_linear_soln.variable_type(u_var);


test(64);

#if INCOMPRESSIBLE
  // Get the Finite Element type for "p".
  FEType fe_pres_type = last_non_linear_soln .variable_type(p_var);
#endif 

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  
#if INCOMPRESSIBLE 
  // Build a Finite Element object of the specified type for
  // the pressure variables.
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
#endif 
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
test(65);
//        AutoPtr<QBase> qrule2(fe_vel_type.default_quadrature_rule(dim));
// fe_vel->attach_quadrature_rule (qrule2.get());

#if INCOMPRESSIBLE 
  fe_pres->attach_quadrature_rule (&qrule);
#endif
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
test(66);
#if INCOMPRESSIBLE 
  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
#endif
  
 const std::vector<Point>& coords = fe_vel->get_xyz();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = last_non_linear_soln .get_dof_map();

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
  // Find out what the timestep size parameter is from the system, and
  // the value of theta for the theta method.  We use implicit Euler (theta=1)
  // for this simulation even though it is only first-order accurate in time.
  // The reason for this decision is that the second-order Crank-Nicolson
  // method is notoriously oscillatory for problems with discontinuous
  // initiaFl data such as the lid-driven cavity.  Therefore,
  // we sacrifice accuracy in time for stability, but since the solution
  // reaches steady state relatively quickly we can afford to take small
  // timesteps.  If you monitor the initial nonlinear residual for this
  // simulation, you should see that it is monotonically decreasing in time.
 // const Real dt    = es.parameters.get<Real>("dt");
  // const Real time  = es.parameters.get<Real>("time");
 // const Real theta = 1.;
    
    DenseMatrix<Real> stiff;
  DenseVector<Real> res;
  VectorValue<Gradient> grad_u_mat;

#if DYNAMIC
  VectorValue<Gradient> grad_u_mat_old;
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");
#endif

#if PORO 
  DenseVector<Real> p_stiff;
  DenseVector<Real> p_res;
  //Real m = 0.0*progress;
  PoroelasticConfig material(dphi,psi);
#endif


#if COMPRESSIBLE && NEO
    NeoHookeCurrentConfig material(dphi);    
#endif

#if INCOMPRESSIBLE && NEO
  DenseVector<Real> p_stiff;
  DenseVector<Real> p_res;
  NeoHookeCurrentConfig material(dphi,psi);
#endif

#if COMPRESSIBLE && MOONEY
  MooneyCurrentConfig material(dphi);
#endif
    
#if INCOMPRESSIBLE && MOONEY
  DenseVector<Real> p_stiff;
  DenseVector<Real> p_res;
  MooneyCurrentConfig material(dphi,psi);
#endif

  // Just calculate jacobian contribution when we need to
  material.calculate_linearized_stiffness = true;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

  for ( ; el != end_el; ++el)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
      
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
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
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
	
      //elem->print_info();

      fe_vel->reinit  (elem);

#if INCOMPRESSIBLE
      fe_pres->reinit (elem);
#endif

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

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

      System& aux_system = es.get_system<System>("Reference-Configuration");
      std::vector<unsigned int> undefo_index;
      #if DYNAMIC
      std::vector<unsigned int> vel_index;
      #endif
           

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

	  #if INCOMPRESSIBLE
	  Number   p = 0.;
	  #endif

#if MOVING_MESH
	grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> u_undefo;
//Fills the vector di with the global degree of freedom indices for the element. :dof_indicies
      aux_system.get_dof_map().dof_indices(elem, undefo_index,d);
      aux_system.current_local_solution->get(undefo_index, u_undefo);
      for (unsigned int l = 0; l != n_u_dofs; l++)
        grad_u_mat(d).add_scaled(dphi[l][qp], u_undefo[l]); 
//std::cout << "u_undefo" << u_undefo[0]<< std::endl;
//std::cout << "dphi[l][qp] " << dphi[1][1] << std::endl;
    }
#endif
   // std::cout << "grad_u_mat" << grad_u_mat(0)<< std::endl;


#if FIXED_MESH
	grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> X_undefo;
      std::vector<Number> X_disp;

      last_non_linear_soln.get_dof_map().dof_indices(elem, undefo_index,d);
      last_non_linear_soln.current_local_solution->get(undefo_index, X_disp);

      for (unsigned int l = 0; l != n_u_dofs; l++)
        grad_u_mat(d).add_scaled(dphi[l][qp], X_disp[l]); // u_current(l)); // -
    }
#endif


#if DYNAMIC & FIXED_MESH
grad_u_mat_old(0) = grad_u_mat_old(1) = grad_u_mat_old(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> X_undefo;
      std::vector<Number> X_disp;
      last_non_linear_soln.get_dof_map().dof_indices(elem, undefo_index,d);
      last_non_linear_soln.old_local_solution->get(undefo_index, X_disp);
      for (unsigned int l = 0; l != n_u_dofs; l++)
        grad_u_mat_old(d).add_scaled(dphi[l][qp],X_disp[l]); 
    }
grad_u_mat(0)=0.5*(grad_u_mat_old(0)+grad_u_mat(0));
grad_u_mat(1)=0.5*(grad_u_mat_old(1)+grad_u_mat(1));
grad_u_mat(2)=0.5*(grad_u_mat_old(2)+grad_u_mat(2));
#endif
          
            #if INCOMPRESSIBLE
            // Compute the current pressure value at this quadrature point.
          for (unsigned int l=0; l<n_p_dofs; l++)
            {
              p += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
            }
	  #endif
	  


#if INCOMPRESSIBLE &&  PORO
Real m=0;
Real p_fluid=0;

#if FLUID
for (unsigned int l=0; l<n_u_dofs; l++)
            {
              m += phi[l][qp]*fluid_system.current_local_solution->el(dof_indices_u[l]);
            }
#endif


#if FLUID_VEL && CHAP
for (unsigned int l=0; l<n_p_dofs; l++)
 {
              p_fluid += psi[l][qp]*fluid_system_vel.current_local_solution->el(dof_indices_p[l]);
 }


//As outlined in Chappel p=(p_curr-p_old)/2
 Real p_fluid_old=0;
for (unsigned int l=0; l<n_p_dofs; l++)
 {
              p_fluid_old += psi[l][qp]*fluid_system_vel.old_local_solution->el(dof_indices_p[l]);
 }

p_fluid=0.5*p_fluid+0.5*p_fluid_old;

#if DECOUPLE
p_fluid=0;
#endif

material.init_for_qp(grad_u_mat, p, qp, m, p_fluid);
#endif

#if !CHAP
material.init_for_qp(grad_u_mat, p, qp, m,p_fluid);
#endif

#endif 


#if INCOMPRESSIBLE && ! PORO
material.init_for_qp(grad_u_mat, p, qp);
#endif 

#if COMPRESSIBLE 
Number p_comp=0;
material.init_for_qp(grad_u_mat,p_comp,qp);
#endif
          for (unsigned int i=0; i<n_u_dofs; i++)
            {
            res.resize(dim);
            material.get_residual(res, i);
            res.scale(JxW[qp]);

       		//Real E=10;
       		//Real nu=0.3;
       		//Real mu = E / (2.0 * (1.0 + nu));
       		//Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

      	Fu(i) += res(0);              
        Fv(i) += res(1) ; 
	      Fw(i) += res(2);  

  #if GRAVITY
  #if STATIC
        Fw(i) += 0.1*JxW[qp]*phi[i][qp];
  #endif       
  #if DYNAMIC
        Real grav=0.0;
	   // Fw(i) += progress*0.1*JxW[qp]*phi[i][qp];
           Fu(i) += progress*grav*JxW[qp]*phi[i][qp];
  #endif        
#endif
	   
 #if DYNAMIC
//Just add (rho/(dt^2))*Mass to the Jacobian matrix
    const Real rho_mix =0.0;
    const Real fac=(2*rho_mix)/(dt*dt);
      for (unsigned int j=0; j<n_u_dofs; j++){
	Real value = 1*fac*JxW[qp]*phi[i][qp]*phi[j][qp];
          Kuu(i,j)+= value;
          Kvv(i,j)+=  value;
          Kww(i,j)+=  value;
      }
 #endif


   
#if DYNAMIC

   std::vector<Number> value_acc(3,0);
    for (unsigned int d = 0; d < dim; ++d) {

   std::vector<Number> old_x;
   std::vector<Number> current_x;
   std::vector<Number> X;

   newton_update.get_dof_map().dof_indices(elem, undefo_index,0);
   last_non_linear_soln.old_local_solution->get(undefo_index, old_x);
   last_non_linear_soln.current_local_solution->get(undefo_index, current_x);
   ref_sys.current_local_solution->get(undefo_index, X);

  #if VELOCITY
  std::vector<Number> velocity_d;
  velocity.current_local_solution->get(undefo_index, velocity_d);
  value_acc[d] = fac*JxW[qp]*phi[i][qp]*((current_x[i]-X[i])-(old_x[i]-X[i])-dt*velocity_d[i]);    
  #endif
 

#if UN_MINUS_ONE
  std::vector<Number> unm1_x;
  unm1.old_local_solution->get(undefo_index, unm1_x);
  value_acc[d] = fac*JxW[qp]*phi[i][qp]*( (current_x[i]-X[i])-(old_x[i]-X[i])-1*((old_x[i]-X[i])-(unm1_x[i]-X[i])) );   //Might want to also change the jaconiamn for this since we have an extra current_x (so multiply M by 0.5, I think). 
#endif

}
      Fu(i) +=  value_acc[0]; 
      Fv(i) +=  value_acc[1];     
      Fw(i) +=  value_acc[2];         
#endif


              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j<n_u_dofs; j++)
                {
                      material.get_linearized_stiffness(stiff, i, j);
		      stiff.scale(JxW[qp]);

		      Kuu(i,j)+=  stiff(u_var, u_var);
#if GRAVITY
                    Kuu(i,j)+= 1*JxW[qp]*phi[i][qp]*phi[j][qp];
#endif

		      Kuv(i,j)+=  stiff(u_var, v_var);
		      Kuw(i,j)+=  stiff(u_var, w_var);
		      
		      Kvu(i,j)+=  stiff(v_var, u_var);
		      Kvv(i,j)+=  stiff(v_var, v_var);
		      Kvw(i,j)+=  stiff(v_var, w_var);

		      Kwu(i,j)+=  stiff(w_var, u_var);
		      Kwv(i,j)+=  stiff(w_var, v_var);
		      Kww(i,j)+=  stiff(w_var, w_var); 
                }
            }

#if INCOMPRESSIBLE && !CHAP
           for (unsigned int i = 0; i < n_p_dofs; i++) {
	  material.get_p_residual(p_res, i);
	  p_res.scale(JxW[qp]);
          Fp(i) += p_res(0);
	  }
    
    for (unsigned int i = 0; i < n_u_dofs; i++) {
          for (unsigned int j = 0; j < n_p_dofs; j++) {
	    material.get_linearized_uvw_p_stiffness(p_stiff, i, j);
	    p_stiff.scale(JxW[qp]);
            Kup(i, j) += p_stiff(0);
	          Kvp(i, j) += p_stiff(1);
            Kwp(i, j) += p_stiff(2);
	  }
    }
    
    for (unsigned int i = 0; i < n_p_dofs; i++) {
          for (unsigned int j = 0; j < n_u_dofs; j++) {
	    material.get_linearized_p_uvw_stiffness(p_stiff, i, j);
	    p_stiff.scale(JxW[qp]);
      Kpu(i, j) += p_stiff(0);
	    Kpv(i, j) += p_stiff(1);
      Kpw(i, j) += p_stiff(2);
       }
    }
#endif

#if CHAP
           for (unsigned int i = 0; i < n_p_dofs; i++) {
         Fp(i) += 0*JxW[qp]*psi[i][qp];
    }
    
    for (unsigned int i = 0; i < n_p_dofs; i++) {
          for (unsigned int j = 0; j < n_p_dofs; j++) {
            Kpp(i, j) += JxW[qp]*psi[i][qp]*psi[j][qp];
    }
    }

#endif



}//end of qp loop

      newton_update.matrix->add_matrix (Ke, dof_indices);
      newton_update.rhs->add_vector    (Fe, dof_indices);
} // end of element loop
test(72);
     //   dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
     newton_update.rhs->close();
     newton_update.matrix->close();


#if LOG_ASSEMBLE_PERFORMANCE
perf_log.pop("assemble stiffness");
#endif 

#if LOG_ASSEMBLE_PERFORMANCE
perf_log.push("assemble bcs");
#endif

//Assemble the boundary conditions.
assemble_bcs(es);

#if LOG_ASSEMBLE_PERFORMANCE
perf_log.pop("assemble bcs");
#endif


  return;
}




