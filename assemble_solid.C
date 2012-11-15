#include "defines.h"
#include "assemble.h"
#include "neohooke_cc.h"
#include "poro_elastic_cc.h"
#include "mooney_cc.h"
#include "anal_neo_cc.h"

#define INERTIA 1

#define DT 0

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
  
  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & newton_update =
   es.get_system<TransientLinearImplicitSystem> ("Newton-update");

  // New
   TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

 TransientLinearImplicitSystem & fluid_system_vel =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system-vel");

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



#if FLUID_P_CONST
    const unsigned int m_var = fluid_system_vel.variable_number ("fluid_M");
  std::vector<unsigned int> dof_indices_m;
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

#if FLUID_P_CONST 
  const DofMap & dof_map_fluid = fluid_system_vel .get_dof_map();
#endif

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
    Fu(Fe), Fv(Fe), Fw(Fe);
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
   

#if INERTIA
test(67);
  const unsigned int a_var = last_non_linear_soln.variable_number ("a");
  const unsigned int b_var = last_non_linear_soln.variable_number ("b");
  const unsigned int c_var = last_non_linear_soln.variable_number ("c");

//B block
  DenseSubMatrix<Number>
  Kua(Ke), Kub(Ke), Kuc(Ke),
  Kva(Ke), Kvb(Ke), Kvc(Ke),
  Kwa(Ke), Kwb(Ke), Kwc(Ke); 

//C block
  DenseSubMatrix<Number>
  Kau(Ke), Kav(Ke), Kaw(Ke),
  Kbu(Ke), Kbv(Ke), Kbw(Ke),
  Kcu(Ke), Kcv(Ke), Kcw(Ke);

//D block
  DenseSubMatrix<Number>
  Kaa(Ke), Kab(Ke), Kac(Ke),
  Kba(Ke), Kbb(Ke), Kbc(Ke),
  Kca(Ke), Kcb(Ke), Kcc(Ke);

  DenseSubVector<Number>
  Fa(Fe), Fb(Fe), Fc(Fe);

  std::vector<unsigned int> dof_indices_a;
  std::vector<unsigned int> dof_indices_b;
  std::vector<unsigned int> dof_indices_c;
test(68);
#endif

    DenseMatrix<Real> stiff;
  DenseVector<Real> res;
  VectorValue<Gradient> grad_u_mat;

  VectorValue<Gradient> grad_u_mat_old;
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");


#if PORO 
  DenseVector<Real> p_stiff;
  DenseVector<Real> p_res;
  PoroelasticConfig material(dphi,psi);
#endif

  // Just calculate jacobian contribution when we need to
  material.calculate_linearized_stiffness = true;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

  for ( ; el != end_el; ++el)
    {  
test(69);  
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

#if FLUID_P_CONST
      dof_map_fluid.dof_indices (elem, dof_indices_m, m_var);
#endif
      //elem->print_info();

      fe_vel->reinit  (elem);

#if INCOMPRESSIBLE
      fe_pres->reinit (elem);
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




#if INERTIA

//B block
   Kua.reposition (u_var*n_u_dofs, 3*n_u_dofs + n_p_dofs, n_u_dofs, n_u_dofs);
   Kub.reposition (u_var*n_u_dofs, 4*n_u_dofs + n_p_dofs, n_u_dofs, n_v_dofs);
   Kuc.reposition (u_var*n_u_dofs, 5*n_u_dofs + n_p_dofs, n_u_dofs, n_w_dofs);
   Kva.reposition (v_var*n_v_dofs, 3*n_u_dofs + n_p_dofs, n_v_dofs, n_u_dofs);
   Kvb.reposition (v_var*n_v_dofs, 4*n_u_dofs + n_p_dofs, n_v_dofs, n_v_dofs);
   Kvc.reposition (v_var*n_v_dofs, 5*n_u_dofs + n_p_dofs, n_v_dofs, n_w_dofs);
   Kwa.reposition (w_var*n_w_dofs, 3*n_u_dofs + n_p_dofs, n_w_dofs, n_u_dofs);
   Kwb.reposition (w_var*n_w_dofs, 4*n_u_dofs + n_p_dofs, n_w_dofs, n_v_dofs);
   Kwc.reposition (w_var*n_w_dofs, 5*n_u_dofs + n_p_dofs, n_w_dofs, n_w_dofs);

test(701);  
//C block
   Kau.reposition (3*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
   Kav.reposition (3*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
   Kaw.reposition (3*n_u_dofs + n_p_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
   Kbu.reposition (4*n_u_dofs + n_p_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
   Kbv.reposition (4*n_u_dofs + n_p_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
   Kbw.reposition (4*n_u_dofs + n_p_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
   Kcu.reposition (5*n_u_dofs + n_p_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
   Kcv.reposition (5*n_u_dofs + n_p_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
   Kcw.reposition (5*n_u_dofs + n_p_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);

//D block
   Kaa.reposition (3*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_u_dofs, n_u_dofs);
   Kab.reposition (3*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_u_dofs, n_v_dofs);
   Kac.reposition (3*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_u_dofs, n_w_dofs);
   Kba.reposition (4*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_v_dofs, n_u_dofs);
   Kbb.reposition (4*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_v_dofs, n_v_dofs);
   Kbc.reposition (4*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_v_dofs, n_w_dofs);
   Kca.reposition (5*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_w_dofs, n_u_dofs);
   Kcb.reposition (5*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_w_dofs, n_v_dofs);
   Kcc.reposition (5*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_w_dofs, n_w_dofs);


Fa.reposition (3*n_u_dofs + n_p_dofs, n_u_dofs);
Fb.reposition (4*n_u_dofs + n_p_dofs, n_v_dofs);
Fc.reposition (5*n_u_dofs + n_p_dofs, n_w_dofs);

  dof_map.dof_indices (elem, dof_indices_a, a_var);
  dof_map.dof_indices (elem, dof_indices_b, b_var);
  dof_map.dof_indices (elem, dof_indices_c, c_var);

test(71);  
#endif


      System& aux_system = es.get_system<System>("Reference-Configuration");
      std::vector<unsigned int> undefo_index;
      std::vector<unsigned int> vel_index;           

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {



 Point rX;
 for (unsigned int l=0; l<n_u_dofs; l++)
            {
rX(0) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
rX(1) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_v[l]);
rX(2) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_w[l]);
            }



#if INERTIA || DT
test(72);  
Real rho_s=15;

Point current_x;
 for (unsigned int l=0; l<n_u_dofs; l++)
 {
current_x(0) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_u[l]);
current_x(1) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_v[l]);
current_x(2) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_w[l]);
}

Point old_x;
 for (unsigned int l=0; l<n_u_dofs; l++)
 {
old_x(0) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_u[l]);
old_x(1) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_v[l]);
old_x(2) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_w[l]);
}
#if INERTIA
Point old_vel;
 for (unsigned int l=0; l<n_u_dofs; l++)
 {
old_vel(0) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_a[l]);
old_vel(1) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_b[l]);
old_vel(2) += phi[l][qp]*last_non_linear_soln.old_local_solution->el(dof_indices_c[l]);
}
Point current_vel;
 for (unsigned int l=0; l<n_u_dofs; l++)
 {
current_vel(0) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_a[l]);
current_vel(1) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_b[l]);
current_vel(2) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_c[l]);
}
#endif

#if UN_MINUS_ONE
Point unm1_x;
 for (unsigned int l=0; l<n_u_dofs; l++)
 {
unm1_x(0) += phi[l][qp]*unm1.old_local_solution->el(dof_indices_u[l]);
unm1_x(1) += phi[l][qp]*unm1.old_local_solution->el(dof_indices_v[l]);
unm1_x(2) += phi[l][qp]*unm1.old_local_solution->el(dof_indices_w[l]);
}
#endif

Point value_acc;
Point value_acc_alt;

#if DT
for (unsigned int d = 0; d < dim; ++d) {
 value_acc_alt(d) = (rho_s)*( ((current_x(d)-rX(d))-(old_x(d)-rX(d)))-((old_x(d)-rX(d))- (unm1_x(d)-rX(d))) );  
value_acc(d) = (rho_s)*((current_x(d))-2*(old_x(d))+ (unm1_x(d)));  
value_acc(d) = (rho_s)*((current_x(d))-(old_x(d)));  
}
#endif

Point res_1;
Point res_2;
#if INERTIA
for (unsigned int d = 0; d < dim; ++d) {
res_1(d) = (rho_s)*((current_vel(d))-(old_vel(d)));
res_2(d) = current_x(d)-dt*current_vel(d)-old_x(d);    
}
/*
std::cout<<" current_vel "<<current_vel<<std::endl;
std::cout<<" res_1 "<<res_1<<std::endl;
std::cout<<" res_2 "<<res_2<<std::endl;
*/
#endif



test(73);  
#endif


Number   p_solid = 0.;

#if MOVING_MESH
	grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
    for (unsigned int d = 0; d < dim; ++d) {
      std::vector<Number> u_undefo;
//Fills the vector di with the global degree of freedom indices for the element. :dof_indicies
      aux_system.get_dof_map().dof_indices(elem, undefo_index,d);
      aux_system.current_local_solution->get(undefo_index, u_undefo);
      for (unsigned int l = 0; l != n_u_dofs; l++)
        grad_u_mat(d).add_scaled(dphi[l][qp], u_undefo[l]); 
 }
#endif


//#include "fixed_mesh_in_solid_assemble_code.txt"
          
      #if INCOMPRESSIBLE
      for (unsigned int l=0; l<n_p_dofs; l++)
            {
              p_solid += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
            }
      #endif
	  


#if INCOMPRESSIBLE 
Real m=0;
Real p_fluid=0;

#if FLUID_VEL 
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


Real m_old=0;

#if FLUID_P_CONST
for (unsigned int l=0; l<n_p_dofs; l++)
 {
   m += psi[l][qp]*fluid_system_vel.current_local_solution->el(dof_indices_m[l]);
 }


for (unsigned int l=0; l<n_p_dofs; l++)
 {
   m_old += psi[l][qp]*fluid_system_vel.old_local_solution->el(dof_indices_m[l]);
 }
#endif

material.init_for_qp(grad_u_mat, p_solid, qp, 1.0*m+0.0*m_old, p_fluid);

#endif
#endif //#if INCOMPRESSIBLE


#if INCOMPRESSIBLE && ! PORO
material.init_for_qp(grad_u_mat, p_solid, qp);
#endif 

          for (unsigned int i=0; i<n_u_dofs; i++)
            {
            res.resize(dim);
            material.get_residual(res, i);
            res.scale(JxW[qp]);
#if INERTIA
            res.scale(dt);
#endif

#if DT
            res.scale(dt);
#endif
//std::cout<< "res "<<res<<std::endl;

      	    Fu(i) += res(0);              
            Fv(i) += res(1) ; 
	    Fw(i) += res(2);  

  	#if GRAVITY
        Real grav=0.0;
        Fu(i) += progress*grav*JxW[qp]*phi[i][qp];
	#endif

#if INERTIA
      Fu(i) +=  JxW[qp]*phi[i][qp]*res_1(0); 
      Fv(i) +=  JxW[qp]*phi[i][qp]*res_1(1);     
      Fw(i) +=  JxW[qp]*phi[i][qp]*res_1(2); 

      Fa(i) +=  JxW[qp]*phi[i][qp]*res_2(0);  
      Fb(i) +=  JxW[qp]*phi[i][qp]*res_2(1);      
      Fc(i) +=  JxW[qp]*phi[i][qp]*res_2(2);  
#endif


// Matrix contributions for the uu and vv couplings.
for (unsigned int j=0; j<n_u_dofs; j++)
   {
    material.get_linearized_stiffness(stiff, i, j);
    stiff.scale(JxW[qp]);

#if DT
            res.scale(dt);
#endif

#if INERTIA 
    stiff.scale(dt);
    Kua(i,j)+=  rho_s*JxW[qp]*phi[i][qp]*phi[j][qp];      
    Kvb(i,j)+=  rho_s*JxW[qp]*phi[i][qp]*phi[j][qp];
    Kwc(i,j)+=  rho_s*JxW[qp]*phi[i][qp]*phi[j][qp];


    Kau(i,j)+=  JxW[qp]*phi[i][qp]*phi[j][qp];      
    Kbv(i,j)+=  JxW[qp]*phi[i][qp]*phi[j][qp];
    Kcw(i,j)+=  JxW[qp]*phi[i][qp]*phi[j][qp];

    Kaa(i,j)+=  -dt*JxW[qp]*phi[i][qp]*phi[j][qp];      
    Kbb(i,j)+=  -dt*JxW[qp]*phi[i][qp]*phi[j][qp];
    Kcc(i,j)+=  -dt*JxW[qp]*phi[i][qp]*phi[j][qp];
#endif




    Kuu(i,j)+=  stiff(u_var, u_var);
    Kuv(i,j)+=  stiff(u_var, v_var);
    Kuw(i,j)+=  stiff(u_var, w_var);	      
    Kvu(i,j)+=  stiff(v_var, u_var);
    Kvv(i,j)+=  stiff(v_var, v_var);
    Kvw(i,j)+=  stiff(v_var, w_var);
    Kwu(i,j)+=  stiff(w_var, u_var);
    Kwv(i,j)+=  stiff(w_var, v_var);
    Kww(i,j)+=  stiff(w_var, w_var); 


#if GRAVITY
    Kuu(i,j)+= 1*JxW[qp]*phi[i][qp]*phi[j][qp];
#endif
                }
            }


#if INCOMPRESSIBLE && FLUID_P_CONST
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

#if CHAP && ! FLUID_P_CONST
           for (unsigned int i = 0; i < n_p_dofs; i++) {
         Fp(i) += 0*JxW[qp]*psi[i][qp];
    }
    
    for (unsigned int i = 0; i < n_p_dofs; i++) {
          for (unsigned int j = 0; j < n_p_dofs; j++) {
            Kpp(i, j) += 1*JxW[qp]*psi[i][qp]*psi[j][qp];
    }
    }
#endif



}//end of qp loop

      newton_update.matrix->add_matrix (Ke, dof_indices);
      newton_update.rhs->add_vector    (Fe, dof_indices);
} // end of element loop

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

std::ofstream lhs_out("lhsoutS3.dat");
Ke.print(lhs_out);
lhs_out.close();
  return;
}




