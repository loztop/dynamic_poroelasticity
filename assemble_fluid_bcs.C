#include "defines.h"
#include "assemble.h"
#include <cmath>


using namespace std;
#define PI 3.14159265
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_fluid_bcs (EquationSystems& es)
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
    // AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_vel_type));
    // QGauss qface (dim-1, fe_vel_type.default_quadrature_order());
    // const std::vector<Real>& JxW_face = fe_face->get_JxW();
    // const std::vector<std::vector<Real> >& psi_face = fe_face->get_phi();	
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
test(68);


  

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
#if NEUMANN_PRESSURE

  const std::vector<std::vector<Real> >&  phi_face =  fe_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
  
  
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Point>& face_normals = fe_face->get_normals();
  fe_face->reinit(elem,s);  

  const std::vector<std::vector<Real> >&  psi_face =  fe_face_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi_face = fe_face_pres->get_dphi();
  fe_face_pres->reinit(elem,s);  


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

#if FLUID_P_CONST
 for (unsigned int qp=0; qp<qface->n_points(); qp++)
{  
Number   lambda = 0.;

for (unsigned int l=0; l<psi_face.size(); l++)
{
              lambda += psi_face[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);

//std::cout<< "last_non_linear_soln.current_local_solution->el(dof_indices_p[l]) "<<last_non_linear_soln.current_local_solution->el(dof_indices_p[l]) <<std::endl;

}
Real val=lambda;
//std::cout<< val<<std::endl;

for (unsigned int i=0; i<psi_face.size(); i++){
Fu(i) +=  -JxW_face[qp]*val*face_normals[qp](0)*phi_face[i][qp];
Fv(i) +=  -JxW_face[qp]*val*face_normals[qp](1)*phi_face[i][qp];
Fw(i) +=  -JxW_face[qp]*val*face_normals[qp](2)*phi_face[i][qp];
}
}  
#endif

    if ((elem->node(n) == side->node(ns)) && (p(0)<0.01)   )

    for (unsigned int qp=0; qp<qface->n_points(); qp++)
    {       
  if ((elem->node(n) == side->node(ns)) && (p(0)<0.01) && (face_normals[qp](0) < -0.99)  )
{
   Real value=0;
    
      //value = 0.1*(1-exp(-pow(time,2.0)/0.25)) ;
      value = 0.0*progress ;

     #if CHAP_SWELL

   Real factor=0.5;

      value = -factor*1000*(1-exp(-pow(((time+0.00000001)/10.0),2.0)/0.25)) ;
      //value = -500.0*progress ;


     #endif

      //value = 1;
      for (unsigned int i=0; i<phi_face.size(); i++){
        
        #if MASS_NEUMANN_PRESSURE //In Chapelle there is no neumann pressure Bc in the mass concservation eqn.
      //  Fp(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        #endif

       

//std::cout<< "p " << p<<std::endl;

//std::cout<< face_normals[qp] <<std::endl;

        #if MOMENTUM_NEUMANN_PRESSURE
        Fu(i) +=  JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        Fv(i) +=  JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
        Fw(i) +=  JxW_face[qp]*value*face_normals[qp](2)*phi_face[i][qp];
        #endif

        

        }
    
    #if KCL
Number   lambda = 0.;
for (unsigned int l=0; l<n_p_dofs; l++)
{
              lambda += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
}

      for (unsigned int i=0; i<psi_face.size(); i++){
        value = lambda
        Fm(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        Fm(i) += - JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
        Fm(i) += - JxW_face[qp]*value*face_normals[qp](2)*phi_face[i][qp];
        } 
        #endif


      for (unsigned int i=0; i<phi_face.size(); i++){
                      for (unsigned int j=0; j<phi_face.size(); j++){
 // Kuu(i,j) += value*face_normals[qp](0)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
 // Kvv(i,j) += value*face_normals[qp](1)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
 // Kww(i,j) += value*face_normals[qp](2)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];

  //  Kuu(i,j) += JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
 // Kvv(i,j) += JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  //Kww(i,j) += JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
          }
        } 
    } //end qp



  } //end if 


/*
  if ((elem->node(n) == side->node(ns)) && (p(2)<0.01)  )
{
    for (unsigned int qp=0; qp<qface->n_points(); qp++)
    {       

   Real value=time*0.00;

      for (unsigned int i=0; i<phi_face.size(); i++){

        #if MOMENTUM_NEUMANN_PRESSURE
        Fu(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        Fv(i) += - JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
        Fw(i) += - JxW_face[qp]*value*face_normals[qp](2)*phi_face[i][qp];
        #endif

      } 
    } //end qp
  } //end if 
  */


  }    
}
#endif //end NEUMANN_PRESSURE


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
    
//Apply Dirichlet Bcs properly
#if CUBE 
  if ((elem->node(n) == side->node(ns))  )
  {
  #if DIRICHLET_VELOCITY

#if SEALED_CUBE

    
if ((elem->node(n) == side->node(ns)) && (p(0)>1.49)  )
    {
      unsigned int source_dof = node->dof_number(system.number(), 0, 0);
       Real value = 0;
       rows.push_back(source_dof);
       system.rhs->set(source_dof,value);
    }  //end if
    

if ((elem->node(n) == side->node(ns)) && ( (p(1)<0.01)  || (p(1)>1.49))  )
    {

       unsigned int source_dof = node->dof_number(system.number(), 1, 0);
       Real value = 0;
       rows.push_back(source_dof);
       system.rhs->set(source_dof,value);
     

    }  //end if

if ((elem->node(n) == side->node(ns)) && ((p(2)<0.01)  || (p(2)>1.49))  )
    {


      unsigned int source_dof = node->dof_number(system.number(), 2, 0);
       Real value = 0;
       rows.push_back(source_dof);
      system.rhs->set(source_dof,value);

    

    }  //end if
    
    
#endif

    
    #if DISK_FLOW
    if ((elem->node(n) == side->node(ns)) && (  (p(0)<0.01 ) ) )
    {
     for (unsigned int d = 1; d < 3; ++d) {
        unsigned int source_dof = node->dof_number(system.number(), d, 0);
       Real value = 0;
       rows.push_back(source_dof);
       system.rhs->set(source_dof,value);
    }  //end dimension loop

//std::cout<<"Disk flow"<<std::endl;
    
       unsigned int source_dof = node->dof_number(system.number(), 0, 0);
    Real value = 0.5* ( (0.75*0.75-(p(1)-0.75)*(p(1)-0.75)) + (0.75*0.75-(p(2)-0.75)*(p(2)-0.75)) ) ;
value = 0.015* ( (0.75*0.75-(p(1)-0.75)*(p(1)-0.75)) + (0.75*0.75-(p(2)-0.75)*(p(2)-0.75)) )*(1-exp(-pow(time,2.0)/0.25)) ; 


value = 0.1* ( (0.75*0.75-(p(1)-0.75)*(p(1)-0.75)) + (0.75*0.75-(p(2)-0.75)*(p(2)-0.75)) ) ; 

if (progress<10.2){
value = 0.10* ( (0.75*0.75-(p(1)-0.75)*(p(1)-0.75)) + (0.75*0.75-(p(2)-0.75)*(p(2)-0.75)) ) ;
}else{
value = 0.0;
}
    rows.push_back(source_dof);
      system.rhs->set(source_dof,value);
       
    }  //end if
   #endif 
   //DISK_FLOW


/*
if ((elem->node(n) == side->node(ns)) && (p(1)<0.0001)  )
    {
       unsigned int source_dof = node->dof_number(system.number(), 1, 0);
       Real value = 0;
       rows.push_back(source_dof);
       system.rhs->set(source_dof,value);

    }  //end if
  */
    
#endif 


#if DIRICHLET_PRESSURE

#if CHAP_SWELL


    if ((elem->node(n) == side->node(ns)) && (p(0)>1.4999 ) )
    {
       unsigned int source_dof = node->dof_number(system.number(), 3, 0);
       if (source_dof <12345678){ //The pressures do not exist everywhere// This is a hack !!
       Real value = 0;
       pressure_rows.push_back(source_dof);
      system.rhs->set(source_dof,value);
     }
    }  //end if


/*
    if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 ) )
    {
       unsigned int source_dof = node->dof_number(system.number(), 3, 0);
       if (source_dof <12345678){ //The pressures do not exist everywhere// This is a hack !!
       Real  value = 1*500*(1-exp(-pow((time/10),2.0)/0.25)) ;
       pressure_rows.push_back(source_dof);
      system.rhs->set(source_dof,value);
     }
    }  //end if
*/

    

#endif


  


/*
    if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 ) )
    {
       unsigned int source_dof = node->dof_number(system.number(), 3, 0);
       if (source_dof <12345678){ //The pressures do not exist everywhere// This is a hack !!
       Real  value = 1*(1-exp(-pow(time,2.0)/0.25)) ;
       pressure_rows.push_back(source_dof);
      system.rhs->set(source_dof,value);
     }
    }  //end if

*/


#endif

} //((elem->node(n) == side->node(ns))  )



#endif //end cube
} // end nodes on side loop

}//for loop
  
 }//end elem->neighbor(s) == NULL 
 }// end s=0; s<elem->n_sides(); s++

test(100);
 dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

#if NEUMANN_PRESSURE
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
#endif

} // end of element loop


#if DIRICHLET_VELOCITY
  test(3);
  system.matrix->close();
  test(4);
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



