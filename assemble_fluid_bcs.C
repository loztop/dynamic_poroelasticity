#include "defines.h"
#include "assemble.h"
#include <cmath>
using namespace std;
#define PI 3.14159265
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_fluid_bcs (EquationSystems& es)
{
  //std::cout<<"Now assembling Boundary Conditions"<<std::endl;  

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
const System & ref_sys = es.get_system("Reference-Configuration"); 

TransientLinearImplicitSystem & fluid_system =
    es.get_system<TransientLinearImplicitSystem> ("fluid-system-vel");
const unsigned int fluid_system_number = fluid_system.number ();
  // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = fluid_system.variable_number ("fluid_U_vel");
    const unsigned int v_var = fluid_system.variable_number ("fluid_V_vel");
    const unsigned int w_var = fluid_system.variable_number ("fluid_W_vel");
    const unsigned int p_var = fluid_system.variable_number ("fluid_P");
    

FEType fe_vel_type = fluid_system.variable_type(u_var);
FEType press_vel_type = fluid_system.variable_type(p_var);

std::vector< int > rows;
std::vector< int > pressure_rows;

 
#if DYNAMIC
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");
    const Real time    = es.parameters.get<Real>("time");

#endif
 test(1);
 //Build face
 #if NEUMANN_PRESSURE
AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));             
AutoPtr<QBase> qface(press_vel_type.default_quadrature_rule(3-1));  
fe_face->attach_quadrature_rule (qface.get());

    // AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_vel_type));
    // QGauss qface (dim-1, fe_vel_type.default_quadrature_order());
    // const std::vector<Real>& JxW_face = fe_face->get_JxW();
    // const std::vector<std::vector<Real> >& psi_face = fe_face->get_phi();	
#endif

  const DofMap & dof_map = fluid_system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), 
    Kvu(Ke), Kvv(Ke), Kvw(Ke), 
    Kwu(Ke), Kwv(Ke), Kww(Ke); 
    
#if INCOMPRESSIBLE
  DenseSubMatrix<Number>  Kup(Ke),Kvp(Ke),Kwp(Ke), Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);
#endif
    
  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);
#if INCOMPRESSIBLE
  DenseSubVector<Number>    Fp(Fe);
#endif
test(2);
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  
#if INCOMPRESSIBLE
  std::vector<unsigned int> dof_indices_p;
#endif

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
test(66);
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

for (unsigned int ns=0; ns<side->n_nodes(); ns++)
    {
       for (unsigned int n=0; n<elem->n_nodes(); n++){
    Node *node = elem->get_node(n);
       Point p;
  for (unsigned int d = 0; d < 3; ++d) {
        unsigned int source_dof = node->dof_number(1, d, 0);
        Real value = ref_sys.current_local_solution->el(source_dof);
        p(d)=value;
    }
    if ((elem->node(n) == side->node(ns)) && (p(0)<0.01)  )
{
    for (unsigned int qp=0; qp<qface->n_points(); qp++)
    {       

   Real value=0;
    
   // value = progress*0.1  ;

     #if CHAP_SWELL

   Real factor=1;

      value = factor*1000*(1-exp(-pow(time,2.0)/0.25)) ;
     // value = 1000*progress ;

     #endif

      //value = 1;
      for (unsigned int i=0; i<phi_face.size(); i++){
        
        #if MASS_NEUMANN_PRESSURE //In Chapelle there is no neumann pressure Bc in the mass concservation eqn.
      //  Fp(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        #endif

        #if MOMENTUM_NEUMANN_PRESSURE
        Fu(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
        Fv(i) += - JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
        Fw(i) += - JxW_face[qp]*value*face_normals[qp](2)*phi_face[i][qp];
        #endif
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
      	Real value = ref_sys.current_local_solution->el(source_dof);
      	p(d)=value;
   	}
    
//Apply Dirichlet Bcs properly
#if CUBE 
  if ((elem->node(n) == side->node(ns))  )
  {
  #if DIRICHLET_VELOCITY

#if SEALED_CUBE
if ((elem->node(n) == side->node(ns)) && (p(0)>1.4999)  )
    {
      unsigned int source_dof = node->dof_number(fluid_system.number(), 0, 0);
       Real value = 0;
       rows.push_back(source_dof);
       fluid_system.rhs->set(source_dof,value);
    }  //end if
    
if ((elem->node(n) == side->node(ns)) && ((p(1)<0.0001)  || (p(1)>1.4999))  )
    {
       unsigned int source_dof = node->dof_number(fluid_system.number(), 1, 0);
       Real value = 0;
       rows.push_back(source_dof);
       fluid_system.rhs->set(source_dof,value);

    }  //end if

if ((elem->node(n) == side->node(ns)) && ((p(2)<0.0001)  || (p(2)>1.4999))  )
    {
      unsigned int source_dof = node->dof_number(fluid_system.number(), 2, 0);
       Real value = 0;
       rows.push_back(source_dof);
      fluid_system.rhs->set(source_dof,value);

    }  //end if
    
#endif
    
    #if DISK_FLOW
    if ((elem->node(n) == side->node(ns)) && (p(0)<0.1 )  )
    {
     for (unsigned int d = 1; d < 3; ++d) {
        unsigned int source_dof = node->dof_number(fluid_system.number(), d, 0);
       Real value = 0;
       rows.push_back(source_dof);
       fluid_system.rhs->set(source_dof,value);
    }  //end dimension loop

//std::cout<<"Disk flow"<<std::endl;
    
       unsigned int source_dof = node->dof_number(fluid_system.number(), 0, 0);
    Real value = progress*0.1* ( (0.75*0.75-(p(1)-0.75)*(p(1)-0.75)) + (0.75*0.75-(p(2)-0.75)*(p(2)-0.75)) ) ;
     //   Real value =0 ;
    //   rows.push_back(source_dof);
      // fluid_system.rhs->set(source_dof,value);
       
    }  //end if
   #endif 
   //DISK_FLOW


if ((elem->node(n) == side->node(ns)) && (p(0)<0.1)  )
    {
       unsigned int source_dof = node->dof_number(fluid_system.number(), 0, 0);

       //Real value = time*0.6;
       //rows.push_back(source_dof);
      // fluid_system.rhs->set(source_dof,value);

             // std::cout<< "value "<<value<<std::endl;


    }  //end if

/*
if ((elem->node(n) == side->node(ns)) && (p(1)<0.0001)  )
    {
       unsigned int source_dof = node->dof_number(fluid_system.number(), 1, 0);
       Real value = 0;
       rows.push_back(source_dof);
       fluid_system.rhs->set(source_dof,value);

    }  //end if
  */
    
#endif 


#if DIRICHLET_PRESSURE

#if CHAP_SWELL
    if ((elem->node(n) == side->node(ns)) && (p(0)>1.4999 ) )
    {
       unsigned int source_dof = node->dof_number(fluid_system.number(), 3, 0);
       if (source_dof <12345678){ //The pressures do not exist everywhere// This is a hack !!
       Real value = 0;
       pressure_rows.push_back(source_dof);
      fluid_system.rhs->set(source_dof,value);
     }
    }  //end if
#endif

/*
        if ((elem->node(n) == side->node(ns)) && (p(0)<0.01 ) )
    {
       unsigned int source_dof = node->dof_number(fluid_system.number(), 3, 0);
       if (source_dof <12345678){ //The pressures do not exist everywhere// This is a hack !!
       Real value = 1*time;
       pressure_rows.push_back(source_dof);
       fluid_system.rhs->set(source_dof,value);
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
      fluid_system.matrix->add_matrix (Ke, dof_indices);
      fluid_system.rhs->add_vector    (Fe, dof_indices);
#endif

} // end of element loop


#if DIRICHLET_VELOCITY
  test(3);
  fluid_system.matrix->close();
  test(4);
	fluid_system.matrix->zero_rows(rows, 1.0);
#endif

#if DIRICHLET_PRESSURE
   fluid_system.matrix->close();
   fluid_system.matrix->zero_rows(pressure_rows, 1.0);
#endif

test(5);
fluid_system.rhs->close();
test(6);

fluid_system.matrix->close();
//std::cout<<"AFTER BCS fluid_sefffffffffffffffsystem.rhs->l2_norm () "<<fluid_system.rhs->l2_norm ()<<std::endl;


//fluid_system.rhs->print(std::cout);
//fluid_system.matrix->print(std::cout);


  return;
}



