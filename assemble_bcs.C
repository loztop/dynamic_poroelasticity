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
  //std::cout<<"Now assembling Boundary Conditions"<<std::endl;  

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
 
#if DYNAMIC
    const Real dt    = es.parameters.get<Real>("dt");
    const Real progress    = es.parameters.get<Real>("progress");
    const Real time    = es.parameters.get<Real>("time");
#endif
 
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
 
  // The value of the linear shape function gradients at the quadrature points
  // const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  
  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
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
 
  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
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

//Now start actually applying the BCS        

for (unsigned int s=0; s<elem->n_sides(); s++){
   if (elem->neighbor(s) == NULL)
   {		
   AutoPtr<Elem> side (elem->build_side(s));





#if TRACTION_BC
  const std::vector<std::vector<Real> >&  phi_face =  fe_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Point>& face_normals = fe_face->get_normals();
  fe_face->reinit(elem,s);  

  const std::vector<Point>& qface_point_ref = fe_face_ref->get_xyz();
  fe_face_ref->reinit(elem,s);  


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

#if ANALNEO
 
   if ((elem->node(n) == side->node(ns)) && (p(0)>0.001)   )
{
for (unsigned int qp=0; qp<qface->n_points(); qp++)
    {


std::vector<unsigned int> undefo_index;

   Point p_qp_ref;
     for (unsigned int d = 0; d < 3; ++d) {
      std::vector<Number> u_undefo;
      ref_sys.get_dof_map().dof_indices(elem, undefo_index,d);
      ref_sys.current_local_solution->get(undefo_index, u_undefo);
      for (unsigned int l = 0; l != n_u_dofs; l++){
         p_qp_ref(d) += phi_face[l][qp]*u_undefo[l]; 

       }
    }






//std::cout<<" p_qp_ref " <<p_qp_ref <<std::endl;




//       if ((elem->node(n) == side->node(ns)) && (p_qp_ref(0)>0.001)  )


/*
Real a= 0.1;
Real b= 0.1;
Real c1= 0.1;

Real lam1 = 1+a*p_qp_ref(0);
Real lam2 = 1+b*p_qp_ref(1);

Real invlam1 = 1.0/lam1;
Real invlam2 = 1.0/lam2;

Real Z = p_qp_ref(2);


RealTensor S;
S(0,0)=lam1-invlam1;
S(0,2)=-a*Z*invlam1*invlam1*invlam2;
S(1,1)=lam2-invlam2;
S(1,2)=-b*Z*invlam1*invlam2*invlam2;
S(2,0)=-a*Z*invlam1*invlam1;
S(2,1)=-b*Z*invlam2*invlam2;
S(2,2)=invlam1*invlam2-lam1*lam2;

RealTensor F;
F(0,0)=lam1;
F(1,1)=lam2;
F(2,0)=-a*Z*invlam1*invlam1*invlam2;
F(2,1)=-b*Z*invlam1*invlam2*invlam2;
F(2,2)=invlam1*invlam2;
RealTensor Ft = F.transpose();
RealTensor sigma;
sigma=(1.0/(F.det()))*F*S*Ft;
//std::cout<<" trac before " <<traction <<std::endl;
DenseVector<Real> trac(3);
//tensor_mult_vector(trac, sigma, face_normals[qp]);
trac.scale(1.0/(F.det()));
  trac.scale(progress);
*/



DenseVector<Real> traction(3);
get_traction(traction,p_qp_ref, progress); 
//std::cout<<" get_traction " <<traction <<std::endl;


DenseVector<Real> traction_current(3);
//get_traction_current(traction_current,p_qp_ref,progress); 
get_traction_current(traction_current,p_qp_ref,progress); 


traction.scale(0);
traction_current.scale(0);


       for (unsigned int i=0; i<phi_face.size(); i++){
     
          Fu(i) += - JxW_face[qp]*traction(0)*phi_face[i][qp];
          Fv(i) += - JxW_face[qp]*traction(1)*phi_face[i][qp];
          Fw(i) += -JxW_face[qp]*traction(2)*phi_face[i][qp];




          Fu(i) += - JxW_face[qp]*traction_current(0)*face_normals[qp](0)*phi_face[i][qp];
          Fv(i) += - JxW_face[qp]*traction_current(1)*face_normals[qp](1)*phi_face[i][qp];
          Fw(i) += - JxW_face[qp]*traction_current(2)*face_normals[qp](2)*phi_face[i][qp];

/*
        Fu(i) += - JxW_face[qp]*traction_current(0)*phi_face[i][qp];
        Fv(i) += - JxW_face[qp]*traction_current(1)*phi_face[i][qp];
        Fw(i) += - JxW_face[qp]*traction_current(2)*phi_face[i][qp];
*/

        }

          for (unsigned int i=0; i<phi_face.size(); i++){
  
            for (unsigned int j=0; j<phi_face.size(); j++){
/*
  Kuu(i,j) += traction(0)*face_normals[qp](0)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  Kvv(i,j) += traction(1)*face_normals[qp](1)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  Kww(i,j) += traction(2)*face_normals[qp](2)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
*/
/*
    Kuu(i,j) += traction(0) *JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  Kvv(i,j) += traction(1)* JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  Kww(i,j) += traction(2) *JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
*/
          }
        }  
   

//iHat*anal_u +jHat*anal_v + kHat*anal_w


    } //end qp
  } //end if 
         #endif  
  }    
}
#endif //end NEUMANN_PRESSURE

/*
#if TRACTION_BC
 	const Real penalty = 1.e10;
        Node *noode = side->get_node(0);
    	unsigned int source_dof = noode->dof_number(1, 0, 0);
        Real x_value = ref_sys.current_local_solution->el(noode->dof_number(1, 0, 0));
        Real y_value = ref_sys.current_local_solution->el(noode->dof_number(1, 1, 0));
        Real z_value = ref_sys.current_local_solution->el(noode->dof_number(1, 2, 0));

        const std::vector<std::vector<Real> >&  phi_face =  fe_face->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
        const std::vector<Real>& JxW_face = fe_face->get_JxW();
        const std::vector<Point>& qface_point = fe_face->get_xyz();
        const std::vector<Point>& face_normals =
        fe_face->get_normals();
        fe_face->reinit(elem,s);  
        Real rsq=    x_value*x_value + y_value*y_value +z_value*z_value;         


	if( (y_value>0.99) && (x_value>3) )
	{
                for (unsigned int qp=0; qp<qface->n_points(); qp++)
                  {
                    const Number value = +10;
                                                         
                    for (unsigned int i=0; i<phi_face.size(); i++){

                     Fu(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
                     Fv(i) += - JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
                     Fw(i) += - JxW_face[qp]*value*face_normals[qp](2)*phi_face[i][qp];
		    }
                    for (unsigned int i=0; i<phi_face.size(); i++){
                      for (unsigned int j=0; j<phi_face.size(); j++){
	Kuu(i,j) += value*face_normals[qp](0)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
	Kvv(i,j) += value*face_normals[qp](1)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
	Kww(i,j) += value*face_normals[qp](2)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
		      }
		    }     
		} //end qp
	  } //end if
#endif //end Traction BCs
      */     





#if DIRICHLET_CLASSIC
//Apply Dirichlet Bcs properly
test(84);
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
	double rsquared = p(0)*p(0) + p(1)*p(1) + p(2)*p(2);
test(84);


#if CUBE 

/*
  if ((elem->node(n) == side->node(ns)) && (p(0)<0.01 )  )
  {
	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, 0, 0);
   	Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);
   	}  //end dimension loop
    }  //end if
*/
  

#if MOVING_DIRICHLET_BCS

  if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 )  )
  {
	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
   	Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);
   	}  //end dimension loop


    }  //end if

    if ((elem->node(n) == side->node(ns)) && (p(0)>1.49 )  )
    {
	for (unsigned int d = 0; d < 3; ++d) {

      	unsigned int source_dof = node->dof_number(1, 0, 0);
	
	Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);

        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);
   	}  //end dimension loop


#if DYNAMIC
 //   if (progress<0.1){

      	unsigned int source_dof = node->dof_number(1, 0, 0);
	Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof) - sin(4*progress*PI)*0.6; //sin(progress*PI)*0.1;
        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);



//}
#endif

     }  //end if
     
#endif
     


#if CHAP_SWELL

    if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 )  )
    {
  unsigned int source_dof = node->dof_number(1, 0, 0);
  Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
  rows.push_back(source_dof);
  newton_update.rhs->set(source_dof,value);

     }  //end if
 
    if ((elem->node(n) == side->node(ns)) && (p(1)<0.0001 )  )
    {
  unsigned int source_dof = node->dof_number(1, 1, 0);
  Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
  rows.push_back(source_dof);
  newton_update.rhs->set(source_dof,value);
     }  //end if


   if ((elem->node(n) == side->node(ns)) && (p(2)<0.0001 )  )
    {
  unsigned int source_dof = node->dof_number(1, 2, 0);
  Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
  rows.push_back(source_dof);
  newton_update.rhs->set(source_dof,value);
     }  //end if 
#endif


#if ANALNEO

// if ((elem->node(n) == side->node(ns)) && ( (p(0)<0.001 )  || (p(2)<0.001 ) || (p(2)>0.999 ) || (p(1)<0.001 ) || (p(1)>0.999 )  )   )


  if ((elem->node(n) == side->node(ns)) && ( p(0)<0.001)  ) 
  {

Real X=p(0);
Real Y=p(1);
Real Z=p(2);

  Real a2= progress*0.0;
  Real b2= progress*0.0;

  DenseVector<Real> new_position(3);
 new_position(0) = X + X*X*a2/2.0;
 new_position(1) = Y + Y*Y*b2/2.0;
 new_position(2) = Z/((1+X*a2)*(1+Y*b2));
//std::cout<< "prog 1 "<<new_position<<std::endl;
// new_position.scale(1);
//std::cout<< "progress "<<progress<<std::endl;

	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
   	Real value = last_non_linear_soln.current_local_solution->el(source_dof) -new_position(d);
        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);
   	}  //end dimension loop
    }  //end if

/*
if ((elem->node(n) == side->node(ns)) && ( p(0)>0.999)  ) 
  {

Real X=p(0);
Real Y=p(1);
Real Z=p(2);

  DenseVector<Real> new_position(3);
 new_position(0) = X + 0.5*progress;
 new_position(1) = Y ;
 new_position(2) = Z;


	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
   	Real value = last_non_linear_soln.current_local_solution->el(source_dof) -new_position(d);
        rows.push_back(source_dof);
	newton_update.rhs->set(source_dof,value);
   	}  //end dimension loop
    }  //end if

*/

#endif


#endif       
    } //end nodes in element lopp
} // end nodes on side loop
#endif


  
 }//end elem->neighbor(s) == NULL 
 }// end s=0; s<elem->n_sides(); s++
test(87);

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p SparseMatrix::add_matrix()
      // and \p NumericVector::add_vector() members do this for us.

#if TRACTION_BC || PENALTY
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



