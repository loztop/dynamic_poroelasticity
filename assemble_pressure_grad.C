#include "defines.h"
#include "assemble.h"
#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"
#include "poro_elastic_cc.h"


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_pressure_grad (EquationSystems& es,
                      const std::string& system_name)
{


    libmesh_assert (system_name == "pressure-grad-system");
    
    const MeshBase& mesh = es.get_mesh();
    
    const unsigned int dim = mesh.mesh_dimension();
    
    LinearImplicitSystem & system =
      es.get_system<LinearImplicitSystem> ("pressure-grad-system");
  
  
 TransientLinearImplicitSystem & pressure_system =
    es.get_system<TransientLinearImplicitSystem> ("pressure-system");
    
    
    
    const System & ref_sys = es.get_system("Reference-Configuration"); 


    const unsigned int u_var = system.variable_number ("p_grad_u");
    const unsigned int v_var = system.variable_number ("p_grad_v");
    const unsigned int w_var = system.variable_number ("p_grad_w");
    const unsigned int p_var = system.variable_number ("p_grad_aux");
    
    FEType fe_vel_type = system.variable_type(u_var);
    
    FEType fe_pres_type = system.variable_type(p_var);
  


    AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
      



    AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  
    fe_vel->attach_quadrature_rule (&qrule);
    fe_pres->attach_quadrature_rule (&qrule);
    
    const std::vector<Real>& JxW = fe_vel->get_JxW();
    
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  
      const std::vector<std::vector<RealTensor> >& d2phi = fe_vel->get_d2phi();



      const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();


    const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
    
    const DofMap & dof_map = system.get_dof_map();
  

AutoPtr<FEBase> fe_face (FEBase::build(3, fe_vel_type));
                
    AutoPtr<QBase> qface(fe_vel_type.default_quadrature_rule(3-1));
  
    fe_face->attach_quadrature_rule (qface.get());

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
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size(); 
      const unsigned int n_w_dofs = dof_indices_w.size(); 
      const unsigned int n_p_dofs = dof_indices_p.size();
      
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Reposition the submatrices...  The idea is this:
      //
      //         -           -          -  -
      //        | Kuu Kuv Kup |        | Fu |
      //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
      //        | Kpu Kpv Kpp |        | Fp |
      //         -           -          -  -
      //
      // The \p DenseSubMatrix.repostition () member takes the
      // (row_offset, column_offset, row_size, column_size).
      // 
      // Similarly, the \p DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
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


std::vector<unsigned int> undefo_index; 
  std::vector<Number> press_defo;
      pressure_system.get_dof_map().dof_indices(elem, undefo_index,0);
      pressure_system.current_local_solution->get(undefo_index, press_defo);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

    VectorValue<Real> grad_p_mat;  
    

for (unsigned int l = 0; l != n_u_dofs; l++){
     grad_p_mat.add_scaled(dphi[l][qp], press_defo[l] ); 
      }
    //  std::cout<<grad_p_mat<<std::endl;

 RealTensor tense_p_mat;  
    
    //      std::cout<<" d2phi[2][qp] "<<d2phi[2][qp]<<std::endl;

for (unsigned int l = 0; l != n_u_dofs; l++){
     tense_p_mat.add_scaled(d2phi[l][qp], press_defo[l] ); 
      }
 //     std::cout<<tense_p_mat<<std::endl;



Real   p_fluid = 0.;
for (unsigned int l=0; l<n_u_dofs; l++)
            {
         p_fluid += phi[l][qp]*pressure_system.current_local_solution->el(dof_indices_u[l]);
            }


Real   u_undef = 0.;
for (unsigned int l=0; l<n_u_dofs; l++)
     {
         u_undef += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
     }
    // std::cout<<u_undef<<std::endl;

          // Assemble the u-velocity row
          // uu coupling
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);

       for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);

             for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
        
      
          for (unsigned int j=0; j<n_u_dofs; j++)
             Fu(j) += JxW[qp]*phi[j][qp]*grad_p_mat(0);

          for (unsigned int j=0; j<n_v_dofs; j++){
             Fv(j) += JxW[qp]*phi[j][qp]*grad_p_mat(1);
          }

           for (unsigned int j=0; j<n_w_dofs; j++){
             Fw(j) += JxW[qp]*phi[j][qp]*grad_p_mat(2);

         //  Fw(j) += JxW[qp]*phi[j][qp]*(tense_p_mat(0,0));

         //   Fw(j) += -JxW[qp]*dphi[j][qp](0)*p_fluid;

             for (unsigned int i=0; i<n_w_dofs; i++){
           //    Fw(j) += -JxW[qp]*dphi[j][qp](2)*p_fluid;
           }

           }
        } // end of the quadrature point qp-loop

        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  


for (unsigned int s=0; s<elem->n_sides(); s++){
   if (elem->neighbor(s) == NULL)
   {    
   AutoPtr<Elem> side (elem->build_side(s));
          const Real penalty = 1;
/*
std::vector<unsigned int> undefo_index; 
      std::vector<Number> press_defo;
      pressure_system.get_dof_map().dof_indices(side_elem, undefo_index,0);
      pressure_system.current_local_solution->get(undefo_index, press_defo);
*/
      const std::vector<std::vector<Real> >&  phi_face =   fe_face->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
        const std::vector<Real>& JxW_face = fe_face->get_JxW();
        const std::vector<Point>& qface_point = fe_face->get_xyz();
        const std::vector<Point>& face_normals =fe_face->get_normals();
        fe_face->reinit(elem,s);  

  for (unsigned int qp=0; qp<qface->n_points(); qp++)
                  {

        //Node *noode = side->get_node(0);
        //unsigned int source_dof = noode->dof_number(pressure_system.number(), 0, 0);
        //Real p_value = pressure_system.current_local_solution->el(source_dof);
        //Real y_value = ref_sys.current_local_solution->el(noode->dof_number(1, 1, 0));
        //Real z_value = ref_sys.current_local_solution->el(noode->dof_number(1, 2, 0));


Real   p_fluid = 0.;
for (unsigned int l=0; l<phi_face.size(); l++)
            {
         p_fluid += phi_face[l][qp]*pressure_system.current_local_solution->el(dof_indices_u[l]);
            }

//std::cout<< "p_fluid "<<p_fluid<<  std::endl;
//std::cout<< "JxW_face[qp] "<<JxW_face[qp]<<  std::endl;
//std::cout<< "JxW[qp] "<<JxW[qp]<<  std::endl;
                    const Number  value = 0;
                                                         
                    for (unsigned int i=0; i<phi_face.size(); i++){

                     //Fu(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
                     //Fv(i) += - JxW_face[qp]*value*face_normals[qp](1)*phi_face[i][qp];
                    // Fw(i) += - JxW_face[qp]*value*face_normals[qp](0)*phi_face[i][qp];
           //Fw(i) += JxW_face[qp]*face_normals[qp](0)*phi_face[i][qp]*p_fluid;


                for (unsigned int l=0; l<phi_face.size(); l++){
            //Fw(i) += JxW[qp]*face_normals[qp](0)*phi_face[i][qp]*phi_face[l][qp]*press_defo[l];
              }


        }
                   for (unsigned int i=0; i<phi_face.size(); i++){
                     for (unsigned int j=0; j<phi_face.size(); j++){
  //Kuu(i,j) += value*face_normals[qp](0)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  //Kvv(i,j) += value*face_normals[qp](1)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  //Kww(i,j) += p_fluid*face_normals[qp](0)*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
         }
      }     
    } //end qp
}}

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      } // end of element loop
    
    return;
  }




