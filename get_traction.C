   // C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "transient_system.h"
#include "perf_log.h"
#include "boundary_info.h"
#include "utility.h"
#include "elem.h"
#include "mesh_data.h"
#include "gmsh_io.h"

// Some (older) compilers do not offer full stream 
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"


#include "assemble.h"
#include "defines.h"

//#include "solid_system.h"
#include <petsc_linear_solver.h>
using namespace std;



void get_traction(DenseVector<Real> & traction, Point rX, const Real progress) {


 Real a= 0.1*progress;
 Real b= 0.1*progress;
 Real c1= 0.1;

Real lam1 = 1+a*rX(0);
Real lam2 = 1+b*rX(1);

Real invlam1 = 1.0/lam1;
Real invlam2 = 1.0/lam2;

Real Z = rX(2);


if ( fabs(rX(0)-1)<1e-6 )
        {
            traction(0) =  lam1 - invlam1;
            traction(1) =  0.0;
            traction(2) =  -a*Z*invlam1*invlam1*invlam2;
        }
        else if ( fabs(rX(1)-0)<1e-6 )
        {
            traction(0) =  0.0;
            traction(1) =  0.0;
            traction(2) =  b*Z*invlam1;
        }
        else if ( fabs(rX(1)-1)<1e-6 )
        {
            traction(0) =  0.0;
            traction(1) =  lam2 - invlam2;
            traction(2) =  -b*Z*invlam2*invlam2*invlam1;
        }
       else if ( fabs(rX(2)-0)<1e-6 )
        {
            traction(0) =  0.0;
            traction(1) =  0.0;
            traction(2) =  lam1*lam2 - invlam1*invlam2;
        }
        else if ( fabs(rX(2)-1)<1e-6 )
        {
            traction(0) =  -a*invlam1*invlam1;
            traction(1) =  -b*invlam2*invlam2;
            traction(2) =  invlam1*invlam2 - lam1*lam2;
        }
        else
        {
        }
traction.scale(2*c1);

}

void get_traction_current(DenseVector<Real> & traction, Point rX, const Real progress) {

//rX(0)=40;
//rX(1)=10.34;
//rX(2)=1;


Real a= 0.1*progress;
Real b= 0.1*progress;
Real c1= 0.1;
Real lam1 = 1+a*rX(0);
Real lam2 = 1+b*rX(1);
Real invlam1 = 1.0/lam1;
Real invlam2 = 1.0/lam2;
Real Z = rX(2);
Real X = rX(0);
Real Y = rX(1);

RealTensor P;
P(0,0)=lam1-invlam1;
P(0,2)=-a*Z*invlam1*invlam1*invlam2;
P(1,1)=lam2-invlam2;
P(1,2)=-b*Z*invlam1*invlam2*invlam2;
P(2,0)=-a*Z*invlam1*invlam1;
P(2,1)=-b*Z*invlam2*invlam2;
P(2,2)=invlam1*invlam2-lam1*lam2;
P= P*2*c1;

RealTensor F;
F(0,0)=lam1;
F(1,1)=lam2;
F(2,0)=-a*Z*invlam1*invlam1*invlam2;
F(2,1)=-b*Z*invlam1*invlam2*invlam2;
F(2,2)=invlam1*invlam2;

RealTensor Pt = P.transpose();
RealTensor Ft = F.transpose();
RealTensor invF = inv(F);
RealTensor invFt = invF.transpose();

RealTensor S;
S=P*invFt;

RealTensor sigma;
sigma=(1.0/(F.det()))*F*S*Ft;


//std::cout<<" sigma " <<sigma <<std::endl;
//std::cout<<" (1.0/(F.det()))*F*P " <<(1.0/(F.det()))*F*P <<std::endl;



Point normal_ref;
Point normal_current;


if ( fabs(rX(0)-19999)<1e-6 )
        {
            normal_ref(0) =  1.0;
            normal_ref(1) =  0.0;
            normal_ref(2) =  0.0;

             normal_current(0) =  1.0;
            normal_current(1) =  0.0;
            normal_current(2) =  0.0;

        }
        else if ( fabs(rX(1)-0)<1e-6 )
        {
            normal_ref(0) =  0.0;
            normal_ref(1) =  -1.0;
            normal_ref(2) =  0.0;

            normal_current(0) =  0.0;
            normal_current(1) =  -1.0;
            normal_current(2) =  0.0;
        }
        else if ( fabs(rX(1)-1)<1e-6 )
        {
             normal_ref(0) =  0.0;
            normal_ref(1) =  1.0;
            normal_ref(2) =  0.0;

             normal_current(0) =  0.0;
            normal_current(1) =  1.0;
            normal_current(2) =  0.0;
        } 
      else if ( fabs(rX(2)-0)<1e-6 )
        {
            normal_ref(0) =  0.0;
            normal_ref(1) =  0.0;
            normal_ref(2) =  -1.0;

            normal_current(0) =  0.0;
            normal_current(1) =  0.0;
            normal_current(2) =  -1.0;
        } 
          else if ( fabs(rX(2)-1)<1e-6 )
        {
           normal_ref(0) =  0.0;
            normal_ref(1) =  0.0;
            normal_ref(2) =  1.0;

    Real fx=pow(1+a*X,-2.0)*pow(1+b*Y,-1.0)*-a;
    Real fy=pow(1+a*X,-1.0)*pow(1+b*Y,-2.0)*-b;

             normal_current(0) =  -fx;
            normal_current(1) =  -fy;
            normal_current(2) =  1.0;

            //normalize normal
          //  normal_current=(1.0)*normal_current*((1.0)/(  pow(fx,2.0) +  pow(fy,2.0) +pow(-1,2.0)     ));

         //    normal_current(0) = 0;
         //   normal_current(1) =  0;
         //   normal_current(2) =  1.0;
        }  
        else
        {
        }





/*
DenseVector<Real> normal_current_vec(3);
Point test_vec;
test_vec(0) =  traction_ref(0);
test_vec(1) =  traction_ref(1);
test_vec(2) =  traction_ref(2) ;

tensor_mult_vector(normal_current_vec,  F.det()*invFt, test_vec);
std::cout<<" normal_current_vec " <<normal_current_vec <<std::endl;


RealTensor V;
V(0,0)=1;
V(1,1)=1;
V(2,0)=0;
V(2,1)=0;
V(2,2)=2;
RealTensor invV=inv(V);
RealTensor invVt=invV.transpose();
std::cout<<" F.det() " <<F.det() <<std::endl;


tensor_mult_vector(normal_current_vec,  F.det()*invFt, test_vec);





  Real fx=pow(1+a*X,-2.0)*pow(1+b*Y,-1.0)*-a;
    Real fy=pow(1+a*X,-1.0)*pow(1+b*Y,-2.0)*-b;

             normal_current(0) =  -fx;
            normal_current(1) =  -fy;
            normal_current(2) =  1.0;

            std::cout<<" normal_current " << normal_current <<std::endl;
*/

/*
//tensor_mult_vector(normal_current_vec,  F.det()*invFt, normal_ref);
//tensor_mult_vector(normal_current_vec,  Ft, normal_ref);

normal_current(0)=normal_current_vec(0);
normal_current(1)=normal_current_vec(1);
normal_current(2)=normal_current_vec(2);
*/

//std::cout<<" normal_current " <<normal_current <<std::endl;
//std::cout<<" normal_ref " << normal_ref <<std::endl;



DenseVector<Real> traction_ref(3);
tensor_mult_vector(traction_ref, Pt, normal_ref);
//std::cout<<" traction_ref " <<traction_ref <<std::endl;

tensor_mult_vector(traction, sigma, normal_ref);
//std::cout<<" traction_current " <<traction <<std::endl;


traction.scale(1.0/(F.det()));

}


void get_bodyforce(DenseVector<Real> & body_force, Point rX, const Real progress) {

static Real a= 0.1*progress;
static Real b= 0.1*progress;
static Real c1= 0.1;

        Real lam1 = 1+a*rX(0);
        Real lam2 = 1+b*rX(1);
        Real invlam1 = 1.0/lam1;
        Real invlam2 = 1.0/lam2;

        body_force(0) = a;
        body_force(1) = b;
        body_force(2) = 2*rX(2)*invlam1*invlam2*( a*a*invlam1*invlam1 + b*b*invlam2*invlam2 );

	body_force.scale(-2*c1);

    }


  void tensor_mult_vector(DenseVector<Real> & ans, RealTensor tens, Point normal){

//std::cout<< " normal "<< normal <<std::endl;
//std::cout<< " tens "<< tens <<std::endl;


ans(0)=tens(0,0)*normal(0)+tens(0,1)*normal(1)+tens(0,2)*normal(2);


ans(1)=tens(1,0)*normal(0)+tens(1,1)*normal(1)+tens(1,2)*normal(2);


ans(2)=tens(2,0)*normal(0)+tens(2,1)*normal(1)+tens(2,2)*normal(2);

  }










  void get_traction_test(DenseVector<Real> & traction, Point rX, const Real progress) {

Real a= 0.01*progress;
Real b= 0.1*progress;
Real c1= 0.1;

Real lam1 = 1+a*rX(0);
Real lam2 = 1+b*rX(1);

Real invlam1 = 1.0/lam1;
Real invlam2 = 1.0/lam2;

Real Z = rX(2);
Real X = rX(0);
Real Y = rX(1);




if ( fabs(rX(0)-1)<1e-6 )
        {
            traction(0) =  0;
            traction(1) =  a;
            traction(2) =  0.0;



        }
        else if ( fabs(rX(1)-0)<1e-6 )
        {
   

        }
        else if ( fabs(rX(1)-1)<1e-6 )
        {


        } 
      else if ( fabs(rX(2)-0)<1e-6 )
        {


        } 
          else if ( fabs(rX(2)-1)<1e-6 )
        {

        }  
        else
        {
        }




}
