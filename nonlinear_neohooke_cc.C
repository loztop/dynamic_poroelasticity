#include "nonlinear_neohooke_cc.h"
#include "defines.h"
#include <math.h>
#include "poro_elastic_cc.h"

//#define INCOMPRESSIBLE 1

/**
 * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
 * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
 */



/*
 #if !PORO
template <typename T> TypeTensor<T> inv(const TypeTensor<T> &A ) {
  double Sub11, Sub12, Sub13;
  Sub11 = A._coords[4]*A._coords[8] - A._coords[5]*A._coords[7];
  Sub12 = A._coords[3]*A._coords[8] - A._coords[6]*A._coords[5];
  Sub13 = A._coords[3]*A._coords[7] - A._coords[6]*A._coords[4];
  double detA = A._coords[0]*Sub11 - A._coords[1]*Sub12 + A._coords[2]*Sub13;
  libmesh_assert( std::fabs(detA)>1.e-15 );

  TypeTensor<T> Ainv(A);

  Ainv._coords[0] =  Sub11/detA;
  Ainv._coords[1] = (-A._coords[1]*A._coords[8]+A._coords[2]*A._coords[7])/detA;
  Ainv._coords[2] = ( A._coords[1]*A._coords[5]-A._coords[2]*A._coords[4])/detA;
  Ainv._coords[3] = -Sub12/detA;
  Ainv._coords[4] = ( A._coords[0]*A._coords[8]-A._coords[2]*A._coords[6])/detA;
  Ainv._coords[5] = (-A._coords[0]*A._coords[5]+A._coords[2]*A._coords[3])/detA;
  Ainv._coords[6] =  Sub13/detA;
  Ainv._coords[7] = (-A._coords[0]*A._coords[7]+A._coords[1]*A._coords[6])/detA;
  Ainv._coords[8] = ( A._coords[0]*A._coords[4]-A._coords[1]*A._coords[3])/detA;

  return Ainv;
}
#endif
*/

void NonlinearNeoHookeCurrentConfig::init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp) {
       this->current_qp = qp;
       
      #if INCOMPRESSIBLE 
      //this->p_current = p_current;
      this->p_solid = p_current;

      #endif
      
       
       F.zero();
       S.zero();

       {
         RealTensor invF;
         invF.zero();
         for (unsigned int i = 0; i < 3; ++i)
           for (unsigned int j = 0; j < 3; ++j) {
#if MOVING_MESH
             invF(i, j) += grad_u(i)(j);
#endif
#if FIXED_MESH
	              F(i, j) += grad_u(i)(j);
#endif
	//std::cout << "invF " << invF<< std::endl;
           }
#if MOVING_MESH
    F.add(inv(invF));
#endif
       }

       if (F.det() < -TOLERANCE) {
               std::cout << "detF < 0" << std::endl;
               libmesh_error();
       }


    this->Ft = F.transpose();
      this->C = Ft*F;

     // this->c_update(C);

       if (this->calculate_linearized_stiffness) {
               this->calculate_tangent();
       }

       this->calculate_stress();
}

void NonlinearNeoHookeCurrentConfig::init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m) {
       this->current_qp = qp;
        this->m = m;

      #if INCOMPRESSIBLE 
      //this->p_current = p_current;
      this->p_solid = p_current;

      #endif
      
       
       F.zero();
       S.zero();

       {
         RealTensor invF;
         invF.zero();
         for (unsigned int i = 0; i < 3; ++i)
           for (unsigned int j = 0; j < 3; ++j) {
#if MOVING_MESH
             invF(i, j) += grad_u(i)(j);
#endif
#if FIXED_MESH
                F(i, j) += grad_u(i)(j);
#endif
  //std::cout << "invF " << invF<< std::endl;
           }
#if MOVING_MESH
    F.add(inv(invF));
#endif
       }

       if (F.det() < -TOLERANCE) {
               std::cout << "detF < 0" << std::endl;
               libmesh_error();
       }


    this->Ft = F.transpose();
      this->C = Ft*F;

     // this->c_update(C);

       if (this->calculate_linearized_stiffness) {
               this->calculate_tangent();
       }

       this->calculate_stress();


}

void NonlinearNeoHookeCurrentConfig::calculate_tangent() {
       Real mu = E / (2 * (1 + nu));
       Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
       Real detF = F.det();

       C_mat.resize(6, 6);
   
#if COMPRESSIBLE    
       for (unsigned int i = 0; i < 3; ++i) {
               for (unsigned int j = 0; j < 3; ++j) {
                       if (i == j) {
                               C_mat(i, j) = 2 * mu + lambda;
                               C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
                       } else {
                               C_mat(i, j) = lambda * detF * detF;
                       }
               }
       }
#endif    

#if INCOMPRESSIBLE_CHEAT    

       for (unsigned int i = 0; i < 3; ++i) {
               for (unsigned int j = 0; j < 3; ++j) {
                       if (i == j) {
                               C_mat(i, j) = 2 * mu + lambda;
                               C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
                       } else {
                               C_mat(i, j) = lambda * detF * detF;
                       }
               }
       }




Real fac=0.5*p_current*detF;
Real delta_a=fac*(1.0/2.0);
DenseMatrix<Real> invCinvC_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_a);
C_mat+=invCinvC_mat;


Real delta_b=fac*(-1.0);
DenseMatrix<Real> Z_mat;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale(delta_b);
C_mat+=Z_mat;

 #endif    
      

#if VERIFY_TEST_FINAL
//THIRD TERM
///////////////////////////

       C_mat.resize(6, 6);

Real fac=1*p_solid*J;

Real delta_a=fac*(1.0/2.0);
DenseMatrix<Real> invCinvC_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_a);
C_mat+=invCinvC_mat;


Real delta_b=fac*(-1.0);
DenseMatrix<Real> Z_mat;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale(delta_b);
C_mat+=Z_mat;
       
///FIRST TERM - Checked- works
fac=2*gamma;
delta_a=fac*(-1.0/3.0)*f*pow(I_3,(-1.0/3.0));
 delta_b=fac*pow(I_3,(-1.0/3.0))*gamma*f;
Real delta_c=fac*(-1.0/3.0)*I_1*gamma*f*pow(I_3,(-1.0/3.0));

DenseMatrix<Real> IinvC_mat;
IinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(b,Identity,IinvC_mat);
IinvC_mat.scale(delta_a+delta_c);

DenseMatrix<Real> II_mat;
II_mat.resize(6, 6);
tensorOtensor_to_voigt(b,b,II_mat);
II_mat.scale(delta_b);

C_mat+=IinvC_mat;
C_mat+=II_mat;
////////////////////////////////SECOND TERM


fac=2*(-1.0/3.0)*gamma*gamma;

delta_a=fac*(-1.0)*I_1*pow(I_3,(-2.0/3.0))*f;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale(delta_a);
C_mat+=Z_mat;

delta_b=fac*I_1*pow(I_3,(-2.0/3.0))*pow(I_3,(-1.0/3.0))*gamma*f;
DenseMatrix<Real> invCI_mat;
invCI_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,b,invCI_mat);
invCI_mat.scale(delta_b);
C_mat+=invCI_mat;

delta_c=fac*I_1*pow(I_3,(-2.0/3.0))*(-1.0/3.0)*pow(I_3,(-1.0/3.0))*I_1*gamma*f;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_c);
C_mat+=invCinvC_mat;

Real delta_d=fac*(-2.0/3.0)*pow(I_3,(-2.0/3.0))*I_1*f;
DenseMatrix<Real> invCinvC_mat2;
invCinvC_mat2.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat2);
invCinvC_mat2.scale(delta_d);
C_mat+=invCinvC_mat2;

Real delta_e=fac*pow(I_3,(-2.0/3.0))*f;
DenseMatrix<Real> invCI_mat2;
invCI_mat2.resize(6, 6);
tensorOtensor_to_voigt(Identity,b,invCI_mat2);
invCI_mat2.scale(delta_e);
C_mat+=invCI_mat2;



//To take into account that C4=2 dS/dC
//C_mat.scale(2.0);


#endif

     
}


void NonlinearNeoHookeCurrentConfig::calculate_stress() {

       double mu = E / (2.0 * (1.0 + nu));
       double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
       Real detF = F.det();
       RealTensor Ft = F.transpose();
       RealTensor C = Ft * F;
       RealTensor b = F * Ft;
       RealTensor identity;
       identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
       RealTensor invC = inv(C);

#if COMPRESSIBLE
	//This is the 2nd Piola-Kirchoff stress tensor, see 5.28 in Bonet. (not sure about log(J) term in Bonet)
       S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);
#endif



#if INCOMPRESSIBLE_CHEAT    
//std::cout<< " p_solid " << p_solid <<std::endl;

S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC) + p_solid*detF*invC;

#endif

#if INCOMPRESSIBLE_STRESS  
       Real I_3 = C.det();
       Real I_1 = C(0,0)+C(1,1)+C(2,2);

//This law is broken !! or maybe the corresponding linearisation is wrong see nonlinear_neohooke_ccBACKUP_NEWVERSIOn for details - need to fix this LATER !
      S = mu*pow(I_3,-1/3)*(identity - (1/3)*I_1*C) + p_current*detF*invC;
#endif	      

      //convert to current configuration using a push forward operation
       tau = (F * S) * Ft;
       sigma = 1/detF * tau;
}

void NonlinearNeoHookeCurrentConfig::get_residual(DenseVector<Real> & residuum, unsigned int & i) {
       B_L.resize(3, 6);
       DenseVector<Real> sigma_voigt(6);
       this->build_b_0_mat(i, B_L);
       tensor_to_voigt(sigma, sigma_voigt);
       B_L.vector_mult(residuum, sigma_voigt);  //should this be B_L^{T} * sigma_voigt ?

      // std::cout<< " sigma " << sigma <<std::endl;

       }

#if INCOMPRESSIBLE
void NonlinearNeoHookeCurrentConfig::get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) {
       Real detF = F.det();
       p_residuum.resize(1);
         //     std::cout<<"neo  "<<std::cout;

       p_residuum(0)=(1/detF) *psi[i][current_qp]*(detF-1);   ///Rp according to chaste fem *psi[i][current_qp] in main code       
}
#endif   


void NonlinearNeoHookeCurrentConfig::tensor_to_voigt(const RealTensor &tensor, DenseVector<Real> &vec) {
  vec(0) = tensor(0, 0);
  vec(1) = tensor(1, 1);
  vec(2) = tensor(2, 2);
  vec(3) = tensor(0, 1);
  vec(4) = tensor(1, 2);
  vec(5) = tensor(0, 2);

}

void NonlinearNeoHookeCurrentConfig::tensorOtensor_to_voigt(const RealTensor &tensorA, const RealTensor &tensorB, DenseMatrix<Real> &mat){
//Top left
  /*
 for (unsigned int i = 0; i < 3; ++i) {
               for (unsigned int j = 0; j < 3; ++j) {
                               mat(i, j) = 2*tensorA(i,i)*tensorB(j,j);               
               }
       }
*/

//Top three rows
for (unsigned int i = 0; i < 3; ++i) {
mat(i,0) = tensorA(i,i)*tensorB(0,0)+tensorA(i,i)*tensorB(0,0);          
mat(i,1) = tensorA(i,i)*tensorB(1,1)+tensorA(i,i)*tensorB(1,1);               
mat(i,2) = tensorA(i,i)*tensorB(2,2)+tensorA(i,i)*tensorB(2,2);
mat(i,3) = tensorA(i,i)*tensorB(0,1)+tensorA(i,i)*tensorB(1,0);          
mat(i,4) = tensorA(i,i)*tensorB(0,2)+tensorA(i,i)*tensorB(2,0);               
mat(i,5) = tensorA(i,i)*tensorB(1,2)+tensorA(i,i)*tensorB(2,1);                 
}



//fourth row
int i=3; int a=0; int b=1;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);                 

//fith row
 i=4;  a=0;  b=2;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);   

//sixth row
 i=5;  a=1;  b=2;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);   

}


void NonlinearNeoHookeCurrentConfig::z_ref_to_voigt(const RealTensor &tensorA, const RealTensor &tensorB, DenseMatrix<Real> &mat){
//Top left
  /*
 for (unsigned int i = 0; i < 3; ++i) {
               for (unsigned int j = 0; j < 3; ++j) {
                               mat(i, j) = 2*tensorA(i,i)*tensorB(j,j);               
               }
       }
*/


//Top three rows
for (unsigned int i = 0; i < 3; ++i) {
mat(i,0) = tensorA(i,0)*tensorB(i,0)*2;          
mat(i,1) = tensorA(i,1)*tensorB(i,1)*2;             
mat(i,2) = tensorA(i,2)*tensorB(i,2)*2;
mat(i,3) = tensorA(i,0)*tensorB(i,1)*2;       
mat(i,4) = tensorA(i,0)*tensorB(i,2)*2;         
mat(i,5) = tensorA(i,1)*tensorB(i,2)*2;                 
}



//fourth row
int i=3; int a=0; int b=1;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);                 

//fith row
 i=4;  a=0;  b=2;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);    

//sixth row
 i=5;  a=1;  b=2;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);     

}



void NonlinearNeoHookeCurrentConfig::get_linearized_stiffness(DenseMatrix<Real> & stiffness, unsigned int & i, unsigned int & j) {
       stiffness.resize(3, 3);

       double G_IK = (sigma * dphi[i][current_qp]) * dphi[j][current_qp];
       stiffness(0, 0) += G_IK;
       stiffness(1, 1) += G_IK;
       stiffness(2, 2) += G_IK;

       B_L.resize(3, 6);
       this->build_b_0_mat(i, B_L);
       B_K.resize(3, 6);
       this->build_b_0_mat(j, B_K);

       B_L.right_multiply(C_mat);
       B_L.right_multiply_transpose(B_K);
       B_L *= 1/F.det();

       stiffness += B_L;
}

#if INCOMPRESSIBLE
void NonlinearNeoHookeCurrentConfig::get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
  // Find and write down the mathematics for this section.
// eulerian tangent matrix at 6.15 - Bonet
RealTensor Ft = F.transpose();
Real detF = F.det();
RealTensor C = Ft * F;
RealTensor invC = inv(C);
RealTensor A = invC*F;
p_stiffness.resize(3);
RealTensor invF = inv(F);

//Build K_B (differentiate u eqn in p direction).
RealTensor DPHIi;

for (unsigned int z = 0; z < 3; ++z) {
	  DPHIi(z,0)=dphi[i][current_qp](z); 
}
//std::cout<< "DPHIi " << DPHIi <<std::endl;
RealTensor invFt = invF.transpose();

RealTensor uvw_p_stiff=psi[j][current_qp]*invFt*DPHIi;
//std::cout<< "DPHIi " << DPHIi <<std::endl;
//std::cout<< "invFt" << invFt <<std::endl;
//std::cout<< "uvw_p_stiff" << uvw_p_stiff <<std::endl;
p_stiffness(0)=uvw_p_stiff(0,0); 
p_stiffness(1)=uvw_p_stiff(1,0); 
p_stiffness(2)=uvw_p_stiff(2,0); 



////////
p_stiffness.resize(3);
p_stiffness(0)=psi[j][current_qp]*dphi[i][current_qp](0);
p_stiffness(1)=psi[j][current_qp]*dphi[i][current_qp](1);
p_stiffness(2)=psi[j][current_qp]*dphi[i][current_qp](2);
/////////


}

void NonlinearNeoHookeCurrentConfig::get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
//Build K_C (differentiate p eqn in u direction).
RealTensor Ft = F.transpose();
Real detF = F.det();
RealTensor C = Ft * F;
RealTensor invC = inv(C);
RealTensor A = invC*F;
p_stiffness.resize(3);
RealTensor invF = inv(F);
RealTensor invFt = invF.transpose();

RealTensor DPHIj;
for (unsigned int z = 0; z < 3; ++z) {
	  DPHIj(z,0)=dphi[j][current_qp](z); 
}
RealTensor p_uvw_stiff= psi[i][current_qp]*invFt*DPHIj;
p_stiffness(0)=p_uvw_stiff(0,0);
p_stiffness(1)=p_uvw_stiff(1,0); 
p_stiffness(2)=p_uvw_stiff(2,0); 

p_stiffness.resize(3);
p_stiffness(0)=psi[i][current_qp]*dphi[j][current_qp](0);
p_stiffness(1)=psi[i][current_qp]*dphi[j][current_qp](1);
p_stiffness(2)=psi[i][current_qp]*dphi[j][current_qp](2);

}




#endif

void NonlinearNeoHookeCurrentConfig::build_b_0_mat(int i, DenseMatrix<Real>& b_0_mat) {
       for (unsigned int ii = 0; ii < 3; ++ii) {
               b_0_mat(ii, ii) = dphi[i][current_qp](ii);
       }
       //See wriggers 4.94
       b_0_mat(0, 3) = dphi[i][current_qp](1);
       b_0_mat(1, 3) = dphi[i][current_qp](0);
       b_0_mat(1, 4) = dphi[i][current_qp](2);
       b_0_mat(2, 4) = dphi[i][current_qp](1);
       b_0_mat(0, 5) = dphi[i][current_qp](2);
       b_0_mat(2, 5) = dphi[i][current_qp](0);
}

/*
void NonlinearNeoHookeCurrentConfig::c_update(RealTensor C) {
     
    // Real m=0;

      this-> C = C;
      this->invC = inv(C);
      this->b = F*Ft;

      this->Identity(0,0)=1.0;
      this->Identity(1,1)=1.0;
      this->Identity(2,2)=1.0;
      this->I_3 = C.det();
      this->I_1 = C(0,0)+C(1,1)+C(2,2);

      RealTensor Csqd = C*C;
      this->I_2 = 0.5*(pow(I_1,2) - (Csqd(0,0)+Csqd(1,1)+Csqd(2,2))) ;
      this->J=pow(I_3,(1.0/2.0));
      this->f = A*exp(D*( pow(I_3,(-1.0/3.0))*I_1*(1+Q*m)-3.0));
      this->gamma = D*(1.0+Q*m);

       #if CHAP
      this->fchap = calc_fchap(J);
      this->fchapd = calc_fchapd(J);
      this->fchapdd = calc_fchapdd(J);
      #endif


}
*/