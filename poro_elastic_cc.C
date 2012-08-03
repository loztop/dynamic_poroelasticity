#include "poro_elastic_cc.h"
#include "defines.h"
#include <math.h>

#if PORO


#if CHAP
Real PoroelasticConfig::calc_fchap(Real J) {
Real ans;
if (abs(J-1)<0.001){
  ans=1;
}else{
ans=2*(J-1-log(J))/(pow((J-1),2));  
}
return ans;
}

Real PoroelasticConfig::calc_fchapd(Real J) {
Real ans;  
if (abs(J-1)<0.001){
  ans=0;
}else{
ans=2*(  -2*(J-1-log(J))/(pow((J-1),3)) + (1-(1.0/J))/(pow((J-1),2)) );  
}
return ans;
}

Real PoroelasticConfig::calc_fchapdd(Real J) {
Real ans;  
if (abs(J-1)<0.001){
  ans=0;
}else{
ans=2*(  1/(J*J*pow((J-1),2)) - 2*(1-(1.0)/J)/(pow((J-1),3)) - 2*(  -3*(J-1-log(J))/(pow((J-1),4)) + (1-(1.0/J))/(pow((J-1),3)) )        );  
}
return ans;
}

Real PoroelasticConfig::calc_mchap() {
Real ans;  

ans= f_density*(J+(p_fluid-p_fluid_zero)/(M*fchap) -1);  

return ans;
}

#endif




void PoroelasticConfig::calculate_stress_poro() {

       double mu = E / (2.0 * (1.0 + nu));
       double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
       Real detF = F.det();
       RealTensor Ft = F.transpose();
       RealTensor C = Ft * F;
       RealTensor b = F * Ft;
       RealTensor identity;
       identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
       RealTensor invC = inv(C);




#if INCOMPRESSIBLE_CHEAT
S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC) + p_solid*detF*invC;
#endif

#if COMPRESSIBLE
  //This is the 2nd Piola-Kirchoff stress tensor, see 5.28 in Bonet. (not sure about log(J) term in Bonet)
       S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);
#endif

#if VERIFY_TEST_FINAL
S = 2*gamma*pow(I_3,(-1.0/3.0))*f*Identity+ 2*( (-1.0/3.0)*pow(gamma,2.0) * I_1*pow(I_3,(-2.0/3.0))*f*invC) +1*(p_solid*J*invC) ;
#endif 


#if CHAP
S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC) ;
S +=  1*K*J*invC + 1*(-1.0)*K*invC ;
S +=  -M*fchap*J*J*invC + M*fchap*J*invC  - 0.5*M*J*J*J*fchapd*invC + M*J*J*fchapd*invC -0.5*M*J*fchapd*invC  - (p_fluid-p_fluid_zero)*J*invC  + 0.5*(pow((p_fluid-p_fluid_zero),2.0)/M)*J*fchapd*pow(fchap,-2.0)*invC;



 //std::cout<<"p_fluid "<<p_fluid<<std::endl;

//       S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);

#endif 

    //   std::cout<< " S " << S <<std::endl;

      //convert to current configuration using a push forward operation
       tau = (F * S) * Ft;
       sigma = 1/detF * tau;
}

#if INCOMPRESSIBLE
void PoroelasticConfig::get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) {
       Real detF = F.det();
       p_residuum.resize(1);
     //  std::cout<<"poro m "<<m<<std::endl;
       p_residuum(0)=(1/detF) *psi[i][current_qp]*(detF-1-m);   ///Rp according to chaste fem *psi[i][current_qp] in main code       
}
#endif

void PoroelasticConfig::calculate_fluid_pressure() {
     this->p_fluid =  1*f*D* pow(I_3,(-1.0/3.0))*I_1*Q - 1*p_solid ;  
}




void PoroelasticConfig::c_update(RealTensor C) {
     
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
      this->mchap = calc_mchap();
      #endif


}

void PoroelasticConfig::init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m, Real p_fluid) {
       this->current_qp = qp;
        this->m = m;
  this->p_fluid=p_fluid;

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

      this->c_update(C);

       if (this->calculate_linearized_stiffness) {
               this->calculate_tangent();
       }

       this->calculate_stress();

#if PORO
              this->calculate_stress_poro();
#endif

}



void PoroelasticConfig::calculate_tangent() {
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

#endif

#if CHAP
/*
Real delta_a=1*2*k2;
DenseMatrix<Real> II_mat;
II_mat.resize(6, 6);
tensorOtensor_to_voigt(b,b,II_mat);
II_mat.scale(delta_a);
C_mat+=II_mat;


Real delta_b=1*(-2)*k2;
DenseMatrix<Real> I_mat;
I_mat.resize(6, 6);
z_ref_to_voigt(b,b,I_mat);
I_mat.scale(delta_b);
C_mat+=I_mat;

Real delta_c=1*J*K;
DenseMatrix<Real> invCinvC_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_c*0.5);
C_mat+=invCinvC_mat;

DenseMatrix<Real> Z_mat;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale(-1.0*delta_c);
C_mat+=Z_mat;


Real delta_d=1*(-1.0)*K;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_d);
C_mat+=Z_mat;
*/
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


Real delta_c=1*J*K;
DenseMatrix<Real> invCinvC_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_c*0.5);
C_mat+=invCinvC_mat;

DenseMatrix<Real> Z_mat;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale(-1.0*delta_c);
C_mat+=Z_mat;


Real delta_d=1*(-1.0)*K;
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_d);
C_mat+=Z_mat;








Real delta_e=1*(-M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_e*J*fchap*J);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_e*( (0.5*pow(J,1.0))*fchap*2*J  +   0.5*pow(J,1.0)*fchapd*J*J  ));
C_mat+=invCinvC_mat;




Real delta_f=1*(M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_f*J*fchap);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_f*( (0.5*pow(J,1.0))*fchap  +   0.5*pow(J,1.0)*fchapd*J  ));
C_mat+=invCinvC_mat;


Real delta_g=1*(-0.5*M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_g*J*fchapd*J*J);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_g*( (0.5*pow(J,1.0))*fchapd*3*J*J  +   0.5*pow(J,1.0)*fchapdd*J*J*J  ));
C_mat+=invCinvC_mat;

Real delta_h=1*(M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_h*J*fchapd*J);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_h*( (0.5*pow(J,1.0))*fchapd*2*J  +   0.5*pow(J,1.0)*fchapdd*J*J  ));
C_mat+=invCinvC_mat;


Real delta_i=1*(-0.5*M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_i*J*fchapd);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_i*( (0.5*pow(J,1.0))*fchapd  +   0.5*pow(J,1.0)*fchapdd*J  ));
C_mat+=invCinvC_mat;


Real delta_j=1*( -1.0*(p_fluid-p_fluid_zero));
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_j*J);
C_mat+=Z_mat;

invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_j*(0.5*pow(J,1.0)));
C_mat+=invCinvC_mat;

/*
Real delta_k=1*0.5*(M*pow(m/f_density,2.0));
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_k*J*fchapd);
C_mat+=Z_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_k*( (0.5*pow(J,1.0))*fchapd  +   0.5*pow(J,1.0)*fchapdd*J  ));
C_mat+=invCinvC_mat;
*/


Real delta_k=1*0.5*(pow((p_fluid-p_fluid_zero),2.0)/M);
Z_mat.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_mat);
Z_mat.scale((-1.0)*delta_k*J*fchapd*pow(fchap,-2.0));
C_mat+=Z_mat;
invCinvC_mat.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
invCinvC_mat.scale(delta_k*( (0.5*pow(J,1.0))*fchapd*pow(fchap,-2.0)  +   0.5*pow(J,1.0)*fchapdd*J*pow(fchap,-2.0)  -1*pow(fchap,-3.0)*pow(fchapd,2.0)*pow(J,2.0) ));
C_mat+=invCinvC_mat;



#endif
     
}

#endif