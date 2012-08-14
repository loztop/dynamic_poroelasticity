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
//S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC) ;
S =  2.0*K1*pow(I_3,(-1.0/3.0))*Identity + 2.0*K2*I_1*pow(I_3,(-2.0/3.0))*Identity - 2.0*K2*pow(I_3,(-2.0/3.0))*C - K1*(2.0/3.0)*I_1*pow(I_3,(-1.0/3.0))*invC - K2*(4.0/3.0)*I_1*pow(I_3,(-2.0/3.0))*invC;
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
      ///this->f = A*exp(D*( pow(I_3,(-1.0/3.0))*I_1*(1+Q*m)-3.0));  Only ned this for the KCL model
      //this->gamma = D*(1.0+Q*m);

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



//Set parameters

K1=2000;
K2=33;
K=2.2 *100000;  
M=2.18*100000;
Kperm=0.0000001;


K1=1;
K2=1;
K=1;  
M=1;
Kperm=1;

K1=2000;
K2=33;
K=2.2 *100000;  
M=2.0 *100000;
Kperm=0.0001;


      this->c_update(C);

       if (this->calculate_linearized_stiffness) {
               this->calculate_tangent();
       }

   //    this->calculate_stress();

#if PORO
              this->calculate_stress_poro();
#endif

}



void PoroelasticConfig::calculate_tangent() {
       Real mu = E / (2 * (1 + nu));
       Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
       Real detF = F.det();

       C_mat.resize(6, 6);
        
#if CHAP

/*
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
*/

Real delta_am=-2*K1*(1.0/3.0);
DenseMatrix<Real> ICinv_matm;
ICinv_matm.resize(6, 6);
tensorOtensor_to_voigt(b,Identity,ICinv_matm);
ICinv_matm.scale(delta_am*pow(I_3,(-1.0/3.0)));
C_mat+=ICinv_matm;

Real delta_bm=2*K2;
DenseMatrix<Real> II_matm;
II_matm.resize(6, 6);
tensorOtensor_to_voigt(b,b,II_matm);
II_matm.scale(delta_bm*pow(I_3,(-2.0/3.0)));
C_mat+=II_matm;

ICinv_matm.resize(6, 6);
tensorOtensor_to_voigt(b,Identity,ICinv_matm);
ICinv_matm.scale(-delta_bm*(2.0/3.0)*I_1*pow(I_3,(-2.0/3.0)));
C_mat+=ICinv_matm;

Real delta_cm=-2*K2;
DenseMatrix<Real> CinvC_matm;
CinvC_matm.resize(6, 6);
tensorOtensor_to_voigt(b*C,Identity,CinvC_matm);
CinvC_matm.scale(delta_cm*(-2.0/3.0)*pow(I_3,(-2.0/3.0)));
C_mat+=CinvC_matm;

DenseMatrix<Real> I_matm;
I_matm.resize(6, 6);
z_ref_to_voigt(b,b,I_matm);
I_matm.scale(delta_cm*pow(I_3,(-2.0/3.0)));
C_mat+=I_matm;


Real delta_dm=-(2.0/3.0)*K1;
DenseMatrix<Real> invCI_matm;
invCI_matm.resize(6, 6);
tensorOtensor_to_voigt(Identity,b,invCI_matm);
invCI_matm.scale(delta_dm*pow(I_3,(-1.0/3.0)));
C_mat+=invCI_matm;

DenseMatrix<Real> invCinvC_matm;
invCinvC_matm.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_matm);
invCinvC_matm.scale(delta_dm*(-1.0/3.0)*I_1*pow(I_3,(-1.0/3.0)));
C_mat+=invCinvC_matm;

DenseMatrix<Real> Z_matm;
Z_matm.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_matm);
Z_matm.scale(delta_dm*(-1.0)*I_1*pow(I_3,(-1.0/3.0)));
C_mat+=Z_matm;

Real delta_em=-(4.0/3.0)*K2;
invCI_matm.resize(6, 6);
tensorOtensor_to_voigt(Identity,b,invCI_matm);
invCI_matm.scale(delta_em*pow(I_3,(-2.0/3.0)));
C_mat+=invCI_matm;

invCinvC_matm.resize(6, 6);
tensorOtensor_to_voigt(Identity,Identity,invCinvC_matm);
invCinvC_matm.scale(delta_em*(-2.0/3.0)*I_1*pow(I_3,(-2.0/3.0)));
C_mat+=invCinvC_matm;

Z_matm.resize(6, 6);
z_ref_to_voigt(Identity,Identity,Z_matm);
Z_matm.scale(delta_em*(-1.0)*I_1*pow(I_3,(-2.0/3.0)));
C_mat+=Z_matm;


    


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