#ifndef poro_elastic_CC_H_
#define poro_elastic_CC_H_

#include <dense_vector.h>
#include <dense_matrix.h>
#include <vector_value.h>
#include <tensor_value.h>
#include <getpot.h>
#include <math.h>


#include "defines.h"
#include "general_material_cc.h"



//#if PORO
/**
 * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
 * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
 */


class PoroelasticConfig : public GeneralMaterialConfig {

public:
#if INCOMPRESSIBLE || (CHAP) 
PoroelasticConfig(const std::vector<std::vector<RealGradient> >& dphi, const std::vector<std::vector<Real> >& psi) :  GeneralMaterialConfig(dphi, psi){

//Stuff I don't really need.
E = 1;
nu = 0.3;
A=1.0; 
D=1.0; 
Q=1.0; 
///

f_phase_ref=0.1;
f_density=1;

#if CHAP
K1=2000;
K2=33;
//nu = 68
K=2.2 *100000;  
//K=1;
//M=10000;
M=2.18*100;
//bzero=1 always, this is hardcoded into the code. 
//kzero =0.01 
rho_solid_zero=1000; //currently only looking at quasi static solid
rho_fluid_zero=1000;
p_fluid_zero=0;      //Doesn't even feature in code at the moment
Kperm=0.0000001;
mchap=0;
#endif
}

PoroelasticConfig(const std::vector<std::vector<RealGradient> >& dphi, const std::vector<std::vector<Real> >& psi, Real p_fluid) : p_fluid(p_fluid),  GeneralMaterialConfig(dphi, psi){

//This constructor is not used
//Stuff I don't really need.
E = 1;
nu = 0.3;
A=1.0; 
D=1.0; 
Q=1.0; 
///

f_phase_ref=0.1;
f_density=1;

#if CHAP
K1=2000;
K2=33;
//nu = 68
K=1.0;  
M=1.0;
//bzero=1 always, this is hardcoded into the code. 
//kzero =0.01 
rho_solid_zero=1000;
rho_fluid_zero=1000;
p_fluid_zero=0;      //Doesn't even feature in code at the moment
Kperm=10.00000001;
mchap=0;
#endif

}
#endif

  Real p_fluid;
  Real f_density;
  Real Kperm;


  Real rho_solid_zero, rho_fluid_zero;


Real m;

Real K1 ;
Real K2 ;
  //Mooney Rivlin Law constants - need to make a seperate Mooney Rivlin law
  // ot just pretend k_1 = 0 , and use the current neo hooken (currentley implemented).
  //Real k_1;
  //Real k_2;


  Real f_phase_ref;
  Real f_phase;
  

//Neo hookean stuff
  Real E;
  Real nu;

  //Old KCL stuff
  Real f ;
Real gamma; 

#if CHAP
Real k1 ;
Real k2 ;
Real K ;
Real p_fluid_zero;
Real M;
Real fchap ;
Real fchapd ;
Real fchapdd ;
Real mchap ;


Real calc_fchap(const Real J);
Real calc_fchapd(const Real J);
Real calc_fchapdd(const Real J);
Real calc_mchap();

#endif


  void calculate_stress_poro();
  void calculate_fluid_pressure();
  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) ;
  void init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m, Real p_fluid);

  void c_update(RealTensor C) ;
  void calculate_tangent();



Real A; 
Real D; 
Real Q;


};

#endif /* NONLINEAR_NEOHOOKE_CC_H_ */

//#endif