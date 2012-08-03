#ifndef poro_elastic_CC_H_
#define poro_elastic_CC_H_

#include <dense_vector.h>
#include <dense_matrix.h>
#include <vector_value.h>
#include <tensor_value.h>
#include <getpot.h>


#include "defines.h"
#include "nonlinear_neohooke_cc.h"

#if PORO
/**
 * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
 * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
 */
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

class PoroelasticConfig : public NonlinearNeoHookeCurrentConfig {

public:
#if INCOMPRESSIBLE || (CHAP) 
PoroelasticConfig(const std::vector<std::vector<RealGradient> >& dphi, const std::vector<std::vector<Real> >& psi) :  NonlinearNeoHookeCurrentConfig(dphi, psi){
f_phase_ref=0.6;
f_density=1;
A=1.0; 
D=1.0; 
Q=1.0; 
#if CHAP
k1=1.0;
k2=1.0;
K=1.0;
M=1.0;
p_fluid_zero=0;
mchap=0;

#endif
}

PoroelasticConfig(const std::vector<std::vector<RealGradient> >& dphi, const std::vector<std::vector<Real> >& psi, Real p_fluid) : p_fluid(p_fluid),  NonlinearNeoHookeCurrentConfig(dphi, psi){
f_phase_ref=0.6;
f_density=1;
A=1.0; 
D=1.0; 
Q=1.0; 
#if CHAP
k1=1.0;
k2=1.0;
K=1.0;
M=1.0;
p_fluid_zero=0;
mchap=0;

#endif

}
#endif

  Real p_fluid;
  Real f_density;

  //Mooney Rivlin Law constants - need to make a seperate Mooney Rivlin law
  // ot just pretend k_1 = 0 , and use the current neo hooken (currentley implemented).
  //Real k_1;
  //Real k_2;


  Real f_phase_ref;
  Real f_phase;
  

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

#endif