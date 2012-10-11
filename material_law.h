#ifndef MATERIAL_LAW_CC_H_
#define MATERIAL_LAW_CC_H_

#include <dense_vector.h>
#include <dense_matrix.h>
#include <vector_value.h>
#include <tensor_value.h>
#include <getpot.h>

//#define LOL 1
//#define INCOMPRESSIBLE 1

#include "defines.h"
#include "defines_two.h"


/**
 * This class implements a constitutive formulation for an Neo-Hookean elastic solid
 * in terms of the current configuration. This implementation is suitable for computing
 * large deformations. See e.g. "Nonlinear finite element methods" (P. Wriggers, 2008,
 * Springer) for details.
 */

/*  v=nu
 * lambda=Ev/((1+v)(1-2v))
 */



class MaterialConfig {
public:
  MaterialConfig(
      const std::vector<std::vector<RealGradient> >& dphi,   const std::vector<std::vector<Real> >& psi) :
      dphi(dphi),psi(psi) {
  }


 


  /**
   * Return the residual vector for the current state.
   */

  
  /**
   * Return the stiffness matrix for the current state.
   */
  void get_linearized_stiffness(DenseMatrix<Real> & stiffness,
      unsigned int & i, unsigned int & j);


  #if INCOMPRESSIBLE  
  void get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i);
  
 #endif
  
  /**
   * Flag to indicate if it is necessary to calculate values for stiffness
   * matrix during initialization.
   */
  bool calculate_linearized_stiffness;
//private:
  void build_b_0_mat(int i, DenseMatrix<Real>& b_l_mat);

 void get_residual(DenseVector<Real> & residuum, unsigned int & i);

void c_update(RealTensor C);
// Should make these functions virtal 
//  void calculate_stress();
//  void calculate_tangent();
   /**
   * Initialize the class for the given displacement gradient at the
   * specified quadrature point.
   */
  //void init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp);
  //void init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m);


  static void tensor_to_voigt(const RealTensor &tensor, DenseVector<Real> &vec);
  static void tensorOtensor_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);
  static void z_ref_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);


  unsigned int current_qp;
  const std::vector<std::vector<RealGradient> >& dphi;

#if INCOMPRESSIBLE || CHAP
  Real p_solid;
#endif

  Number p_current;
  const std::vector<std::vector<Real> >& psi;


RealTensor Identity;
Real I_1, I_2, I_3, J;
DenseMatrix<Real> C_mat;
RealTensor F, S, tau, sigma, C, invC, b, Ft;
DenseMatrix<Real> B_L;
DenseMatrix<Real> B_K;

};

#endif /* MATERIAL_LAW_CC_H_ */

