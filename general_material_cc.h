#ifndef GENERAL_MATERIAL_H_
#define GENERAL_MATERIAL_H_

#include <dense_vector.h>
#include <dense_matrix.h>
#include <vector_value.h>
#include <tensor_value.h>
#include <getpot.h>


#include "defines.h"

<<<<<<< HEAD
 

class GeneralMaterialConfig {
public:

/*
=======
/**
 * This class implements a constitutive formulation for an Neo-Hookean elastic solid
 * in terms of the current configuration. This implementation is suitable for computing
 * large deformations. See e.g. "Nonlinear finite element methods" (P. Wriggers, 2008,
 * Springer) for details.
 */

/*  v=nu
 * lambda=Ev/((1+v)(1-2v))
 */

//#if COMPRESSIBLE 
class GeneralMaterialConfig {
public:
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
  GeneralMaterialConfig(
      const std::vector<std::vector<RealGradient> >& dphi) :
      dphi(dphi) {

  }
<<<<<<< HEAD
*/

  GeneralMaterialConfig(
      const std::vector<std::vector<RealGradient> >& dphi,   const std::vector<std::vector<Real> >& psi) :
      dphi(dphi),psi(psi) {
  }

=======
//#endif

//#if INCOMPRESSIBLE || (CHAP) 
//class GeneralMaterialConfig {
//public:
  GeneralMaterialConfig(
      const std::vector<std::vector<RealGradient> >& dphi,   const std::vector<std::vector<Real> > psi) :
      dphi(dphi),psi(psi) {
  }

//#endif

 

  /**
   * Initialize the class for the given displacement gradient at the
   * specified quadrature point.
   */
  //void init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp);
  //void init_for_qp(VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m);

  /**
   * Return the residual vector for the current state.
   */
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76

  void get_residual(DenseVector<Real> & residuum, unsigned int & i);

  /**
   * Return the stiffness matrix for the current state.
   */
  void get_linearized_stiffness(DenseMatrix<Real> & stiffness,
      unsigned int & i, unsigned int & j);

<<<<<<< HEAD

  void get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i);

=======
  #if INCOMPRESSIBLE || (CHAP) 
  void get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i);
  #endif
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
  
  /**
   * Flag to indicate if it is necessary to calculate values for stiffness
   * matrix during initialization.
   */
  bool calculate_linearized_stiffness;
//private:
  void build_b_0_mat(int i, DenseMatrix<Real>& b_l_mat);

 // void calculate_stress();
 // void calculate_tangent();
  static void tensor_to_voigt(const RealTensor &tensor, DenseVector<Real> &vec);
  static void tensorOtensor_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);
    static void z_ref_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);
 
  void c_update(RealTensor C) ;


  unsigned int current_qp;
<<<<<<< HEAD



  Number p_current;
  Real p_solid;


  const std::vector<std::vector<RealGradient> >& dphi;
  const std::vector<std::vector<Real> >& psi;


  
=======
  const std::vector<std::vector<RealGradient> >& dphi;







  
  //#if INCOMPRESSIBLE || (CHAP)  
  Number p_current;
  Real p_solid;
  const std::vector<std::vector<Real> > psi;
 // #endif     
  
>>>>>>> 2401f45059a0293beb4a22be9a802b731c757b76
    DenseMatrix<Real> C_mat;

  Real I_1, I_2, I_3, J;
  RealTensor F, S, tau, sigma, C, invC, b, Ft, Identity;
  DenseMatrix<Real> B_L;
  DenseMatrix<Real> B_K;
};

#endif /* GENERAL_MATERIAL_H_ */

