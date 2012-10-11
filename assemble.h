#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_


// The definition of a geometric element
#include "elem.h"
#include "assemble.h"
//#include "nonlinear_neohooke_cc.h"
//#include "solid_system.h"

#include "math.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void assemble_newton (EquationSystems& es,
                      const std::string& system_name);

void assemble_bcs (EquationSystems& es);
void assemble_fluid_bcs (EquationSystems& es);

void assemble_pressure_gard (EquationSystems& es,
                      const std::string& system_name);


void test(int a);

void verify_jack(EquationSystems& es);

void setup_es(EquationSystems& equation_systems);

void read_options(unsigned int &  n_timesteps, std::string& result_file_name,const char* & plot_out_file_name, int& argc, char**& argv) ;


Point get_expanding_sphere_bcs(EquationSystems& es, const Elem* elem, int n,double scale);
Point constrain_tet_nodes(EquationSystems& es, const Elem* elem, int n);



  void get_traction(DenseVector<Real> & traction, Point rX, const Real progress);

  void get_traction_current(DenseVector<Real> & traction, Point rX, const Real progress);

  void get_traction_test(DenseVector<Real> & traction, Point rX, const Real progress);


  void get_bodyforce(DenseVector<Real> & body_force, Point rX, const Real progress);

  void tensor_mult_vector(DenseVector<Real> & ans, RealTensor tens, Point normal);


#endif 
