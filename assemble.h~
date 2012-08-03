#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_


// The definition of a geometric element
#include "elem.h"
#include "assemble.h"
//#include "nonlinear_neohooke_cc.h"
#include "solid_system.h"

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


Point get_expanding_sphere_bcs(EquationSystems& es, const Elem* elem, int n,double scale);
Point constrain_tet_nodes(EquationSystems& es, const Elem* elem, int n);

#endif 
