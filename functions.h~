#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

void assemble_newton (EquationSystems& es,
                      const std::string& system_name);

void test(int a);

void test2(int a);

void get_expanding_sphere_bcs(int n);

// save the undeformed mesh to an auxiliary system
//void save_initial_mesh(EquationSystems& es);


double diffclock(clock_t clock1,clock_t clock2);


void compute_fluid_input(EquationSystems& es, unsigned int qp,const Elem* elem, Real m_old, Real& p_fluid, Real& J);

void print_resluts(EquationSystems& es,Real time,Real l);


#endif 
