//#include "functions.h"
#include <iostream>
#include "elem.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "solid_system.h"
#include "linear_implicit_system.h"
#include "transient_system.h"

using namespace libMesh;
  //We need to find the correct displacemnt for nodes that are on an edge (i.e are not a vertex)
  //To compute the correct poswition we take the two vertex nodes (a,c) and compute the midpoint (b) which will be the desired boundary condition.
  
Point constrain_tet_nodes(EquationSystems& es, const Elem* elem, int n){
  //const MeshBase& aux_mesh = es.get_mesh();
  //Elem* aux_el=aux_mesh.elem(elem->id());

TransientLinearImplicitSystem&  last_non_linear_soln =
        es.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");

std::vector<Point> refNode_points(10);

for (unsigned int n=0; n<elem->n_nodes(); n++){

	Node *node = elem->get_node(n);
	Point p;
	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
      	Real value = last_non_linear_soln.current_solution(source_dof);
      	p(d)=value;
   	}

refNode_points.at(n)=p;
}

  if(n==7){

    Point a = 1*refNode_points[0];
    Point c = 1*refNode_points[3];
    return c+(0.5*a-0.5*c);
  }
  else if(n==4){
    Point a = 1*refNode_points[1];
    Point c = 1*refNode_points[0];
    return c+(0.5*a-0.5*c);
  }
  else if(n==8){
 //        std::cout<<"  n in func " << n << std::endl;

    Point a = 1*refNode_points[1];
    Point c = 1*refNode_points[3];
    Point b=c+(0.5*a-0.5*c);
    return c+(0.5*a-0.5*c);
  }
  else if(n==5){
    Point a = 1*refNode_points[1];
    Point c = 1*refNode_points[2];
    return c+(0.5*a-0.5*c);
  }
  else if(n==9){
    Point a = 1*refNode_points[3];
    Point c = 1*refNode_points[2];
    return c+(0.5*a-0.5*c);
  }
  else if(n==6){
    Point a = 1*refNode_points[0];
    Point c = 1*refNode_points[2];
    return c+(0.5*a-0.5*c);
  }
 else{
   return 1*refNode_points[n];
 }
  
}
