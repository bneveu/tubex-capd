//created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"
#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;

TFunction f("y1", "y2", "(-0.7*y1 ; 0.7*y1 - (ln(2)/5.)*y2)");
TFunction f1("y1", "y2", "(0.7*y1 ; -0.7*y1 + (ln(2)/5.)*y2)");





void contract(TubeVector& x, double t0, bool incremental)
{
  CtcCapd ctccapd(f,f1);
  if (x.volume() < DBL_MAX && x.nb_slices() > 1)
    ctccapd.preserve_slicing(true);
  else
    ctccapd.preserve_slicing(false);
  ctccapd.contract (x, t0, incremental);

}

int main()

{   
  TFunction f("y1", "y2", "(-0.7*y1 ; 0.7*y1 - (ln(2)/5.)*y2)");
 

    int n = 2;
    Interval domain(0.,6.);
    TubeVector x(domain,n,IntervalVector(2,Interval(-1e5,1e5))); // initial tube in CP paper
    //TubeVector x(domain,n);


    // Boundary condition:
    IntervalVector init = x(x.tdomain().lb());
    init[0] = Interval(1.25);
    // init[1] = Interval(1.1,1.3);
    //init[1]=Interval(-0.413732, -0.372333);
    x.set(init, 0.);

    // Additional restriction (maximum value):
    Interval domain_restriction(1.,3.);
    IntervalVector max_restriction(2);
    max_restriction[0] = Interval(-1e5,1e5);  // restriction in CP paper
    //    max_restriction[0] = Interval(-1.e307,1.e307); // the diameter must be bounded somewhere
    max_restriction[1] = Interval(1.1,1.3);


    x.set(max_restriction,domain_restriction);
    /* =========== SOLVER =========== */

//    contract(x);

    Vector epsilon(n);
    epsilon[0] = 0.04;
    epsilon[1] = 0.04;
//
//
    tubex::Solver solver(epsilon);
//
    solver.set_refining_fxpt_ratio(2);
    //    solver.set_propa_fxpt_ratio(0.99);
    solver.set_propa_fxpt_ratio(0.);

    solver.set_var3b_fxpt_ratio(0.9);
    //solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.9);
    solver.set_var3b_external_contraction(true);
//
  //  solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(10000);
    //    solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(0);

    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    // list<TubeVector> l_solutions = solver.solve(x,&contract);
    //list<TubeVector> l_solutions = solver.solve(x,f);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    
    return 0;
    
  
}
