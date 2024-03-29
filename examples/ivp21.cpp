//created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"
#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;

TFunction f("x1", "x2" ,"(-x2+0.1*x1*(1-x1^2-x2^2);x1+0.1*x2*(1-x1^2-x2^2))");    
TFunction f1("x1", "x2" ,"(x2-0.1*x1*(1-x1^2-x2^2);-x1-0.1*x2*(1-x1^2-x2^2))");    

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

    TFunction f("x1", "x2" ,"(-x2+0.1*x1*(1-x1^2-x2^2);x1+0.1*x2*(1-x1^2-x2^2))");    
    Interval domain(0.,5.);

    TubeVector x(domain, 2);
    IntervalVector v(2);
    v[0]=Interval(0.7,1.3);
    v[1]=Interval(0.,0.);

    x.set(v, 0.); // ini

    double eps=0.05;

   /* =========== SOLVER =========== */
    Vector epsilon(2, eps);


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.8);
    //    solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.999);
   // solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(1000);
    //solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(-1);
    solver.set_contraction_mode(4);
    solver.set_var3b_external_contraction(true);
    solver.set_stopping_mode(0);
    cout << " avant solve " << x[0](0) << endl;
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x,  &contract);
    // list<TubeVector> l_solutions = solver.solve(x, f);
    cout << l_solutions.front()(domain.ub()) << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam) <<" , "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<< " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


   
    return 0;
    

}
