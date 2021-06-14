//created by neveu dec 11 2020
// problem Bvpsolve27 (xi=0.1)


#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"
#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;

  TFunction f("x1", "x2" ,"(x2;10*x1*(1-x2))");
  TFunction f1("x1", "x2" ,"(-x2;-10*x1*(1-x2))");
void contract(TubeVector& x, double t0, bool incremental)
{CtcCapd ctccapd(f,f1);
  if (x.volume() < DBL_MAX && x.nb_slices() > 1)
    ctccapd.preserve_slicing(true);
  else
    ctccapd.preserve_slicing(false);
  ctccapd.contract (x, t0, incremental);
}

int main() {

   
    /* =========== PARAMETERS =========== */

    Interval domain(0.,1.);

    //TubeVector x(domain,2);
    TubeVector x(domain,IntervalVector(2,Interval(-1000,1000)));

    IntervalVector v(2);
    v[0]=Interval(1.);
    v[1]=Interval(-100.,100.);

    
    x.set(v, 0.); // ini
    v[0]=Interval(1./3.);
    v[1]=Interval(-100.,100.);

    x.set(v,1.);
    
    
    
    double eps0=0.1;
    double eps1=0.1;
    
   

    /* =========== SOLVER =========== */
    Vector epsilon(2);
    epsilon[0]=eps0;
    epsilon[1]=eps1;

    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(0.99);
    //    solver.set_var3b_fxpt_ratio(0.9);

    solver.set_var3b_propa_fxpt_ratio(0.99);
    

    solver.set_var3b_timept(2);
    solver.set_trace(1);
    solver.set_max_slices(2000);

    solver.set_bisection_timept(3);

    solver.set_refining_mode(2);
    solver.set_stopping_mode(0);
    solver.set_contraction_mode(4);
    solver.set_var3b_external_contraction(true);
    
    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, f);
    std::cerr.rdbuf(OldBuf);
    
    cout << "nb sol " << l_solutions.size() << endl;
    if(l_solutions.size()>0){
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max gate diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

    }

    
    return 0;
}
