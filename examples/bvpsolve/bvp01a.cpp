//created by bedouhene 09/12/2019
// BvpSolve 2 (4 variants : ksi=0.2, ksi=0.1, ksi=0.01, ksi=0.001

#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"

#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;


TFunction f("x1", "x2" ,"(x2;x1)");
TFunction f1("x1", "x2" ,"(-x2;-x1)");


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

  Interval domain(0.,1);

    TubeVector x(domain, 2);

    IntervalVector v(2);
    v[0]=Interval(1.,1.);

    v[1]=Interval(-100.,100.);
    //    v[1]=Interval(-1.e300,1.e300);
    x.set(v, 0.); // ini
    v[0]=Interval(0.);

    v[1]=Interval(-100,100.);
    //v[1]=Interval(-1.e300,1.e300);


    x.set(v,1.);
    
    double eps0=0.01;
    double eps1=0.01;
    
    /* =========== SOLVER =========== */
    Vector epsilon(2);
    epsilon[0]=eps0;
    epsilon[1]=eps1;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.);
    //solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_fxpt_ratio(0.9);
    solver.set_var3b_propa_fxpt_ratio(0.9);
    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(2000);

    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(0);
    solver.set_var3b_external_contraction(true);
    cout.precision(6);
    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
    //list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x,f);
    std::cerr.rdbuf(OldBuf);
    cout << "nb sol " << l_solutions.size() << endl;
    cout << l_solutions.front() << endl;

    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    
    
    return 0;
}
