//Created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"
#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;

TFunction f("x1", "x2", "(x2;-x1)");
TFunction f1("x1", "x2", "(-x2;x1)");
    

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

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution

    double pi=M_PI;
    int n = 2;

    Interval domain(0.,pi/2);
    TubeVector x(domain, n);
    IntervalVector x0(2);
    IntervalVector x1(2);
    x0[0]=Interval(0.,0.);
    x0[1]=Interval(-1.e8,1.e8);

    x.set(x0,0.);
    x1[0]=Interval(2.,2.);
    x1[1]=Interval(-1.e8,1.e8);
    x.set(x1,pi/2);


    /* =========== SOLVER =========== */

    Vector epsilon(n, 0.0005);
    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.);
    // solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    //solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(40000);
    //    solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(0);
    solver.set_var3b_external_contraction(true);
    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);

        cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<" , " <<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

        
    return 0;
}
