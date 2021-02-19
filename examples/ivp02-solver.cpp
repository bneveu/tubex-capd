//created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"

#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;
#include <capd/capdlib.h>
#include <tubex_capd2tubex.h>
#include <tubex_TubeVectorODE.h>

TFunction f("x","(-sin(x))");
TFunction f1("x","(sin(x))");



void contract(TubeVector& x, double t0, bool incremental)
{
  Interval domain= x[0].tdomain();

  double timestep = 0.;
  IntervalVector a0(x(domain.lb()));
  TubeVector a = TubeVectorODE(domain,f,a0,timestep,CAPD_MODE);

  if (x.volume() < DBL_MAX)
    x&=a;
  else
    x=x&a;

  Interval domain1= Interval(-domain.ub(),-domain.lb());
  a0=x(domain.ub());
  TubeVector b = TubeVectorODE(domain1,f1,a0,timestep,CAPD_MODE);
  if (x.volume() < DBL_MAX)
    x&=reversetube(b);
  else
    x=x&reversetube(b);
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    

    Interval domain(0.,10.);
    
    TubeVector x(domain, 1);
    IntervalVector v(1);
    //    double init = 2.*atan(exp(-10)*tan(0.5));  //initial condition for backtrack 
    double init = 2.*atan(exp(0)*tan(0.5));
    v[0]=Interval(init);

    //    x.set(v, 10.); // init for backward propagation
    x.set(v, 0.); // ini

    double eps=0.002;

    /* =========== SOLVER =========== */
    Vector epsilon(1, eps);


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2);
    solver.set_propa_fxpt_ratio(0.9);
    solver.set_var3b_fxpt_ratio(-1);

   // solver.set_var3b_timept(0);
    solver.set_trace(1);
    //    solver.set_max_slices(1);
    solver.set_max_slices(10000);
    //solver.set_refining_mode(0.9);
    solver.set_bisection_timept(0);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(0);
    
    //list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, f);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;

    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

   
    return 0;
}
