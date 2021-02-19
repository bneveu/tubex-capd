

#include <iostream>
#include <vector>
#include "ibex.h"
#include "tubex.h"

#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;


#include <capd/capdlib.h>
#include <tubex_capd2tubex.h>
#include <tubex_TubeVectorODE.h>


 TFunction f("x1", "(x1+0)");
 TFunction f1("x1", "(-x1+0)");

void contract(TubeVector& x, double t0, bool incremental)
{
  
 
    // Boundary constraints

    Variable vx0, vx1;
    SystemFactory fac;
    fac.add_var(vx0);
    fac.add_var(vx1);
    fac.add_ctr(sqr(vx0) + sqr(vx1) = 1);
    System sys(fac);
    ibex::CtcHC4 hc4(sys);
    IntervalVector bounds(2);
    bounds[0] = x[0](0.);
    bounds[1] = x[0](1.);
    hc4.contract(bounds);
    x.set(IntervalVector(bounds[0]), 0.);
    x.set(IntervalVector(bounds[1]), 1.);
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
{   float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    /* =========== PARAMETERS =========== */
   
    int n = 1;

    Vector epsilon(n,0.0005);
    //Vector epsilon(n,1.e-12);
    Interval domain(0.,1.);

    TubeVector x(domain, n);
    TrajectoryVector truth1(domain, TFunction("exp(t)/sqrt(1+exp(2))"));
    TrajectoryVector truth2(domain, TFunction("-exp(t)/sqrt(1+exp(2))"));

    /* =========== SOLVER =========== */


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.99);
    // solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.99);

  //  solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(5000);
    //solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(0);
    solver.set_var3b_external_contraction(false);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    cout << l_solutions.front()(domain.ub()).diam() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;
    cout << l_solutions.back()<<" ti-> " <<l_solutions.back()(domain.lb()) << " tf -> "<< l_solutions.back()(domain.ub()) <<" max diam : (" <<l_solutions.back()[0].max_gate_diam(t_max_diam)<<")"<< " volume :  "<< l_solutions.back().volume()<<" ti (diam) -> " <<l_solutions.back()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.back()(domain.ub()).diam() << endl;

    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
}
