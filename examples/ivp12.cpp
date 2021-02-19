//created by bedouhene 09/12/2019

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
TFunction f("x", "(-x^2)");
TFunction f1("x", "(x^2)");

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
  if (!(x.is_empty())){
    Interval domain1= Interval(-domain.ub(),-domain.lb());
    a0=x(domain.ub());

    TubeVector b = TubeVectorODE(domain1,f1,a0,timestep,CAPD_MODE);
    if (x.volume() < DBL_MAX)
      x&=reversetube(b);
    else
      x=x&reversetube(b);
  }
   
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution

    Interval domain(0.,5.);
    //    TubeVector x(domain,0.1, 1);
    TubeVector x(domain, 1);
    IntervalVector v(1);
    v[0]=Interval(0.1,0.4);

    x.set(v, 0.); // ini

    double epsilon=0.1;
    Vector eps(1,epsilon);

    /* =========== SOLVER =========== */
    tubex::Solver solver(eps);

    solver.set_refining_fxpt_ratio(2);
    //    solver.set_propa_fxpt_ratio(0.);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(0.9);
   
    solver.set_var3b_propa_fxpt_ratio(0.9);

    solver.set_var3b_timept(1);
    solver.set_trace(1);
    solver.set_max_slices(5000);
    //    solver.set_max_slices(1);
    solver.set_refining_mode(2);
    solver.set_bisection_timept(-1);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(0);
    solver.set_var3b_external_contraction(true);
    cout << "avant solve " << endl;
    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, f);    
    //list<TubeVector> l_solutions = solver.solve(x, &contract);
    std::cerr.rdbuf(OldBuf); 

    cout << "nb sol " << l_solutions.size() << endl;
    cout << l_solutions.front() << endl;

    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;



    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;

}
