//created by neveu sept 30 2020
// problem Bvpsolve32 (xi=1)

#include <iostream>
#include <vector>
#include "tubex.h"
#include "tubex-solve.h"
#include <tubex_CtcCapd.h>

using namespace std;
using namespace ibex;
using namespace tubex;
TFunction f("x1", "x2" ,"x3", "x4", "(x2;x3;x4;10*(x2*x3-x1*x4))");
TFunction f1("x1", "x2" ,"x3", "x4","(-x2;-x3;-x4;-10*(x2*x3-x1*x4))");

void contract(TubeVector& x, double t0, bool incremental)
{CtcCapd ctccapd(f,f1);
  if (x.volume() < DBL_MAX && x.nb_slices() > 1)
    ctccapd.preserve_slicing(true);
  else
    ctccapd.preserve_slicing(false);
  ctccapd.contract (x, t0, incremental);
}
  
    
int main() {
  
  //  TFunction f("x1", "x2" ,"x3", "x4", "(x2;x3;x4;4.6415888336128*(x2*x3-x1*x4))");

   

    /* =========== PARAMETERS =========== */

    Interval domain(0.,1.);
   

    IntervalVector bounds (4);

    bounds[0]=Interval(0,2);
    bounds[1]=Interval(0,2);

    bounds[2]=Interval(-20,20);
    bounds[3]=Interval(-50,0);
    TubeVector x(domain,bounds);

    //    TubeVector x(domain,IntervalVector(4,Interval(-100,100)));
    //TubeVector x(domain,IntervalVector(4,Interval(-50,50)));

    IntervalVector v(4);
    v[0]=Interval(0);
    v[1]=Interval(0);
    v[2]=Interval(0,100);
    v[3]=Interval(-100,0);
    x.set(v, 0.); // ini
    v[0]=Interval(1);
    v[1]=Interval(0);
    v[2]=Interval(-100,0);
    v[3]=Interval(-100,0);
    
    x.set(v,1.);
    
    
    
    double eps0=0.1; 
    double eps1=0.1;    

    /* =========== SOLVER =========== */
    Vector epsilon(4);
    epsilon[0]=eps0;
    epsilon[1]=eps0;
    epsilon[2]=eps1;
    epsilon[3]=eps1;

    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.0);

    solver.set_propa_fxpt_ratio(0.99);
    //solver.set_var3b_fxpt_ratio(-1);

    solver.set_var3b_fxpt_ratio(0.9999);

    solver.set_var3b_propa_fxpt_ratio(0.9999);

    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(10000);
    

    solver.set_bisection_timept(3);

    solver.set_refining_mode(2);
    solver.set_stopping_mode(0);
    solver.set_contraction_mode(4);
    solver.set_var3b_external_contraction(true);
    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
    //    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x, f);
    std::cerr.rdbuf(OldBuf);
    
    cout << "nb sol " << l_solutions.size() << endl;
    if(l_solutions.size()>0){
    double t_max_diam;

    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max gate diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<<", "<< l_solutions.front()[2].max_gate_diam(t_max_diam)<< ", "<<l_solutions.front()[3].max_gate_diam(t_max_diam)<<")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

    }

    
    return 0;
}
