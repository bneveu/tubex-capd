/** 
 *  tubex-lib - Examples

 */

#include <tubex.h>
//#include <tubex-rob.h>
#include <capd/capdlib.h>
#include <tubex_capd2tubex.h>
#include <tubex_TubeVectorODE.h>

using namespace std;
using namespace tubex;
using namespace vibes;


int main()
{
    // ----- Generate reference tube thanks to ODE integration -----

  Interval domain(0,10.);

  
  /* 
  TFunction f("x","y","(-sin(x);1)"); //function to be integrated
  
  IntervalVector a0(2);
  a0[0]=Interval(1.0,1.0); // initial condition for reference tube
  a0[1]=Interval(0.0,0.0);
  */
  
  TFunction f("x","(-sin(x))"); //function to be integrated
  
  IntervalVector a0(1,Interval(1.0,1.0)); // initial condition for reference tube
  
  
  double timestep = 0.;
  TubeVector a = TubeVectorODE(domain,f,a0,timestep,CAPD_MODE);
  cout << a << " " << a[0].last_slice()->output_gate() << " volume " << a[0].volume() <<   endl;
}


   
   
