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


   Interval domain(-10,0.);
  double init = 2.*atan(exp(-10)*tan(0.5));
  

  TFunction f("x","(sin(x))"); //function to be integrated

  IntervalVector a0(1,Interval(init,init+1.e-8)); // initial condition for reference tube
  
  double timestep = 0.0;
  TubeVector a = TubeVectorODE(domain,f,a0,timestep,CAPD_MODE);

  TubeVector b = reversetube(a);
  cout << b << " " << b[0].first_slice()->input_gate() << " volume " << b[0].volume() << endl;
    
}


   
   
