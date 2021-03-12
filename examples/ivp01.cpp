/** 
 *  tubex-lib - Examples

 */

#include <tubex.h>
#include <tubex_capd2tubex.h>
#include <tubex_TubeVectorODE.h>
#include <tubex_CtcCapd.h>

using namespace std;
using namespace tubex;
using namespace vibes;


int main()
{
    // ----- Generate reference tube thanks to ODE integration -----

  
 
  TFunction f("x","(x+0)"); //function to be integrated

  Interval domain (-10,0);
  IntervalVector a0(1,Interval(exp(-10),exp(-10))); // initial condition for reference tube
  
  
  double timestep = 0;
  TubeVector a = TubeVectorODE(domain,f,a0,timestep,CAPD_MODE);
  TubeVector b=	reversetube(a);	    
  cout << b << " " << b[0].first_slice()->input_gate() << " volume " << a[0].volume() <<   endl;
}


   
   
