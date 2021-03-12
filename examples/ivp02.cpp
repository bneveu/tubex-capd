/** 
 *  tubex-lib - Examples

 */

#include <tubex.h>
#include <tubex_capd2tubex.h>
#include <tubex_TubeVectorODE.h>

using namespace std;
using namespace tubex;


int main()
{


  Interval domain(0,10.);

  
 
  TFunction f("x","(-sin(x))"); //function to be integrated
  
  IntervalVector a0(1,Interval(1.0,1.0)); // initial condition for reference tube
  
  
  double timestep = 0.;
  TubeVector a = TubeVectorODE(domain,f,a0,timestep,CAPD_MODE);
  cout << a << " " << a[0].last_slice()->output_gate() << " volume " << a[0].volume() <<   endl;
}


   
   
