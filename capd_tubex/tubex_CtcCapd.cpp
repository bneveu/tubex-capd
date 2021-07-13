/**
 *  CtcCapd.cpp
 * ----------------------------------------------------------------------------
 *  \date       2021
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2021 Tubex
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */



#include "tubex_Exception.h"
#include "tubex_CtcCapd.h"

using namespace std;
using namespace ibex;
using namespace tubex;


namespace tubex
{
  CtcCapd::CtcCapd ( TFunction& ffwd,  TFunction& fbwd) : ffwd(ffwd), fbwd(fbwd){preserve_slicing(false); capdstringfwd=tubexFnc2capdString(ffwd),capdstringbwd=tubexFnc2capdString(fbwd);}

 
  void CtcCapd::contract (TubeVector& x,double t0, bool incremental, double timestep) {
    if (!(x.is_empty())){
      Interval domain= x[0].tdomain();
      if (! incremental || domain.lb()==t0 || domain.ub()==t0  ){
	if (! incremental || domain.lb()==t0)
	try{
	  IntervalVector a0 = x(domain.lb());

	  const TubeVector& a = capd2tubex(domain,domain, capdstringfwd,a0,timestep);
	  //  cout << " a " << endl;
	  if (m_preserve_slicing || incremental)
	    x&=a;
	  else
	    x=x&a;
	}
	catch(const char* s){ //cout << s << endl;
	}
	catch (exception& e) {// cout << e.what() << endl;  cout << " exit " << endl;
	}
    	//  cout << " x after fwd " << x << endl;
	if (!(x.is_empty()) && (! incremental || domain.ub()==t0)){
	  Interval domain1= Interval(-domain.ub(),-domain.lb());
	  IntervalVector a0 = x(domain.ub());
	  try {
	    const TubeVector& b = capd2tubex(domain1,domain1,capdstringbwd,a0,timestep);
	    if (m_preserve_slicing|| incremental )
	      x&=reversetube(b);
	    else
	      x=x&reversetube(b);
	  }
	  catch(const char* s){ //cout << s << endl;
	  }
	  catch (exception& e) {//cout << e.what() << endl; cout << " exit2 " << endl; }
	  }
	}
      }
      else{

	if (t0 < domain.ub())
	  {
	    Interval domaint0 (t0,domain.ub());    
	    IntervalVector a0 = x(t0);
	    try{
	      //	      cout << " appel capd 1 " << a0 << " domain " << domaint0 << endl;
	      x&=capd2tubex(domain, domaint0,capdstringfwd,a0,timestep);
	      }
	      catch(const char* s){ //cout << s << endl;
	      }
	      catch (exception& e) {//cout << e.what() << endl; cout << " exit2 " << endl; }
	      }
	    }
    
	  if ((!(x.is_empty())) && t0> domain.lb()){
	    Interval domain1= Interval(-domain.ub(),-domain.lb());
	    Interval domain1t0= Interval(-t0,-domain.lb());
	    
	    IntervalVector a0 = x(t0);
	    try{
	      //	      cout << " appel capd 2 " << a0 << " domain " << domain1t0 << endl;

	      x&= reversetube(capd2tubex(domain1,domain1t0,capdstringbwd,a0,timestep));

	    }
	    catch(const char* s){ //cout << s << endl;
	    }
	    catch (exception& e) {//cout << e.what() << endl; cout << " exit2 " << endl; }
	    }
	  }
      }
    }
  }
  
}
