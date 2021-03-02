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
  CtcCapd::CtcCapd ( TFunction& ffwd,  TFunction& fbwd) : ffwd(ffwd), fbwd(fbwd){preserve_slicing(false);}

 
  void CtcCapd::contract (TubeVector& x ,double t0, bool incremental, double timestep) {
    if (!(x.is_empty())){
      Interval domain= x[0].tdomain();
      if (! incremental || domain.lb()==t0 || domain.ub()==t0  ){

	if (! incremental || domain.lb()==t0)
	try{
	  IntervalVector a0(x(domain.lb()));
	  TubeVector a = capd2tubex(domain,ffwd,a0,timestep);
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
	  IntervalVector a0=x(domain.ub());
	  try {
	    TubeVector b = capd2tubex(domain1,fbwd,a0,timestep);
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
	    //    cout << " a0 " << a0 << endl;
	    try{
	      TubeVector tubea = capd2tubex(domaint0,ffwd,a0,timestep);
	      TubeVector tube1 (domain ,x.size());
	      for (int i=0; i<tube1.size(); i++){
		Slice* s1= tube1[i].first_slice();
		if (t0> domain.lb()){
		  tube1[i].sample(t0,a0[i]); s1=s1->next_slice();
		}
		for(const Slice *s = tubea[i].first_slice() ; s != NULL ; s = s->next_slice()){
		  tube1[i].sample(s->tdomain().ub(),s->output_gate());
		    s1->set_envelope(s->codomain());
		    s1=s1->next_slice();
		  }
		}
		x&=tube1;
		//    cout << " x after fwd " << x << endl;
	      }
	      catch(const char* s){ //cout << s << endl;
	      }
	      catch (exception& e) {//cout << e.what() << endl; cout << " exit2 " << endl; }
	      }
	    }
    
	  if ((!(x.is_empty())) && t0> domain.lb()){
	    Interval domain1= Interval(-t0,-domain.lb());
	    IntervalVector a0 = x(t0);
	    try{
	      //    cout << " a0 bwd " << a0 <<  " domain1 " << domain1 << endl;

	      TubeVector tubeb= reversetube(capd2tubex(domain1,fbwd,a0,timestep));
	      //cout << "t0" << t0 <<  " tubeb " << tubeb << endl;
	      TubeVector tube1 (domain,x.size());

	      for (int i=0; i<tube1.size(); i++){
		Slice* s1=tube1[i].first_slice();
		s1->set_input_gate( tubeb[i].first_slice()->input_gate());
		for(const Slice *s = tubeb[i].first_slice() ; s != NULL ; s = s->next_slice()){
		  tube1[i].sample(s->tdomain().ub(),s->output_gate());
		  s1->set_envelope(s->codomain());
		  s1=s1->next_slice();
		}
	      }
	      x&=tube1;
	    }
	    catch(const char* s){ //cout << s << endl;
	    }
	    catch (exception& e) {//cout << e.what() << endl; cout << " exit2 " << endl; }
	    }
	    //    cout << " x after bwd " << x << endl;
	  }
      }
    }
  }
  
}
