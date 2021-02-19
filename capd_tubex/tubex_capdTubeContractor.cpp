/**
 *  TubeVectorODE class
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2021 Tubex
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#include "tubex_capdTubeContractor.h"
#include "tubex_Exception.h"

using namespace std;
using namespace ibex;
using namespace tubex;


namespace tubex
{
    
  void capdcontract (TubeVector & x, const  TFunction& f,  const TFunction& fbwd ,double t0, bool incremental, double timestep) {
    if (!(x.is_empty())){
      if (! incremental){
	Interval domain= x[0].tdomain();
	try{
	  IntervalVector a0(x(domain.lb()));
	  TubeVector a = capd2tubex(domain,f,a0,timestep);
	  //  cout << " a " << endl;
	  if (x.volume() < DBL_MAX && x.nb_slices()>=2 )
	    x&=a;
	  else
	    x=x&a;
	}
	catch(const char* s){ //cout << s << endl;
	}
	catch (exception& e) {// cout << e.what() << endl;  cout << " exit " << endl;
	}
      
  
	//  cout << " x after fwd " << x << endl;
	if (!(x.is_empty())){
	  Interval domain1= Interval(-domain.ub(),-domain.lb());
	  IntervalVector a0=x(domain.ub());
	  try {
	    TubeVector b = capd2tubex(domain1,fbwd,a0,timestep);
	    if (x.volume() < DBL_MAX)
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
	Interval domain (t0, x[0].tdomain().ub());
	if (t0 < x[0].tdomain().ub())
	  {
    
	    IntervalVector a0 = x(t0);
	    //    cout << " a0 " << a0 << endl;
	    try{
	      TubeVector tubea = capd2tubex(domain,f,a0,timestep);
	      TubeVector tube1 (x[0].tdomain(),x.size());
	      for (int i=0; i<tube1.size(); i++){
		Slice* s1= tube1[i].first_slice();
		if (t0> x[0].tdomain().lb()){
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
    
	  if ((!(x.is_empty())) && t0> x[0].tdomain().lb()){
	    Interval domain1= Interval(-t0,-x[0].tdomain().lb());
	    IntervalVector a0 = x(t0);
	    try{
	      //    cout << " a0 bwd " << a0 <<  " domain1 " << domain1 << endl;
	      TubeVector tubec = capd2tubex(domain1,fbwd,a0,timestep);
	      TubeVector tubeb= reversetube(tubec);
	      //cout << "t0" << t0 << "tubec " << tubec << " tubeb " << tubeb << endl;
	      TubeVector tube1 (x[0].tdomain(),x.size());

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
