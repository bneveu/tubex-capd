/** 
 *  CAPD to Tubex
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Julien Damers  Bertrand Neveu
 *  \copyright  Copyright 2020 Tubex
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "tubex_capd2tubex.h"
#include "tubex_Exception.h"
#include "capd/capdlib.h"

using namespace std;
using namespace ibex;
using namespace tubex;

namespace tubex
{

    string tubexFnc2capdString(const TFunction& f)
    {

        int a_capd_dim = f.nb_vars();
        int a_ibex_dim = a_capd_dim+1;

	ibex::IntervalVector a_ibex(a_ibex_dim);
        vector<ibex::IntervalVector> ibex_curve;

        // Generate the string that CAPD will process to compute
        string capd_string ="time:t;var:";
        string function_string = f.expr().substr(1, f.expr().size() - 2); // removing outside parentheses
        const char comma = ',';
        const char semicolon = ';';
        replace(function_string.begin(),function_string.end(),semicolon,comma);
        for(int i=0;i<a_capd_dim-1; i++)
        {
            capd_string.append(f.arg_name(i));
            capd_string.append(",");

        }
        capd_string.append(f.arg_name(a_capd_dim-1));
        capd_string.append(";fun:");
        try
        {
            if (function_string.size()>2)
            {
                capd_string.append(function_string);
                capd_string.append(";");

            }
        }
        catch(exception& e)
        {
	  cout << "\n\nException caught!\n" << e.what() << endl << endl;
        }

        return (capd_string);

    }



  vector<ibex::IntervalVector> capd2ibex(const ibex::Interval& domain, capd::IMap& vectorField, const ibex::IntervalVector& x0, vector<ibex::IntervalVector>& gates, vector<double>& gatetimes,const double& timestep )
    {
        int a_capd_dim = x0.size();

        vector<ibex::IntervalVector> ibex_curve;
	gates.push_back(x0);
	gatetimes.push_back(domain.lb());

        try
        {
            // CAPD processing

            // The solver uses high order enclosure method to verify the existence
            // of the solution.
            // The order will be set to 20.
	  capd::dynsys::ILastTermsStepControl SCP(1, 1e-6);
	  capd::IOdeSolver solver(vectorField,20,SCP);

            // Set a fixed integration step if needed
            if (timestep!=0)
            {
                solver.setStep(timestep);
            }
            capd::ITimeMap timeMap(solver);

            // This is our initial condition
            capd::IVector a_capd(a_capd_dim);
            for (int i = 0; i<a_capd_dim; i++)
            {
                a_capd[i] = capd::interval(x0[i].lb(),x0[i].ub());
            }
            // define a doubleton representation of the interval vector x
	    capd::C0Rect2Set s(a_capd,domain.lb());
	    //capd::C0HORect2Set s(a_capd,domain.lb());

            // Here we start to integrate.
            double T=domain.ub();
            timeMap.stopAfterStep(true);
	    
            do
            {
	      capd::IVector result= timeMap(T,s);
	      ibex::IntervalVector vresult(x0.size());
	      for (int i=0;i<x0.size(); i++)
		vresult[i]=Interval (result[i].left().leftBound(),result[i].right().rightBound() );
	      gates.push_back(vresult);
	      capd::IVector enclosure = s.getLastEnclosure();
	      //	      cout << " enclosure " << enclosure << endl;
	      ibex::IntervalVector ibexenclosure(x0.size());
	      for (int i=0;i<x0.size(); i++)
		ibexenclosure[i]=Interval(enclosure[i].left().leftBound(),enclosure[i].right().rightBound() );
	      //	      cout << " ibex enclosure " << ibexenclosure << endl;
	      ibex_curve.push_back(ibexenclosure);
	      if (timeMap.getCurrentTime().right().rightBound() < domain.ub())
		gatetimes.push_back(timeMap.getCurrentTime().right().rightBound());
	      else
		gatetimes.push_back(domain.ub());
	      //      cout << " result " << result << " time " << timeMap.getCurrentTime() << endl;
	    
            }while(!timeMap.completed());


        }

        catch(exception& e)
        {
	  // cout << "e1 " <<"\n\nException caught!\n" << e.what() << endl << endl;
	  ibex_curve.clear();
	  throw Exception("capd2ibex", e.what());

        }
	
        return(ibex_curve);
    }
  

  TubeVector ibex2tubex(vector<ibex::IntervalVector>& ibex_curve, vector<ibex::IntervalVector> & gates , vector<double> & gatetimes){
  // Creating our Tube Vector object thanks to the info previously stored
    vector<ibex::Interval> vdomains;
    for (int i=0 ; i< gatetimes.size()-1; i++)
      vdomains.push_back(Interval(gatetimes[i],gatetimes[i+1]));
    
    TubeVector output_tube(vdomains,ibex_curve);
    // Sampling the tubevector with the gates calculated by capd
    for (int i =0; i< gates.size(); i++){
      //cout << i << "  " << gatetimes[i] << "  " << gates[i] << endl;
    output_tube.sample(gatetimes[i],gates[i]);
  }
  return(output_tube);
  }
  

  TubeVector capd2tubex(const Interval& domain, const TFunction& f, const IntervalVector& x0, const double timestep){
    string s =tubexFnc2capdString(f);
    return capd2tubex(domain , s , x0, timestep);
  }

  TubeVector capd2tubex(const Interval& domain, const string& capd_string, const IntervalVector& x0, const double timestep)
  {
    try
        {
	  capd::IMap vectorField(capd_string);
	  //cout << " after capd_string " << endl;
	  vector<double> gatetimes;
	  vector<IntervalVector> gates;
	  vector<IntervalVector> ibex_curve = capd2ibex(domain, vectorField, x0, gates,gatetimes,timestep);
	  //	    cout << " ibex curve " <<  ibex_curve.size() << endl;

	  if (!ibex_curve.empty())
	    return(ibex2tubex(ibex_curve,gates, gatetimes));
	  else
	    {
	      TubeVector empt(domain,x0.size());
	      empt.set_empty();
	      return empt;
	    }
	  
        }
    catch(exception& e)
      {
	throw Exception("capdt2ubex", e.what());
      }
    
  }



  TubeVector reversetube(const TubeVector & tubevector){
    Interval domain1 (-tubevector[0].tdomain().ub(), -tubevector[0].tdomain().lb());
    TubeVector tubevector1(domain1, tubevector.size());
    for (int i=0;i< tubevector.size(); i++)
      tubevector1[i]=reversetube(tubevector[i]);
    return tubevector1;  

  }

  Tube reversetube(const Tube & tube){
    Tube tube1(Interval(-tube.tdomain().ub(), -tube.tdomain().lb()));
    for(const Slice *s = tube.first_slice() ; s != NULL ; s = s->next_slice()){
      tube1.sample(-(s->tdomain().lb()),s->input_gate());
      tube1.first_slice()->set_envelope(s->codomain());
    
    }
    tube1.first_slice()->set_input_gate(tube.last_slice()->output_gate());
    return tube1;
  }

}



  


