/** 
 *  CAPD to Tubex
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Julien Damers
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
        int a_ibex_dim = a_capd_dim+1;

	ibex::IntervalVector a_ibex(a_ibex_dim);
        vector<ibex::IntervalVector> ibex_curve;


        try
        {
            // CAPD processing

            // The solver uses high order enclosure method to verify the existence
            // of the solution.
            // The order will be set to 20.
	  capd::dynsys::ILastTermsStepControl SCP(1, 1e-4);
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
	    //	    capd::C0HORect2Set s(a_capd,domain.lb());

            // Here we start to integrate.
            double T=domain.ub();
            timeMap.stopAfterStep(true);
            capd::interval prevTime(domain.lb());
	    gates.push_back(x0);
	    gatetimes.push_back(domain.lb());
	    
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
  

  TubeVector ibex2tubex(vector<ibex::IntervalVector>& ibex_curve, vector<ibex::IntervalVector> & gates
		      , vector<double> & gatetimes){
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
  
  TubeVector ibex2tubex2(vector<ibex::IntervalVector>& ibex_curve, vector<ibex::IntervalVector> & gates
			, vector<double> & gatetimes)
    {
        /*
         * To create a tube we need to assert that for each slice, the end of the slice
         * is equal to the beginning of the next slice. This function first verify if
         * for each box of the ibex_curve the condition mentioned above is verified. If
         * it is the box is directly converted into a slice. If not, the current box and the next
         * one are processed to verify this condition.
         */

        /*
         * Initialize two vectors that will store time domain and the interval on the other dimensions
         * for each slice. v_domains[i] + v_codomains[i] = slice[i] of the tube
         */
      vector<ibex::Interval> v_domains; // Store the time interval of each slice of the tube
      vector<ibex::IntervalVector> v_codomains; // Store the other dimensions of the slice of the tube
      vector<ibex::IntervalVector> v_gates; // Store the other dimensions of the slice of the tube
      // cout << ibex_curve[0].size() << endl;
      ibex::IntervalVector current_box(ibex_curve[0].size());
      ibex::IntervalVector mid_box(ibex_curve[0].size());
      
        for(size_t i=0; i<ibex_curve.size()-1;i++)
        {   
	  ibex::IntervalVector current_box = ibex_curve[i];
	  //cout << i << " current_box " << current_box << endl;
	  //    if((current_box[0].ub() == ibex_curve[i+1][0].lb()))
	  //            {
                /*
                 * |###############||###########|
                 * | current_box[0]||next box[0]|
                 * |###############||###########|
                 *
                 *  condition verified? Ok
                 *
                 * --> convert to slice directly
                 */
                v_domains.push_back(current_box[0]);
                v_codomains.push_back(current_box.subvector(1,current_box.size()-1));

		//            }
		//            else
		//            {
                /*
                 * |##########################|
                 * |      current_box[0]      |
                 * |##########################|
                 *                 |#######################|
                 *                 |      next box[0]      |
                 *                 |#######################|
                 *
                 *  condition verified? No
                 *
                 * --> Truncating boxes and creating a link in the middle
                 * |##############||##########||###########|
                 * |current_box[0]||mid_box[0]||next box[0]|
                 * |##############||##########||###########|
                 *
                 * --> convert current_box and mid_box to slices



                mid_box[0] = current_box[0]&ibex_curve[i+1][0]; // storing time dimension of the mid box
                for (int k=1;k<mid_box.size();k++)
                {
                    // Ensure that we keep the convex hull so that our computation is guaranteed (may add pessimism)
                    mid_box[k] = current_box[k]|ibex_curve[i+1][k];

                }
                current_box[0] = Interval(current_box[0].lb(), mid_box[0].lb());
                ibex_curve[i+1][0] = Interval(mid_box[0].ub(),ibex_curve[i+1][0].ub());

                v_domains.push_back(current_box[0]);
                v_codomains.push_back(current_box.subvector(1,current_box.size()-1));
                v_domains.push_back(mid_box[0]);
                v_codomains.push_back(mid_box.subvector(1,current_box.size()-1));
            }
		*/
        }

	
	
        v_domains.push_back(ibex_curve[ibex_curve.size()-1][0]);

	if (v_domains[v_domains.size()-1].ub() > gatetimes[gatetimes.size()-1])
	  v_domains[v_domains.size()-1]=Interval(v_domains[v_domains.size()-1].lb(),gatetimes[gatetimes.size()-1]);
        v_codomains.push_back(ibex_curve[ibex_curve.size()-1].subvector(1,ibex_curve[ibex_curve.size()-1].size()-1));

        // Creating our Tube Vector object thanks to the info previously stored
        TubeVector output_tube(v_domains,v_codomains);
	// Sampling the tubevector with the gates calculated by capd
	for (int i =0; i< gates.size(); i++){
	  //cout << i << "  " << gatetimes[i] << "  " << gates[i] << endl;
	  output_tube.sample(gatetimes[i],gates[i]);
	}
        return(output_tube);
    }


    TubeVector capd2tubex(const Interval& domain, const TFunction& f, const IntervalVector& x0, const double timestep)
    {
        string capd_string = tubexFnc2capdString(f);
        try
        {
	  //	  cout << " before capd_string " << capd_string << endl;
	  
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
            throw Exception("capd2tubex", e.what());
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



  


