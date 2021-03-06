/**
 *  \file
 *  TubeVectorODE class
 * ----------------------------------------------------------------------------
 *  \date       2021
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2020 Tubex
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_CAPDTUBECONTRACTOR_H__
#define __TUBEX_CAPDTUBECONTRACTOR_H__

#include "tubex.h"
#include "tubex_capd2tubex.h"

/*
 A contractor for ODE using CAPD
 It contracts a tube vector x, using the ODE differential function f forward and the ODE fbwd function backward (mirror of f)
 with timestep (by default 0)
 For example  if f  is  TFunction f("x1", "x2" ,"(x2;x1)")
                 fbwd is TFunction fbwd ("x1", "x2" ,"(-x2;-x1)");
 If incremental=false, a forward and a backward integration are performed along the tubevector
 If incremental=true, a forward integration is performed from t0 to the final time of tubevector x 
               and a backward integration is performed from t0 to the initial  time of tubevector x 
 
 */
namespace tubex
{

  class CtcCapd : public DynCtc {
  public : 
    CtcCapd (const TFunction& f, const TFunction& fbwd );
    
  void contract(TubeVector& x, double t0, bool incremental, double timestep=0);
  private :
  void init();
  const TFunction& f;
  const TFunction& fbwd;
  };
  
}


#endif
