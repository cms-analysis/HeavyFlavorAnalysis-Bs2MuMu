/*
 *  mlp_effLinkDef.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.05.12.
 *
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pargma link C++ function uncertainty_mlp(TTree*,pair<double,double>,double);
#pragma link C++ function test_uncertainty();
#pragma link C++ function test_uncertainty(double);

#endif
