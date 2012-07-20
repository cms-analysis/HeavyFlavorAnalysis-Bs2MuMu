/*
 *  NCFunctionsLinkDef.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 8.3.10.
 *
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ function f_p1(double*,double*);
#pragma link C++ function f_p2(double*,double*);

#pragma link C++ function f_expo(double*,double*);
#pragma link C++ function f_Gauss(double*,double*);
#pragma link C++ function f_2G(double*,double*);

#pragma link C++ function f_gauss(double*,double*);
#pragma link C++ function f_2g(double*,double*);
#pragma link C++ function f_3g(double*,double*);
#pragma link C++ function f_double_gauss(double*,double*);

#pragma link C++ function f_cb(double*,double*);
#pragma link C++ function f_p1acb(double*,double*);
#pragma link C++ function f_p2acb(double*,double*);

#pragma link C++ function f_fnov(double*,double*);

#pragma link C++ function f_argus(double*,double*);
#pragma link C++ function f_aag(double*,double*);
#pragma link C++ function f_exparggau(double*,double*);
#pragma link C++ function f_aacb(double*,double*);

#pragma link C++ function f_p0ag(double*,double*);
#pragma link C++ fucntion f_p1ag(double*,double*);
#pragma link C++ fucntion f_p2ag(double*,double*);
#pragma link C++ function f_p0a2g(double*,double*);

#pargma link C++ function f2_chi2ellipsis(double*,double*);

#pragma link C++ function f_gauss_expo(double*,double*);
#pragma link C++ function f_gauss_linear(double*,double*);
#pragma link C++ function f_double_gauss_linear(double*,double*);

#pragma link C++ function f_boltzmann(double*,double*);

#pragma link C++ function f_charact(double*,double*);
#pragma link C++ function f_skewnormal(double*,double*);

#endif
