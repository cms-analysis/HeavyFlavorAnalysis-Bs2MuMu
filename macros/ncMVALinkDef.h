/*
 *  ncMVALinkDef.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.04.12.
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ function ncEvalAll(TTree*);
#pragma link C++ function ncEvalAll(const char*);

#pragma link C++ function ncRunTraining(TTree*,double,TTree*,double);
#pragma link C++ function ncRunDefaultTraining(const char*, const char*);

#endif
