#ifndef BAYESIANLIMIT

double posterior(double s,int n,double e0,double esig,
                 double b0,double bsig,double alpha);

double posteriorb0(double s,int n,double e0,double esig,double b,double alpha);

double postint(double su,int n,double e0,double esig,
               double b0,double bsig,double alpha);

double postintb0(double su,int n,double e0,double esig,double b,double alpha);

double postinte0b0(double su,int n,double e,double b,double alpha);

double incompletebeta(double a, double b, double x2);

double incompletegamma(double a, double x);

double blimit(double beta,int n,double e0,double esig,
              double b0,double bsig,double alpha);

double blimitb0(double beta,int n,double e0,double esig,double b,double alpha);

double blimite0b0(double beta,int n,double e,double b,double alpha);

#define BAYESIANLIMIT 1
#endif
