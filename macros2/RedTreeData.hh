#ifndef REDTREEDATA
#define REDTREEDATA

struct RedTreeData {
  Long64_t run, evt; 
  int ls, tm, pr, procid, pvn, rr;
  bool json, hlt, cb;
  double bdt, bdt2, pvw8;
  //  npv

  bool gmuid, gmupt, gmueta, gtqual, gtpt, gteta;
  double w8mu, w8tr;

  double pvlip, pvlips, pvlip2, pvlips2, pvip, pvips;

  int q, type;
  double pt, eta, phi, tau, m, cm, cosa, alpha, iso;
  int isotrk, closetrk; 
  double chi2, dof, pchi2dof, fls3d, fl3d, flxy, fl3dE, flsxy, docatrk, docatrkbdt, maxdoca, lip, lipE, tip, tipE;
  
  double osiso, osreliso, osmpt, osmptrel, osmdr; 

  int m1q, m2q;  
  double m1pt, m1eta, m1phi, m1ip, m1chi2;
  double m2pt, m2eta, m2phi, m2ip, m2chi2;
  double kpt, keta, kphi; 
  double k1pt, k1eta, k1phi, k2pt, k2eta, k2phi; 

  bool m1id, m2id; 
  int m1gt, m2gt, k1gt, k2gt; 
  int m1pix, m1bpix, m1bpixl1, m1pv; 
  int m2pix, m2bpix, m2bpixl1, m2pv; 

  double mudist, mudeltar; 
  double g1pt, g2pt, g3pt, g4pt, g1eta, g2eta, g3eta, g4eta, gtau;
  double t1pt, t1eta, t2pt, t2eta, t3pt, t3eta, t4pt, t4eta; 

  double mpsi, mkk;
  double psipt, psieta, psiphi, phipt, phieta, phiphi, dr;
  double md0, dm, ptd0; 

  double hm1pt, hm1eta, hm1phi, hm2pt, hm2eta, hm2phi; 
};

#endif
