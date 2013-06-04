#ifndef REDTREEDATA
#define REDTREEDATA

struct RedTreeData {
  Long64_t run, evt; 
  int ls, tm, pr, procid, pvn, rr;
  bool json, hlt, hltm, hltm2, cb;
  double bdt, bdt2, pvw8;

  bool gmuid, gmutmid, gmumvaid, gmupt, gmueta, gtqual, gtpt, gteta;

  double pvlip, pvlips, pvlip2, pvlips2, pvip, pvips, pvip3d, pvips3d;

  int q, type;
  double pt, eta, phi, tau, m, me, cm, m3, m4, cosa, alpha, iso;
  int isotrk, closetrk, closetrks1, closetrks2, closetrks3; 
  double chi2, dof, chi2dof, pchi2dof, fls3d, fl3d, flxy, fl3dE, flsxy, docatrk, docatrkbdt, maxdoca, lip, lipE, tip, tipE;
  
  double osiso, osreliso, osmpt, osmptrel, osmdr; 

  int m1q, m2q;  
  double m1pt, m1eta, m1phi, m1ip, m1chi2;
  double m2pt, m2eta, m2phi, m2ip, m2chi2;
  double kpt, keta, kphi; 
  double k1pt, k1eta, k1phi, k2pt, k2eta, k2phi; 

  bool m1id, m2id, m1tmid, m2tmid, m1mvaid, m2mvaid, m1rmvaid, m2rmvaid; 
  double m1rmvabdt, m2rmvabdt; 
  double m1trigm, m2trigm; 
  int m1gt, m2gt, k1gt, k2gt; 
  int m1pix, m1bpix, m1bpixl1, m1pv; 
  int m2pix, m2bpix, m2bpixl1, m2pv; 

  double mudist, mudeltar; 
  double g1pt, g2pt, g3pt, g4pt, g1eta, g2eta, g3eta, g4eta, gtau;
  int    g1id, g2id; 
  double t1pt, t1eta, t2pt, t2eta, t3pt, t3eta, t4pt, t4eta; 

  double mpsi, mkk;
  double psipt, psieta, psiphi, phipt, phieta, phiphi, dr;
  double md0, dm, ptd0; 

  double hm1pt, hm1eta, hm1phi, hm2pt, hm2eta, hm2phi; 

  double pvdchi2, m1iso, m1xpdist, m2iso, m2xpdist, othervtx;

};

#endif
