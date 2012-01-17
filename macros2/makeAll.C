void makeAll(std::string files = "anaBmm.scratch.files", int channels = 31) {

  if (channels & 16) {
    plotEfficiencies a3(files.c_str()); 
    a3.makeAll(7); 
  }

  if (channels & 8) {
    plotPU a2(files.c_str()); 
    a2.makeAll(); 
  }

  if (channels & 4) {
    plotOverlays a1(files.c_str()); 
    a1.makeAll(); 
  }

  if (channels & 1) {
    plotResults a0(files.c_str()); 
    a0.fNormProcessed = false; 
    a0.fDoUseBDT = false; 
    a0.fDoApplyCowboyVeto = false;   
    a0.fDoApplyCowboyVetoAlsoInSignal = false;   
    a0.computeNormUL();
    a0.computeCsBF();
    a0.acceptancePerProcess();
    
    if (channels & 2) {
      a0.fNormProcessed = false; 
      a0.fDoUseBDT = true; 
      a0.fDoApplyCowboyVeto = false;   
      a0.fDoApplyCowboyVetoAlsoInSignal = false;   
      a0.computeNormUL();
      a0.computeCsBF();
      //    acceptancePerProcess();
    }
  }
  
}
