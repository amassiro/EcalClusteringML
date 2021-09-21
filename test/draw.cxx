

void draw() {

  TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);
  TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);
  
  TCanvas* ccEB_PU = new TCanvas ("ccEB_PU","",1600,600);
  TCanvas* ccEE_PU = new TCanvas ("ccEE_PU","",1600,600);
  
  
  TTree* tree = (TTree*) _file0->Get("SimDigiTreeProducer/tree");
  
  
  UInt_t          run;
  UShort_t        lumi;
  UShort_t        bx;
  UInt_t          event;
  Float_t         digi_ped_subtracted_EB[612000];
  Float_t         simenergy_EB[61200*3];
  Int_t           ieta[61200];
  Int_t           iphi[61200];
  Float_t         digi_ped_subtracted_EE[146480];
  Float_t         simenergy_EE[14648*3];
  Int_t           ix[14648];
  Int_t           iy[14648];
  Int_t           iz[14648];
  
  
//   *5 + bunch_cross + 3
  
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("lumi", &lumi);
  tree->SetBranchAddress("bx", &bx);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("digi_ped_subtracted_EB", digi_ped_subtracted_EB);
  tree->SetBranchAddress("simenergy_EB", simenergy_EB);
  tree->SetBranchAddress("ieta", ieta);
  tree->SetBranchAddress("iphi", iphi);
  tree->SetBranchAddress("digi_ped_subtracted_EE", digi_ped_subtracted_EE);
  tree->SetBranchAddress("simenergy_EE", simenergy_EE);
  tree->SetBranchAddress("ix", ix);
  tree->SetBranchAddress("iy", iy);
  tree->SetBranchAddress("iz", iz);


  int MAXEVENTS = 10;
  
  TH2F* histoEB_SimEnergy = new TH2F ("histoEB_SimEnergy", "SimEnergy" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
  TH2F* histoEE_SimEnergy = new TH2F ("histoEE_SimEnergy", "SimEnergy" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
  
  histoEB_SimEnergy->GetXaxis()->SetTitle("i#phi");
  histoEB_SimEnergy->GetYaxis()->SetTitle("i#eta");
  
  histoEE_SimEnergy->GetXaxis()->SetTitle("x");
  histoEE_SimEnergy->GetYaxis()->SetTitle("y");
  
  

  TH2F* histoEB_SimEnergy_PU = new TH2F ("histoEB_SimEnergy_PU", "SimEnergy_PU" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
  TH2F* histoEE_SimEnergy_PU = new TH2F ("histoEE_SimEnergy_PU", "SimEnergy_PU" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
  
  histoEB_SimEnergy_PU->GetXaxis()->SetTitle("i#phi");
  histoEB_SimEnergy_PU->GetYaxis()->SetTitle("i#eta");
  
  histoEE_SimEnergy_PU->GetXaxis()->SetTitle("x");
  histoEE_SimEnergy_PU->GetYaxis()->SetTitle("y");


  
  
  for (int ievent = 0; ievent<MAXEVENTS; ievent++) {
    tree->GetEntry(ievent);
    
    for (int iEBchannel = 0; iEBchannel<61200; iEBchannel++) {
      if (simenergy_EB[iEBchannel*3+0] > 1) {
        histoEB_SimEnergy->Fill( iphi[iEBchannel], ieta[iEBchannel], simenergy_EB[iEBchannel*3+0] );   
      }
      if (simenergy_EB[iEBchannel*3+1] > 1) {
        histoEB_SimEnergy_PU->Fill( iphi[iEBchannel], ieta[iEBchannel], simenergy_EB[iEBchannel*3+1] );   
      }
    }
    
    for (int iEEchannel = 0; iEEchannel<14648; iEEchannel++) {
      if (simenergy_EE[iEEchannel*3+0] > 1) {
        histoEE_SimEnergy->Fill(ix[iEEchannel] + 100*(iz[iEEchannel]>0), iy[iEEchannel] , simenergy_EE[iEEchannel*3+0] );
      }
      if (simenergy_EE[iEEchannel*3+1] > 1) {
        histoEE_SimEnergy_PU->Fill(ix[iEEchannel] + 100*(iz[iEEchannel]>0), iy[iEEchannel] , simenergy_EE[iEEchannel*3+1] );
      }
    }
    
  }
  
  
  
  
  ccEE->cd();
  histoEE_SimEnergy->Draw("colz");
  
  
  ccEB->cd();
  histoEB_SimEnergy->Draw("colz");
  
  
  
  
  
  
  ccEE_PU->cd();
  histoEE_SimEnergy_PU->Draw("colz");
  
  
  ccEB_PU->cd();
  histoEB_SimEnergy_PU->Draw("colz");
  
  


  
}


