#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>

//int main(int argc, char *argv[])
int main()
{
//    if (argc < 3){
//     std::cout<<"Syntax is [Command] [Same_Event_root_file] [Mixed_Event_root_file]"<<std::endl;
//     exit(EXIT_FAILURE);
//   }

  //TFile *corr = TFile::Open((TString)argv[1]);
  TFile *corr = TFile::Open("fout_frixione.root");

  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
  TH1F* h_dPhi_iso_same[nztbins];    
  TH1F* h_dPhi_iso_mix[nztbins];

  TFile *mix = TFile::Open("fout_mixed_frixione.root");
  
  if (mix == NULL) {
    std::cout << " file 2 fail" << std::endl;
    exit(EXIT_FAILURE);
  }

  for (int izt = 0; izt<nztbins; izt++){
    h_dPhi_iso_same[izt] = (TH1F*)corr->Get(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));

   if (h_dPhi_iso_same[izt] == NULL) {
      std::cout << "tree 1 fail" << std::endl;
      exit(EXIT_FAILURE);}

    h_dPhi_iso_mix[izt] = (TH1F*)mix->Get(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));

    if (h_dPhi_iso_mix[izt] == NULL) {
      std::cout << "tree 2 fail" << std::endl;
      exit(EXIT_FAILURE);}
    //TAxis *mix_xaxis = h_dPhi_iso_mix->GetXaxis();
    //Int_t mix_bin_phi = xaxis->FindBin(0.0);
    Int_t mix_bin_phi = h_dPhi_iso_mix[izt]->FindBin(0.0);
    //TAxis *mix_yaxis = h_dPhi_iso_mix->GetYaxis();
    //Int_t mix_bin_eta = yaxis->FindBin(0.0);
    //Double_t mix_integ = h_dPhi_iso_mix->Integral(mix_bin_phi,mix_bin_phi,mix_bin_eta,mix_bin_eta);
    Double_t mix_integ = h_dPhi_iso_mix[izt]->Integral(mix_bin_phi,mix_bin_phi);
    Double_t norm = 1/mix_integ;
    h_dPhi_iso_mix[izt]->Scale(norm);
    h_dPhi_iso_same[izt]->Divide(h_dPhi_iso_mix[izt]);
  }
  
//   TAxis *xaxis = h_dPhi_iso_same->GetXaxis();
//   Int_t bin_phi = xaxis->FindBin(1.57);
//   TAxis *yaxis = h_dPhi_iso_same->GetYaxis();
//   Int_t bin_eta_neg = yaxis->FindBin(-0.6);
//   Int_t bin_eta_pos = yaxis->FindBin(0.6);
//   std::cout << bin_eta_neg << std::endl;
//   std::cout << bin_eta_pos << std::endl;
//   corr->Print();
//   Double_t same_integ = h_dPhi_iso_same->Integral(bin_phi,bin_phi);
  
  //TFile *mix = TFile::Open((TString)argv[2]);
  //TFile *mix = TFile::Open("bin_gamma_hadron.root");
  
  //TH2D  *h2 = (TH2D*)mix->Get("Correlation"); 
  //corr->Close();       
  //mix->Close();
  
  //TH1D *PhiProj = h_dPhi_iso_same->ProjectionX("phi",bin_eta_neg,bin_eta_pos);
  //PhiProj->Write();
  //TH1OAD *EtaProj = h_dPhi_iso_same->ProjectionY("eta",0,-1);
  //EtaProj->Write();

  TFile *MyFile = new TFile("ztbin_Ratio_Correlation.root","RECREATE");
  for (int izt = 0; izt<nztbins; izt++){
    //h_dPhi_iso_same[izt]->SetTitle("#gamma-H Correlation [GS-Mix Corrected]");
    //h_dPhi_iso_same[izt]->GetYaxis()->SetTitle("#Delta #eta");
    h_dPhi_iso_same[izt]->GetXaxis()->SetTitle("#Delta #phi");
    h_dPhi_iso_same[izt]->Write();
  }
  MyFile->Close();
}
