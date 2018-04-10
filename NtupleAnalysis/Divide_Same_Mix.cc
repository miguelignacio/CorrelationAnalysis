#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLegend.h>
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

  TFile *corr = TFile::Open("Same_Event_Correlation.root");
  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

  TFile *mix = TFile::Open("Mix_Event_Correlation.root");
  if (mix == NULL) {
    std::cout << " file 2 fail" << std::endl;
    exit(EXIT_FAILURE);
  }
  //FIXME add Corr_config.yaml support
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.05; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
  //FIXME

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  TH2F* Corr_Iso_same[nztbins*nptbins]; 
  TH2F* Mix_Corr[nztbins*nptbins];

  TH2F* Corr_AntiIso_same[nztbins*nptbins];
  
  TH1D* N_Iso_Triggers[nptbins];
  TH1D* N_BKGD_Triggers[nptbins];

  TFile *MyFile = new TFile("Same_Mix_Ratio.root","RECREATE");

  for(int ipt = 0; ipt < nptbins; ipt++){
    for (int izt = 0; izt<nztbins; izt++){

      Corr_Iso_same[izt+ipt*nztbins] = (TH2F*)corr->Get(
	  Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	  1,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

      if (Corr_Iso_same[izt+ipt*nztbins] == NULL) {
	std::cout << "tree 1 fail" << std::endl;
	exit(EXIT_FAILURE);}
      
//       Mix_Corr[izt+ipt*nztbins] = (TH2F*)mix->Get(
//       Form("Correlation_ptmin%1.0f_ptmax%1.0f_ztmin%1.0f_ztmax%1.0f",
//       ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));   

//       if (Mix_Corr[izt+ipt*nztbins] == NULL) {
// 	std::cout << "tree 2 fail" << std::endl;
// 	exit(EXIT_FAILURE);
//       }
      
      Corr_AntiIso_same[izt+ipt*nztbins] = (TH2F*)corr->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
          2,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

      if (Corr_AntiIso_same[izt+ipt*nztbins] == NULL) {
	std::cout << "tree 2 fail" << std::endl;
	exit(EXIT_FAILURE);}

    //Normalize Mixing to 1 at ∆eta∆phi = 0,0
//     TAxis *mix_xaxis = Mix_Corr[izt+ipt*nztbins]->GetXaxis();
//     Int_t mix_bin_phi = mix_xaxis->FindBin(0.0);
//     TAxis *mix_yaxis = Mix_Corr[izt+ipt*nztbins]->GetYaxis();
//     Int_t mix_bin_eta = mix_yaxis->FindBin(0.0);
//     //Double_t mix_Iso_integ = Mix_Corr[izt+ipt*nztbins]->Integral(mix_bin_phi,mix_bin_phi,mix_bin_eta,mix_bin_eta);
//     Double_t mix_Iso_integ = Mix_Corr[izt+ipt*nztbins]->GetBinContent(mix_bin_phi,mix_bin_eta,mix_bin_eta);
//     std::cout<<"Iso Mix Integ: "<<mix_Iso_integ<<std::endl;
//     Double_t Iso_norm = 1.0/mix_Iso_integ;
//     Mix_Corr[izt+ipt*nztbins]->Scale(Iso_norm);

    //DIVIDE MIXING
//      Corr_Iso_same[izt+ipt*nztbins]->Divide(Mix_Corr[izt+ipt*nztbins]);
//      Corr_AntiIso_same[izt+ipt*nztbins]->Divide(Mix_Corr[izt+ipt*nztbins]);

      Corr_Iso_same[izt+ipt*nztbins]->Write();
    }//zt bin
    N_Iso_Triggers[ipt] = (TH1D*)corr->Get(Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]));
    N_BKGD_Triggers[ipt] = (TH1D*)corr->Get(Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]));
    
    for (int izt = 0; izt < nztbins; izt++){
      Corr_AntiIso_same[izt+ipt*nztbins]->Write();
    }
    N_Iso_Triggers[ipt]->Write();
    N_BKGD_Triggers[ipt]->Write();
  }//ipt

  TCanvas *canvas = new TCanvas("canv","canv", 1200,1000);
  corr->Close();
  mix->Close();
  MyFile->Close();
  //  theApp.Run();
  std::cout<<"Success"<<std::endl;

  return 0;
}
    
