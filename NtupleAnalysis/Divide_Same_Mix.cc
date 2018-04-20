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

  //TFile *corr = TFile::Open("Same_Event_Correlation.root");
  TFile *corr = TFile::Open("Same_Event_Correlation_TPC_ISO.root");
  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

  TFile *mix_LpT = TFile::Open("Mix_Event_Correlation.root");
  if (mix_LpT == NULL) {
    std::cout << " file 2 fail" << std::endl;
    exit(EXIT_FAILURE);
  }

  //TFile *mix = TFile::Open("Mix_Event_Correlation.root");
  TFile *mix_HpT = TFile::Open("Mix_Correlation_TPCISO_13C_4GeV.root");
  if (mix_HpT == NULL) {
    std::cout << " file 3 fail" << std::endl;
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

  TH2F* Same_DNN1_Corr[nztbins*nptbins]; 
  TH2F* Same_DNN2_Corr[nztbins*nptbins];

  TH2F* Mix_DNN1_Corr[nztbins*nptbins];
  TH2F* Mix_DNN2_Corr[nztbins*nptbins];

  TH1D* N_Iso_Triggers[nptbins];
  TH1D* N_BKGD_Triggers[nptbins];

  TFile *MyFile = new TFile("Same_Mix_Ratio.root","RECREATE");

  for(int ipt = 0; ipt < nptbins; ipt++){
    for (int izt = 0; izt<nztbins; izt++){

      Same_DNN1_Corr[izt+ipt*nztbins] = (TH2F*)corr->Get(
	  Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	  1,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

      if (Same_DNN1_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << "Same 1 TH2D fail" << std::endl;
	exit(EXIT_FAILURE);}

      Same_DNN2_Corr[izt+ipt*nztbins] = (TH2F*)corr->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
          2,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

      if (Same_DNN2_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << "Same 2 TH2D fail" << std::endl;
	exit(EXIT_FAILURE);}

      //Implement due to Track pT Skimming min bias
      if (ptbins[ipt]*ztbins[izt] < 4){

	Mix_DNN1_Corr[izt+ipt*nztbins] = (TH2F*)mix_LpT->Get(
        Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
        1,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));   

	if (Mix_DNN1_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << " mix Unskimmed TH2D file 1 fail" << std::endl;
	exit(EXIT_FAILURE);
	}
      
	Mix_DNN2_Corr[izt+ipt*nztbins] = (TH2F*)mix_LpT->Get(
        Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	2,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

	if (Mix_DNN2_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << " mix Unskimmed  TH2D file 2 fail" << std::endl;
        exit(EXIT_FAILURE);}

      }

      else{
	Mix_DNN1_Corr[izt+ipt*nztbins] = (TH2F*)mix_HpT->Get(
        Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
        1,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));   

	if (Mix_DNN1_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << " mix 4GeV Track Skimmed TH2D 1 fail" << std::endl;
	exit(EXIT_FAILURE);
	}
      
	Mix_DNN2_Corr[izt+ipt*nztbins] = (TH2F*)mix_HpT->Get(
        Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	2,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

	if (Mix_DNN2_Corr[izt+ipt*nztbins] == NULL) {
	std::cout << " mix 4GeV Track Skimmed TH2D 2 fail" << std::endl;
        exit(EXIT_FAILURE);
	}
      }
      
    //Normalize MixOAing to 1 at ∆eta∆phi = 0,0
    TAxis *mix_xaxis = Mix_DNN1_Corr[izt+ipt*nztbins]->GetXaxis();
    TAxis *mix_yaxis = Mix_DNN1_Corr[izt+ipt*nztbins]->GetYaxis();
    Double_t mix_DNN1_intgrl = Mix_DNN1_Corr[izt+ipt*nztbins]->GetBinContent
      (mix_xaxis->FindBin(0.0),mix_yaxis->FindBin(0.0));
    std::cout<<"Mix Event DNN1 Value at (0,0) = "<<mix_DNN1_intgrl<<std::endl;
    Double_t mix_DNN2_intgrl = Mix_DNN2_Corr[izt+ipt*nztbins]->GetBinContent
      (mix_xaxis->FindBin(0.0),mix_yaxis->FindBin(0.0)); //binning is the same
//     Double_t DNN1_norm = 1.0/mix_DNN1_intgrl;
//     Double_t DNN2_norm = 1.0/mix_DNN2_intgrl;
    Mix_DNN1_Corr[izt+ipt*nztbins]->Scale(1.0/mix_DNN1_intgrl);
    Mix_DNN2_Corr[izt+ipt*nztbins]->Scale(1.0/mix_DNN2_intgrl);

    //DIVIDE MIXING
    Same_DNN1_Corr[izt+ipt*nztbins]->Divide(Mix_DNN1_Corr[izt+ipt*nztbins]);
    Same_DNN2_Corr[izt+ipt*nztbins]->Divide(Mix_DNN2_Corr[izt+ipt*nztbins]);

      Same_DNN1_Corr[izt+ipt*nztbins]->Write();
    }//zt bin
    N_Iso_Triggers[ipt] = (TH1D*)corr->Get(Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]));
    N_BKGD_Triggers[ipt] = (TH1D*)corr->Get(Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]));
    
    for (int izt = 0; izt < nztbins; izt++){
      Same_DNN2_Corr[izt+ipt*nztbins]->Write();
    }
    N_Iso_Triggers[ipt]->Write();
    N_BKGD_Triggers[ipt]->Write();
    std::cout<<N_Iso_Triggers[ipt]->GetEntries()<<std::endl;
  }//ipt

  TCanvas *canvas = new TCanvas("canv","canv", 1200,1000);
  corr->Close();
  mix_LpT->Close();
  mix_HpT->Close();
  MyFile->Close();
  //  theApp.Run();
  std::cout<<"Success"<<std::endl;

  return 0;
}
    
