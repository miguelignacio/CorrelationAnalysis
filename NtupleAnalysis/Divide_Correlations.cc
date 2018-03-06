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
  //TFile *corr = TFile::Open("Same_gamma_hadron.root");
  TFile *corr = TFile::Open("Correlation_def.root");

  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

  //TH2D * h1 = (TH2D*)corr->Get("Correlation");
  TH2D * h1 = (TH2D*)corr->Get("Iso_Correlation");
    if (h1 == NULL) {
      std::cout << "tree 1 fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    //h1->Sumw2(); 

    //    std::cout<<h1->Integral();
    TAxis *xaxis = h1->GetXaxis();
    Int_t bin_phi = xaxis->FindBin(1.57);
    TAxis *yaxis = h1->GetYaxis();
    Int_t bin_eta_neg = yaxis->FindBin(-0.6);
    Int_t bin_eta_pos = yaxis->FindBin(0.6);
    //integrate along eta at phi=pi/2
    corr->Print();
    Double_t same_integ = h1->Integral(bin_phi,bin_phi);

    //TFile *mix = TFile::Open((TString)argv[2]);
    //TFile *mix = TFile::Open("bin_gamma_hadron.root");
    TFile *mix = TFile::Open("GS_gamma_hadron.root");
    
    if (mix == NULL) {
      std::cout << " file 2 fail" << std::endl;
      exit(EXIT_FAILURE);
    }

    //TH2D  *h2 = (TH2D*)mix->Get("Correlation"); 
    TH2D  *h2 = (TH2D*)mix->Get("Iso_Correlation");
    if (h2 == NULL) {
      std::cout << "tree 2 fail " << std::endl;
      exit(EXIT_FAILURE);
    }
    TAxis *mix_xaxis = h1->GetXaxis();
    //Int_t mix_bin_phi = xaxis->FindBin(1.57);
    Int_t mix_bin_phi = xaxis->FindBin(0.0);
    TAxis *mix_yaxis = h1->GetYaxis();
    Int_t mix_bin_eta = yaxis->FindBin(0.0);
    mix->Print();
    //Double_t mix_integ = h2->Integral(mix_bin_phi,mix_bin_phi);
    Double_t mix_integ = h2->Integral(mix_bin_phi,mix_bin_phi,mix_bin_eta,mix_bin_eta);
    Double_t norm = 1/mix_integ;

    std::cout<<norm<<std::endl;

    h2->Scale(norm);     
    //h2->Scale(1/17.07);
    h1->Divide(h2);
    h1->Draw("SURF2");

    //corr->Close();       
    //mix->Close();


    TFile *MyFile = new TFile("Ratio_Correlation.root","RECREATE");
    h1->SetTitle("#gamma-H Correlation [GS-Mix Corrected]");
    h1->GetYaxis()->SetTitle("#Delta #eta");
    h1->GetXaxis()->SetTitle("#Delta #phi");
    h1->Write();
    TH1D *PhiProj = h1->ProjectionX("phi",bin_eta_neg,bin_eta_pos);
    PhiProj->Write();
    TH1D *EtaProj = h1->ProjectionY("eta",0,-1);
    EtaProj->Write();
}
