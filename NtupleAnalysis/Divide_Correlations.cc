#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
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
  TFile *corr = TFile::Open("fout_largeDNNBkgrnd.root");
  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

  //TFile *mix = TFile::Open("fout_mixed_frixione.root");  
  //TFile *mix = TFile::Open("fout_minbiasfd_2GeV_mix.root");
TFile *mix = TFile::Open("fout_fdc_mixed_frixione.root");
  if (mix == NULL) {
    std::cout << " file 2 fail" << std::endl;
    exit(EXIT_FAILURE);
  }

  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;

  TH2F* Corr_Iso_same[nztbins]; 
  TH2F* Corr_Iso_mix[nztbins];

  TH2F* Corr_AntiIso_same[nztbins];
  TH2F* Corr_AntiIso_mix[nztbins];

  TH1D* IsoPhiProj[nztbins];
  TH1D* AntiIsoPhiProj[nztbins];

  TH1D* IsoPhiProj_LargeEta[nztbins];
  TH1D* AntiIsoPhiProj_LargeEta[nztbins];

  //Grab Histos, correct with mixed, project, get pedestal
  for (int izt = 0; izt<nztbins; izt++){
    Corr_Iso_same[izt] = (TH2F*)corr->Get(Form("IsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
   if (Corr_Iso_same[izt] == NULL) {
      std::cout << "tree 1 fail" << std::endl;
      exit(EXIT_FAILURE);}

   //Corr_Iso_mix[izt] = (TH2F*)mix->Get(Form("IsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
   Corr_Iso_mix[izt] = (TH2F*)mix->Get(Form("Correlation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
    if (Corr_Iso_mix[izt] == NULL) {
      std::cout << "tree 2 fail" << std::endl;
      exit(EXIT_FAILURE);
    }

    Corr_AntiIso_same[izt] = (TH2F*)corr->Get(Form("AntiIsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
    if (Corr_AntiIso_same[izt] == NULL) {
      std::cout << "tree 1 fail" << std::endl;
      exit(EXIT_FAILURE);}

    //Corr_AntiIso_mix[izt] = (TH2F*)mix->Get(Form("AntiIsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
    Corr_AntiIso_mix[izt] = (TH2F*)mix->Get(Form("Correlation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]));
    if (Corr_AntiIso_mix[izt] == NULL) {
      std::cout << "tree 2 fail" << std::endl;
      exit(EXIT_FAILURE);
    }

    //Normalize Mixing to 1 at ∆eta∆phi = 0,0
    TAxis *mix_xaxis = Corr_Iso_mix[izt]->GetXaxis();
    Int_t mix_bin_phi = mix_xaxis->FindBin(0.0);
    TAxis *mix_yaxis = Corr_Iso_mix[izt]->GetYaxis();
    Int_t mix_bin_eta = mix_yaxis->FindBin(0.0);
    Double_t mix_Iso_integ = Corr_Iso_mix[izt]->Integral(mix_bin_phi,mix_bin_phi,mix_bin_eta,mix_bin_eta);
    std::cout<<"Iso Mix Integ: "<<mix_Iso_integ<<std::endl;
    Double_t Iso_norm = 1.0/mix_Iso_integ;
    Corr_Iso_mix[izt]->Scale(Iso_norm);

    Double_t mix_AntiIso_integ = Corr_AntiIso_mix[izt]->Integral(mix_bin_phi,mix_bin_phi,mix_bin_eta,mix_bin_eta);
    std::cout<<"Anti Iso Mix Integral: "<<mix_AntiIso_integ<<std::endl;
    Double_t AntiIso_norm = 1.0/mix_AntiIso_integ;
    Corr_AntiIso_mix[izt]->Scale(AntiIso_norm);

    //Divide
    //Corr_Iso_same[izt]->Divide(Corr_Iso_mix[izt]);
    //Corr_AntiIso_same[izt]->Divide(Corr_AntiIso_mix[izt]);

    //Phi Projection
    TAxis *same_Iso_yaxis = Corr_Iso_same[izt]->GetYaxis();
    Int_t bin_eta_neg = same_Iso_yaxis->FindBin(-0.6);
    Int_t bin_eta_pos = same_Iso_yaxis->FindBin(0.6);
    IsoPhiProj[izt] =Corr_Iso_same[izt]->ProjectionX(Form("Iso_Phi_Projection_ztmin%1.0f_ztmax%1.0f",
				    10*ztbins[izt],10*ztbins[izt+1]),bin_eta_neg,bin_eta_pos);

    AntiIsoPhiProj[izt] = Corr_AntiIso_same[izt]->ProjectionX(Form("Anti_Phi_Projection_ztmin%1.0f_ztmax%1.0f",
				     10*ztbins[izt],10*ztbins[izt+1]),bin_eta_neg,bin_eta_pos);

    //Phi Projection at large Eta
    Int_t Large_bin_eta_neg = same_Iso_yaxis->FindBin(-1.4);
    Int_t Large_bin_eta_pos = same_Iso_yaxis->FindBin(1.4);
    IsoPhiProj_LargeEta[izt] =Corr_Iso_same[izt]->ProjectionX(Form("Iso_Phi_LargeEta_ztmin%1.0f_ztmax%1.0f",
					    10*ztbins[izt],10*ztbins[izt+1]),Large_bin_eta_neg,Large_bin_eta_pos);

    AntiIsoPhiProj_LargeEta[izt] = Corr_AntiIso_same[izt]->ProjectionX(Form("Anti_Phi_LargeEta_ztmin%1.0f_ztmax%1.0f",
							  10*ztbins[izt],10*ztbins[izt+1]),Large_bin_eta_neg,Large_bin_eta_pos);

    Int_t Large_bin_eta_temp1 = same_Iso_yaxis->FindBin(-0.6);
    Int_t Large_bin_eta_temp2 = same_Iso_yaxis->FindBin(0.6);

    TH1D* TempIso = Corr_Iso_same[izt]->ProjectionX("temprI",Large_bin_eta_temp1,Large_bin_eta_temp2);
    IsoPhiProj_LargeEta[izt]->Add(TempIso,-1);

    TH1D* TempAntiIso = Corr_AntiIso_same[izt]->ProjectionX("temprA",Large_bin_eta_temp1,Large_bin_eta_temp2);
    AntiIsoPhiProj_LargeEta[izt]->Add(TempAntiIso,-1);

    std::cout<<"Entries: "<<Corr_Iso_same[izt]->GetEntries()<<std::endl;
    TempIso->Delete();
    TempAntiIso->Delete();

    float IsoZYAM[nztbins];
    float AntiIsoZYAM[nztbins];

  }
  
  TFile *MyFile = new TFile("Ratio_Correlation.root","RECREATE");
  TCanvas *test = new TCanvas("","test");

  for (int izt = 0; izt<nztbins; izt++){
    if(izt == 0) continue;
    Corr_Iso_same[izt]->SetTitle(Form("Iso #gamma-h: 10 < p_{T}^{Clus.} <16 GeV, %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85",
                                       ztbins[izt],ztbins[izt+1]));
    Corr_Iso_same[izt]->GetXaxis()->SetTitle("#Delta #phi");
    Corr_Iso_same[izt]->GetYaxis()->SetTitle("#Delta #eta");
    Corr_Iso_same[izt]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N}{d#phi d#eta}");
    Corr_Iso_same[izt]->GetZaxis()->SetTitleOffset(1.6);
    Corr_Iso_same[izt]->Write();

    TCanvas *Isocanvas = new TCanvas("Isocanv","Isocanvas",800,700); 
    gStyle->SetOptStat("");
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.1);
    Corr_Iso_same[izt]->GetZaxis()->SetTitleOffset(2.2);
    //Corr_Iso_same[izt]->GetZaxis()->CenterTitle();
    Corr_Iso_same[izt]->Draw("SURF2");
    Isocanvas->SaveAs(Form("pics/IsoCorr%i.png",izt));
    Isocanvas->Clear();
  }

  for (int izt = 0; izt<nztbins; izt++){
    if (izt == 0) continue;
//     Corr_AntiIso_same[izt]->SetTitle(Form("Anti Iso #gamma-h: 10 < p_{T}^{Clus.} <16 GeV, %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85",
// 				      ztbins[izt],ztbins[izt+1]));
    Corr_AntiIso_same[izt]->SetTitle(Form("~Background #gamma-h: 10 < p_{T}^{Clus.} <16 GeV, %1.1f < z_{T} < %1.1f, 0.0 < DNN < 0.3",
					  ztbins[izt],ztbins[izt+1]));
    Corr_AntiIso_same[izt]->GetXaxis()->SetTitle("#Delta #phi");
    Corr_AntiIso_same[izt]->GetYaxis()->SetTitle("#Delta #eta");
    Corr_AntiIso_same[izt]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N}{d#phi d#eta}");
    Corr_AntiIso_same[izt]->GetZaxis()->SetTitleOffset(1.6);
    Corr_AntiIso_same[izt]->Write();

    TCanvas *AntiIsocanvas = new TCanvas("AntiIsocanv","AntiIsocanvas",800,700);
    gStyle->SetOptStat("");
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.1);
    Corr_AntiIso_same[izt]->GetZaxis()->SetTitleOffset(2.2);
    //Corr_AntiIso_same[izt]->GetZaxis()->CenterTitle();
    Corr_AntiIso_same[izt]->Draw("SURF2");
    AntiIsocanvas->SaveAs(Form("pics/AntiIsoCorr%i.png",izt));
    AntiIsocanvas->Clear();
  }

  for (int izt = 0; izt<nztbins; izt++){
    IsoPhiProj[izt]->SetTitle(Form("Isolated #Delta#phi Projection  %1.0f < z_{T} < %1.0f, 0.55 < DNN < 0.85, |#Delta#eta| < 0.6",
				   10*ztbins[izt],10*ztbins[izt+1]));
    IsoPhiProj[izt]->SetMinimum(0.0);
    IsoPhiProj[izt]->Write();
  }

  for (int izt = 0; izt<nztbins; izt++){
//     AntiIsoPhiProj[izt]->SetTitle(Form("Anti Isolated #Delta#phi %1.0f < z_{T} < %1.0f, 0.55 < DNN < 0.85, |#Delta#eta| < 0.6", 
// 				       10*ztbins[izt],10*ztbins[izt+1]));
    AntiIsoPhiProj[izt]->SetTitle(Form("~Background #Delta#phi %1.0f < z_{T} < %1.0f, 0.0 < DNN < 0.3, |#Delta#eta| < 0.6", 
				       10*ztbins[izt],10*ztbins[izt+1]));
    AntiIsoPhiProj[izt]->SetMinimum(0.0);
    AntiIsoPhiProj[izt]->Write();
  }

  for (int izt = 0; izt<nztbins; izt++){
    IsoPhiProj_LargeEta[izt]->SetTitle(Form("Isolated #Delta#phi %1.0f < z_{T} < %1.0f, 0.55 < DNN < 0.85, 0.8 < |#Delta#eta| < 1.4",
					    10*ztbins[izt],10*ztbins[izt+1]));
    IsoPhiProj_LargeEta[izt]->SetMinimum(0.0);
    //IsoPhiProj_LargeEta[izt]->Write();
  }

  for (int izt = 0; izt<nztbins; izt++){
    AntiIsoPhiProj_LargeEta[izt]->SetTitle(Form("~Background #Delta#phi %1.0f < z_{T} < %1.0f, 0.0 < DNN < 0.3, 0.8 < |#Delta#eta| < 1.4",
						10*ztbins[izt],10*ztbins[izt+1]));
    AntiIsoPhiProj_LargeEta[izt]->SetMinimum(0.0);
    //AntiIsoPhiProj_LargeEta[izt]->Write();
  }

  //float scales[7] = {0.1, 0.114, 0.39, 0.13, 0.046, 0.03, 0.0188};
  //float scales[7] = {0.1, 0.114, 0.39, 0.13, 0.046, 0.035, 0.024};
  float scales[7] = {0.2, 0.228, 0.78, 0.26, 0.092, 0.07, 0.048};

  for (int izt = 0; izt<nztbins; izt++){
    if(izt == 0) continue;

    TCanvas *canvas = new TCanvas("canv","canvas",1600,800);
    gStyle->SetOptStat("");
    canvas->Divide(2);
    
    canvas->cd(1);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.0);

    IsoPhiProj_LargeEta[izt]->SetLineColor(kRed);
    IsoPhiProj[izt]->SetMarkerStyle(20);
    IsoPhiProj[izt]->SetMarkerSize(1);
    IsoPhiProj[izt]->SetMarkerColor(kBlue);
    IsoPhiProj[izt]->SetTitle(Form("Iso #gamma-h: 10 < p_{T}^{Clus.} <16 GeV,  %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85", ztbins[izt],ztbins[izt+1]));

    IsoPhiProj[izt]->SetXTitle("d#phi [rad]");
    IsoPhiProj[izt]->SetYTitle("#frac{1}{N_{trig}} #frac{d^{2}N}{d#phid#eta}");
    IsoPhiProj[izt]->GetXaxis()->CenterTitle();
    IsoPhiProj[izt]->GetYaxis()->CenterTitle();
    IsoPhiProj[izt]->GetYaxis()->SetTitleOffset(2.2);

    IsoPhiProj[izt]->Scale(1.0/1.2);
    IsoPhiProj_LargeEta[izt]->Scale(1.0/1.6);
    IsoPhiProj[izt]->GetYaxis()->SetRangeUser(0,scales[izt]);
    IsoPhiProj_LargeEta[izt]->GetYaxis()->SetRangeUser(0,scales[izt]);
    IsoPhiProj[izt]->Draw();
    //IsoPhiProj_LargeEta[izt]->Draw("same");

    TLegend *Iso_legend = new TLegend(0.75,0.75,0.95,0.85);
    Iso_legend->AddEntry(IsoPhiProj[izt],"|#Delta#eta| < 0.6","p");
    //Iso_legend->AddEntry(IsoPhiProj_LargeEta[izt],"0.6 < |#Delta#eta| < 1.4","l");


    //ZYAM
    float Iso_ZYAM_Int = (IsoPhiProj[izt]->Integral(11,13))/3;
    TLine *Iso_ZYAM = new TLine(-M_PI/2,Iso_ZYAM_Int,3*M_PI/2,Iso_ZYAM_Int);
    Iso_ZYAM->SetLineColorAlpha(kRed, 0.9);
    Iso_ZYAM->SetLineWidth(4);
    Iso_ZYAM->Draw("same");
    IsoPhiProj[izt]->Draw("same");

    //Iso_legend->AddEntry(Iso_ZYAM,"ZYAM","l");
    Iso_legend->Draw();

    canvas->cd(2);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.07);

    AntiIsoPhiProj_LargeEta[izt]->SetLineColor(kRed);
    AntiIsoPhiProj[izt]->SetMarkerStyle(20);
    AntiIsoPhiProj[izt]->SetMarkerSize(1);
    AntiIsoPhiProj[izt]->SetMarkerColor(kBlue);
//     AntiIsoPhiProj[izt]->SetTitle(Form("Anti Iso #gamma-h: 10 < p_{T}^{Clus.} <16 GeV, %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85",
// 				       ztbins[izt],ztbins[izt+1]));
    AntiIsoPhiProj[izt]->SetTitle(Form("~Background #gamma-h: 10 < p_{T}^{Clus.} <16 GeV, %1.1f < z_{T} < %1.1f, 0.0 < DNN < 0.3",
				       ztbins[izt],ztbins[izt+1]));

    AntiIsoPhiProj[izt]->SetXTitle("d#phi [rad]");
    AntiIsoPhiProj[izt]->SetYTitle("1/N_{trigger} dN/d#phi");
    AntiIsoPhiProj[izt]->SetYTitle("");
    AntiIsoPhiProj[izt]->GetXaxis()->CenterTitle();
    AntiIsoPhiProj[izt]->GetYaxis()->CenterTitle();

    AntiIsoPhiProj[izt]->Scale(1/1.2);
    AntiIsoPhiProj_LargeEta[izt]->Scale(1/1.6);
    AntiIsoPhiProj[izt]->GetYaxis()->SetRangeUser(0,scales[izt]);
    AntiIsoPhiProj_LargeEta[izt]->GetYaxis()->SetRangeUser(0,scales[izt]);
    AntiIsoPhiProj[izt]->Draw("same");
    //AntiIsoPhiProj_LargeEta[izt]->Draw("same");

    Double_t AntiIso_ZYAM_Int = (AntiIsoPhiProj[izt]->Integral(11,13))/3;
    TLine *AntiIso_ZYAM = new TLine(-M_PI/2,AntiIso_ZYAM_Int,3*M_PI/2,AntiIso_ZYAM_Int);
    AntiIso_ZYAM->SetLineColorAlpha(kRed, 0.9);
    AntiIso_ZYAM->SetLineWidth(4);
    AntiIso_ZYAM->Draw("same");
    AntiIsoPhiProj[izt]->Draw("same");

    TLegend *AntiIso_legend = new TLegend(0.68,0.75,0.88,0.85);
    AntiIso_legend->AddEntry(IsoPhiProj[izt],"|#Delta#eta| < 0.6","p");
    //AntiIso_legend->AddEntry(IsoPhiProj_LargeEta[izt],"0.6 < |#Delta#eta| < 1.4","l");
    AntiIso_legend->AddEntry(AntiIso_ZYAM, "ZYAM","l");
    AntiIso_legend->Draw();
    
    canvas->SaveAs(Form("pics/phi_proj%i.png",izt));
    canvas->Clear();
  }

  TH1D* Isoped[nztbins];
  TH1D* Antiped[nztbins];

  for (int izt = 0; izt<nztbins; izt++){
    if(izt == 0) continue;

    TCanvas *pedcanvas = new TCanvas("pedcanv","pedcanvas",1600,800);
    gStyle->SetOptStat("");
    pedcanvas->Divide(2);

    pedcanvas->cd(1);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.0);

    Double_t Iso_ZYAM_Int = (IsoPhiProj[izt]->Integral(11,13))/3;
    for (int i = 1; i <25; i++) {
      Double_t y = IsoPhiProj[izt]->GetBinContent(i);
      Double_t y_error = IsoPhiProj[izt]->GetBinError(i);
      Double_t new_y = y - Iso_ZYAM_Int;
      Double_t new_y_error = (y-Iso_ZYAM_Int)*y_error/y;
      IsoPhiProj[izt]->SetBinContent(i,new_y);
      IsoPhiProj[izt]->SetBinError(i,new_y_error);
    }

//     Isoped[izt] = new TH1D(Form("Isoped_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),"Isoped",24,-M_PI/2,3*M_PI/2);
//     float Iso_ZYAM_Int = (IsoPhiProj[izt]->Integral(11,13))/3;
//     for (int i = 1; i <25; i++) Isoped[izt]->SetBinContent(i,Iso_ZYAM_Int);
//     IsoPhiProj[izt]->Add(Isoped[izt],-1);
//     IsoPhiProj[izt]->SetTitle(Form("Iso #gamma-h: PEDESTAL SUBTRACTED 10 < p_{T}^{Clus.} < 16 GeV, %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85",
// 				   ztbins[izt],ztbins[izt+1]));
    IsoPhiProj[izt]->SetTitle(Form("Iso #gamma-h: PEDESTAL SUBTRACTED, %1.1f < z_{T} < %1.1f",ztbins[izt],ztbins[izt+1]));
    IsoPhiProj[izt]->Draw();

    pedcanvas->cd(2);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.07);

    Double_t AntiIso_ZYAM_Int = (AntiIsoPhiProj[izt]->Integral(11,14))/4;
    for (int i = 1; i <25; i++) {
      Double_t y = AntiIsoPhiProj[izt]->GetBinContent(i);
      Double_t y_error = AntiIsoPhiProj[izt]->GetBinError(i);
      Double_t new_y = y - Iso_ZYAM_Int;
      Double_t new_y_error = (y-AntiIso_ZYAM_Int)*y_error/y;
      AntiIsoPhiProj[izt]->SetBinContent(i,new_y);
      AntiIsoPhiProj[izt]->SetBinError(i,new_y_error);
    }

//     Antiped[izt] = new TH1D(Form("Antiped_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),"Antiped",24,-M_PI/2,3*M_PI/2);
//     Double_t AntiIso_ZYAM_Int = (AntiIsoPhiProj[izt]->Integral(11,13))/3;
//     for (int i = 1; i <25; i++) Antiped[izt]->SetBinContent(i,AntiIso_ZYAM_Int);
//     AntiIsoPhiProj[izt]->Add(Antiped[izt],-1);
//     AntiIsoPhiProj[izt]->SetTitle(Form("AntiIso #gamma-h: PEDESTAL SUBTRACTED 10 < p_{T}^{Clus.} < 16 GeV, %1.1f < z_{T} < %1.1f, 0.55 < DNN < 0.85",
//  			      ztbins[izt],ztbins[izt+1]));

    AntiIsoPhiProj[izt]->SetTitle(Form("~Background #gamma-h: PEDESTAL SUBTRACTED, %1.1f < z_{T} < %1.1f",ztbins[izt],ztbins[izt+1]));
    AntiIsoPhiProj[izt]->Draw();

    pedcanvas->SaveAs(Form("pics/PEDEST_phi_proj%i.png",izt));
    pedcanvas->Clear();

  }

  for (int izt = 0; izt<nztbins; izt++){
    if (izt ==0) continue;
    TCanvas *OLcanvas = new TCanvas("OLcanv","OLcanvas",1000,800);
    gStyle->SetOptStat("");
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    IsoPhiProj[izt]->GetYaxis()->SetTitleOffset(1.6);
    IsoPhiProj[izt]->Scale(77305);
    AntiIsoPhiProj[izt]->Scale(77305*0.6);
    TLegend *Overlay_legend = new TLegend(0.8,0.83,0.95,0.93);
    Overlay_legend->AddEntry(IsoPhiProj[izt], "Isolated","p");
    Overlay_legend->AddEntry(AntiIsoPhiProj[izt], "Anti Isolated","l");
    AntiIsoPhiProj[izt]->SetTitle(Form("Iso & ~Background #gamma-h: PEDESTAL SUBTRACTED, %1.1f < z_{T} < %1.1f",ztbins[izt],ztbins[izt+1]));
    AntiIsoPhiProj[izt]->SetMarkerStyle(2);
    AntiIsoPhiProj[izt]->SetMarkerColor(kRed);
    AntiIsoPhiProj[izt]->SetLineColor(kRed);
    AntiIsoPhiProj[izt]->Draw();
    IsoPhiProj[izt]->Draw("same");
    Overlay_legend->Draw("same");
    OLcanvas->SaveAs(Form("pics/Overlay_phi_proj%i.png",izt));
  }

  MyFile->Close();
}
