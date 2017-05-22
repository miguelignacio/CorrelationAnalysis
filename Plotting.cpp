#include "TList.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include <iostream>

#include "/root/atlasstyle-00-03-05/AtlasStyle.h"
#include "/root/atlasstyle-00-03-05/AtlasStyle.C"
#include "/root/atlasstyle-00-03-05/AtlasUtils.h"
#include "/root/atlasstyle-00-03-05/AtlasUtils.C"
#include "/root/atlasstyle-00-03-05/AtlasLabels.h"
#include "/root/atlasstyle-00-03-05/AtlasLabels.C"

int axis_pionCen = 0;
int axis_pionZvtx = 1; 
int axis_pionMass = 2;
int axis_pionPt = 3;

void Plot(double minpT, double maxpT){
    auto fIn = new TFile("Ntuple.root","READ");
    THnSparse *h_PionTrack=0;
    fIn->GetObject("h_PionTrack",h_PionTrack);
    THnSparse *  h_PionTrack_Mixed = 0;
    fIn->GetObject("h_PionTrack_Mixed",h_PionTrack_Mixed);
      
    ///// pion--track correlation plots
    int axis_corr_dphi    = 10;
    int axis_corr_deta    = 11;
    int axis_corr_dabseta = 12;
    int axis_corr_pionpT  = 3; 
    int axis_corr_mass    = 2;

    double binwidthMass = h_PionTrack->GetAxis(axis_corr_mass)->GetBinWidth(1);
    double Width_PT  = h_PionTrack->GetAxis(axis_corr_pionpT)->GetBinWidth(1);
    
    //Pi0 mass
    h_PionTrack->GetAxis(axis_corr_mass)->SetRange(0.100/binwidthMass, 0.500/binwidthMass);
    h_PionTrack_Mixed->GetAxis(axis_corr_mass)->SetRange(0.100/binwidthMass, 0.500/binwidthMass);
    //Pi0 pT
    h_PionTrack->GetAxis(axis_corr_pionpT)->SetRange(minpT/Width_PT, maxpT/Width_PT);
    h_PionTrack_Mixed->GetAxis(axis_corr_pionpT)->SetRange(minpT/Width_PT, maxpT/Width_PT);

    TH1F* h_pT   = (TH1F*)h_PionTrack->Projection(axis_pionPt);
    TH1F* h_mass = (TH1F*)h_PionTrack->Projection(axis_pionMass);
    
    h_pT->Sumw2();
    h_mass->Sumw2();
    
    TCanvas* c = new TCanvas();
    h_pT->Draw();
    c->SaveAs("h_pT.pdf");
    c->Clear();
    h_mass->Draw();
    h_mass->SetTitle("; m_{#gamma#gamma} [GeV]; Entries");
    c->SaveAs(Form("h_mass_%2.f_%2.f.pdf",minpT,maxpT));
   
   

    TH2F* h_dphi_deta = (TH2F*)h_PionTrack->Projection(axis_corr_dphi, axis_corr_deta);
    TH2F* h_dphi_deta_Mixed = (TH2F*)h_PionTrack_Mixed->Projection(axis_corr_dphi, axis_corr_deta);
  
    std::cout << " maximum " << h_dphi_deta->GetMaximum() << std::endl;
    c->Clear();
    //gPad->SetLogz();
    gPad->SetPhi(-60);
    gPad->SetTheta(45);
    //h_dphi_deta->SetMinimum(2);
    //h_dphi_deta->SetMaximum(h_dphi_deta->GetMaximum()/2.0);
    h_dphi_deta->GetXaxis()->SetNdivisions(4);
    h_dphi_deta->GetYaxis()->SetNdivisions(4);
    h_dphi_deta->GetZaxis()->SetNdivisions(4);
    h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
    h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);    
    h_dphi_deta->Draw("SURF1 Z FB BB");
    c->SaveAs("h_dphi_deta.pdf");
    
    c->Clear();
    h_dphi_deta_Mixed->Draw("SURF1 Z FB BB");
    h_dphi_deta_Mixed->GetXaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetYaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetZaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetXaxis()->CenterTitle(kTRUE);
    h_dphi_deta_Mixed->GetYaxis()->CenterTitle(kTRUE);
    c->SaveAs("h_dphi_deta_Mixed.pdf");
    
    h_dphi_deta->Divide(h_dphi_deta_Mixed);
    c->Clear();
    //h_dphi_deta->SetMaximum(h_dphi_deta->GetMaximum()/3.0);
    h_dphi_deta->Draw("SURF1 Z FB BB");
    h_dphi_deta->GetXaxis()->SetNdivisions(4);
    h_dphi_deta->GetYaxis()->SetNdivisions(4);
    h_dphi_deta->GetZaxis()->SetNdivisions(4);
    h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
    h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
    h_dphi_deta->SetTitle("; #Delta#eta; #Delta#phi/#pi");
    c->SaveAs("h_dphi_deta_Corrected.png");
    c->SaveAs("h_dphi_deta_Corrected.bmp");
    c->Close();  
    }


void PionMass(THnSparse* h_Pion){
    
    h_Pion->Sumw2();
    
    h_Pion->GetAxis(axis_pionMass)->SetTitle("m_{#gamma#gamma} [GeV]");
    
    double width_Mass = h_Pion->GetAxis(axis_pionMass)->GetBinWidth(1);
    h_Pion->GetAxis(axis_pionMass)->SetRange(0.100/width_Mass, 0.25/width_Mass);
    
    auto h_pT = h_Pion->Projection(axis_pionPt);
    auto c = new TCanvas();
    gPad->SetLogy(1);
    h_pT->Draw();
    c->SaveAs("pionpt.pdf");
    gPad->SetLogy(0);
    
    
    for(int n = 1; n<10; n++){
    h_Pion->GetAxis(axis_pionPt)->SetRange(5*n,5*(n+1) );
    auto h_Mass = h_Pion->Projection(axis_pionMass);
    auto h_pTtemp = h_Pion->Projection(axis_pionPt);
    h_Mass->Draw();
    c->SaveAs(Form("mass_%i.pdf",n));
    c->Clear();
    h_pTtemp->Draw();
    c->SaveAs(Form("pt_%i.pdf",n));
    }
    
    
}


void Plotting(){
    SetAtlasStyle();
    TFile* fIn = new TFile("Ntuple.root","READ");
    fIn->Print();
    fIn->ls();

    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion);
     
    PionMass(h_Pion);
    //Pion histograms:
   
     
    //Plot(0,50);
    //Plot(14,18);
    //Plot(18,22);
    //Plot(22,26);
    //Plot(26,30);
}