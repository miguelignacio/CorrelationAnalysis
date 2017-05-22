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
int axis_pionTrackEta = 12;
int axis_pionTrackPhi = 14;

int axis_pion_PhPt_1 = 10;
int axis_pion_PhPt_2 = 11;

int axis_pion_PhM02_1 = 16;
int axis_pion_PhM02_2 = 17;

int axis_corr_dphi    = 10;
int axis_corr_deta    = 11;
int axis_corr_dabseta = 12;
int axis_corr_pionpT  = 3; 
int axis_corr_mass    = 2;
int axis_corr_trackpT = 7;



void PlotCorrelation(){
    auto fIn = new TFile("Ntuple.root","READ");
    THnSparse *h_PionTrack=0;
    fIn->GetObject("h_PionTrack", h_PionTrack);
    THnSparse *  h_PionTrack_Mixed = 0;
    fIn->GetObject("h_PionTrack_Mixed", h_PionTrack_Mixed);

    double Width_Mass = h_PionTrack->GetAxis(axis_corr_mass)->GetBinWidth(1);
    double Width_PT   = h_PionTrack->GetAxis(axis_corr_pionpT)->GetBinWidth(1);
    double Width_TrackpT   = h_PionTrack->GetAxis(axis_corr_trackpT)->GetBinWidth(1);
    double Width_dabsEta = h_PionTrack->GetAxis(axis_corr_dabseta)->GetBinWidth(1);
    //Mass cut
    h_PionTrack->GetAxis(axis_corr_mass)->SetRange(0.100/Width_Mass, 0.200/Width_Mass);
    h_PionTrack_Mixed->GetAxis(axis_corr_mass)->SetRange(0.100/Width_Mass, 0.200/Width_Mass);
    //Pi0 pT cut
    h_PionTrack->GetAxis(axis_corr_pionpT)->SetRange(8.0/Width_PT, 20.0/Width_PT);
    h_PionTrack_Mixed->GetAxis(axis_corr_pionpT)->SetRange(8.0/Width_PT, 20.0/Width_PT);
    
    auto h_TrackpT = h_PionTrack->Projection(axis_corr_trackpT);
    auto h_Mass =  h_PionTrack->Projection(axis_corr_mass);
    auto h_Pionpt = h_PionTrack->Projection(axis_corr_pionpT);
    auto c = new TCanvas();
    gPad->SetLogy();
    h_TrackpT->Draw();
    c->SaveAs("Track_pt.pdf");
    c->Clear();
    gPad->SetLogy(0);
    
    h_Mass->Draw();
    c->SaveAs("Mass.pdf");
    c->Clear();
    
    h_Pionpt->Draw();
    c->SaveAs("PionpT.pdf");
    c->Clear(); 
    
    
    
    TH2F* h_dphi_deta = 0x0;
    TH2F* h_dphi_deta_Mixed = 0x0;
    TH1F* h_dphi_near = 0x0;
    TH1F* h_dphi_far = 0x0;
    TH1F* h_trackpt = 0x0;
    for(int n=2; n<5; n++){    
        //Track pT cut
        double track_ptmin = (n)*Width_TrackpT;
        double track_ptmax = (n+2)*Width_TrackpT;
        
        std::cout << "track min " << track_ptmin << " " << track_ptmax << std::endl;
        std::cout << " Track pt width" << Width_TrackpT << std::endl;
        //selecting track
        h_PionTrack->GetAxis(axis_corr_trackpT)->SetRange(n+1,n+2); //(track_ptmin, track_ptmax);
        h_PionTrack_Mixed->GetAxis(axis_corr_trackpT)->SetRange(n+1,n+2);//track_ptmin, track_ptmax);

        h_trackpt = (TH1F*)h_PionTrack->Projection(axis_corr_trackpT);
        h_trackpt->Draw();
        c->SaveAs(Form("Trackpt_%2.2f_%2.2f.pdf",track_ptmin, track_ptmax));
        
        c->Clear();
        h_dphi_deta = (TH2F*)h_PionTrack->Projection(axis_corr_dphi,  axis_corr_dabseta);
        h_dphi_deta_Mixed = (TH2F*)h_PionTrack_Mixed->Projection(axis_corr_dphi,  axis_corr_dabseta);
    
        //Divide 
        h_dphi_deta->Divide(h_dphi_deta_Mixed);
        gPad->SetPhi(-60);
        gPad->SetTheta(45);
    
        
        
        
        h_dphi_deta->Draw("SURF1 Z FB BB");
        h_dphi_deta->GetXaxis()->SetNdivisions(4);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->GetZaxis()->SetNdivisions(4);
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetTitle("; #Delta#eta; #Delta#phi/#pi");
        c->SaveAs(Form("2D_%2.2f_%2.2f.pdf",track_ptmin, track_ptmax));
        c->Clear();
        //Projecting dphi
        h_dphi_near = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_near",1, 5);
        h_dphi_far = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_far",6, 10);
    
        h_dphi_near->SetLineColor(kAzure-3);
        h_dphi_far->SetLineColor(kRed);
        h_dphi_near->SetMarkerColor(kAzure-3);
        h_dphi_far->SetMarkerColor(kRed);
        
        h_dphi_near->Sumw2();
        h_dphi_far->Sumw2();
        h_dphi_near->Draw();
        h_dphi_near->GetYaxis()->SetNdivisions(5);
        h_dphi_far->SetLineColor(2);
        h_dphi_far->SetMarkerColor(2);
        h_dphi_far->SetMarkerStyle(24);
        h_dphi_far->SetLineStyle(kDashed);
        h_dphi_far->Draw("same");
        h_dphi_near->GetYaxis()->SetTitle("arb units");
        myText(.63,.85,kBlack, Form(" %2.1f < p_{T}^{#pi0} < %2.0f GeV",8.0,20.0));
        myText(.63,.77,kBlack, Form(" %2.1f < p_{T}^{h} < %2.1f GeV", track_ptmin, track_ptmax));
        myText(.63,.69,kAzure-3, "|#Delta#eta| <0.6");
        myText(.63,.61,kRed, "0.8< |#Delta#eta| <1.4");
        c->SaveAs(Form("h_dphi_%2.2f_%2.2f.png",track_ptmin, track_ptmax));
        
        //undo the cuts on track pt for next iteration
        h_PionTrack->GetAxis(axis_corr_trackpT)->SetRange(0,100); //(track_ptmin, track_ptmax);
        h_PionTrack_Mixed->GetAxis(axis_corr_trackpT)->SetRange(0,100);//track_ptmin, track_ptmax);
    }
    
    delete h_PionTrack;
    delete h_PionTrack_Mixed;
    delete fIn;
    c->Close();
    delete c;
    return;
    
}









void Plot(double minpT, double maxpT){
    auto fIn = new TFile("Ntuple.root","READ");
    THnSparse *h_PionTrack=0;
    fIn->GetObject("h_PionTrack",h_PionTrack);
    THnSparse *  h_PionTrack_Mixed = 0;
    fIn->GetObject("h_PionTrack_Mixed",h_PionTrack_Mixed);
      


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
       ////
        
    auto c = new TCanvas();  
    h_Pion->GetAxis(axis_pionMass)->SetTitle("m_{#gamma#gamma} [GeV]");
    h_Pion->GetAxis(axis_pionPt)->SetTitle("p^{#gamma#gamma}_{T} [GeV]");
       
    auto h_MassAll = h_Pion->Projection(axis_pionMass);
    h_MassAll->Draw();
    c->SaveAs("MassAll.pdf");
    c->Clear();
       
    double width_Mass = h_Pion->GetAxis(axis_pionMass)->GetBinWidth(1);
    h_Pion->GetAxis(axis_pionMass)->SetRange(0.100/width_Mass, .20/width_Mass);
    
    double width_Pt = h_Pion->GetAxis(axis_pionPt)->GetBinWidth(1);
    h_Pion->GetAxis(axis_pionPt)->SetRange(5/width_Pt,50/width_Pt);
   

    
    
    auto h_2D_masspT = h_Pion->Projection(axis_pionPt, axis_pionMass);
    h_2D_masspT->Draw("colz");
    gPad->SetLogz(1);
    c->SaveAs("2DMassPT.pdf");
    c->Clear();
    
    auto h_2D_TrackEtaPhi = h_Pion->Projection( axis_pionTrackEta ,  axis_pionTrackPhi );
    h_2D_TrackEtaPhi->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DTrackEtaPhi.pdf");
    c->Clear();
    
    auto h_2D_M02 = h_Pion->Projection( axis_pion_PhM02_1 ,  axis_pion_PhM02_2 );
    h_2D_M02->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DM02.pdf");
    c->Clear();
    
    auto h_2D_Pt = h_Pion->Projection( axis_pion_PhPt_1  , axis_pion_PhPt_2);
    h_2D_Pt->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DPt.pdf");
    c->Clear();
    
    
   
    auto h_pT = h_Pion->Projection(axis_pionPt);
    gPad->SetLogy(1);
    h_pT->Draw();
    c->SaveAs("pionpt.pdf");
    gPad->SetLogy(0);
    
    
    
    
    for(int n = 1; n<1; n++){
    h_Pion->GetAxis(axis_pionPt)->SetRange(5*n,5*(n+1) );
    auto h_Mass = h_Pion->Projection(axis_pionMass);
    auto h_pTtemp = h_Pion->Projection(axis_pionPt);
    h_Mass->Draw();
    c->SaveAs(Form("mass_%i.pdf",n));
    c->Clear();
    h_pTtemp->Draw();
    c->SaveAs(Form("pt_%i.pdf",n));
    }
    c->Close();
    
}


void Plotting(){
    SetAtlasStyle();

    TFile* fIn = new TFile("Ntuple.root","READ");
    fIn->Print();
    fIn->ls();

    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion);
     
    PlotCorrelation();
    //PionMass(h_Pion);
    //Pion histograms:
   
     
    //Plot(0,50);
    //Plot(14,18);
    //Plot(18,22);
    //Plot(22,26);
    //Plot(26,30);
}