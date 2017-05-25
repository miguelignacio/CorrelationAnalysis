#include "TList.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TVirtualFitter.h"
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

const int axis_corr_dphi    = 10;
const int axis_corr_deta    = 11;
const int axis_corr_dabseta = 12;
const int axis_corr_pionpT  = 3; 
const int axis_corr_mass    = 2;
const int axis_corr_trackpT = 7;
const int axis_corr_zt      = 13;
const int axis_corr_xi      = 14;

   // define a function with 3 parameters
double fitf(Double_t *x,Double_t *par) {
  Double_t arg1 = 0;
  if (par[1]!=0) arg1 = (x[0] - 0.0)/par[1];
  Double_t arg2 = 0;
  if (par[3]!=0) arg2 = (x[0] - 1.0)/par[3];
  
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg1*arg1)+par[2]*TMath::Exp(-0.5*arg2*arg2)+par[4];
  return fitval;
}


void Looping(THnSparse* h_PionTrack, THnSparse* h_PionTrack_Mixed, const int axisNumber, std::vector<double> bins, char* name)
{
    
    TH2F* h_dphi_deta = 0x0;
    TH2F* h_dphi_deta_Mixed = 0x0;
    TH1F* h_dphi_near = 0x0;
    TH1F* h_dphi_far = 0x0;
    TH1F* h_trackpt = 0x0;    
    TH1F* h_deta = 0x0;
    TH1F* h_deta_Mixed = 0x0;
    
    auto c = new TCanvas();

    //double bins[3] = {0.5, 2.0, 10.0}; //bins of pt hadron.
    
    double Width   = h_PionTrack->GetAxis(axisNumber)->GetBinWidth(1);
    
    
    
    for(int n=0; n<bins.size()-1; n++){    //loop over bins of hadron pT
        //Track pT cut
        double min = bins.at(n);
        double max = bins.at(n+1);
        
        std::cout << "Min " << min << " " << max << std::endl;
        std::cout << " Width" << Width << std::endl;
        ///////////////////////////////selecting track
        h_PionTrack->GetAxis(axisNumber)->SetRange(min/Width +1, max/Width);
        h_PionTrack_Mixed->GetAxis(axisNumber)->SetRange(min/Width+1, max/Width);
        ///////////////////////////////
        c->cd();
        h_trackpt = (TH1F*)h_PionTrack->Projection(axisNumber);
        h_trackpt->Draw();
        c->SaveAs(Form("PDFOUTPUT/%s_%2.2f_%2.2f.png", name, min, max));
        c->Clear();
                
    
        h_dphi_deta = (TH2F*)h_PionTrack->Projection(axis_corr_dphi,  axis_corr_dabseta);
        h_dphi_deta->Sumw2();
        
        auto c2 = new TCanvas("c","c",1200,600);
        c2->Divide(2);
        c2->cd(1);
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);



        c2->cd(1);
        h_dphi_deta->Draw("SURF1 FB BB");
        
        
        myText(0.4,0.92,kBlack, "Same Events");
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        
        
                
        h_dphi_deta_Mixed = (TH2F*)h_PionTrack_Mixed->Projection(axis_corr_dphi,  axis_corr_dabseta);
        h_dphi_deta_Mixed->Sumw2();
        
        c->cd();    
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0,200);
        //h_deta->Sumw2();
        h_deta->Draw();
        double scale = h_deta->GetMaximum()/h_dphi_deta_Mixed->GetNbinsY();
        std::cout << "SCALE " << h_deta->GetMaximum() << " scale:" << scale << std::endl;
        c->SaveAs(Form("PDFOUTPUT/%s_Deta_%2.2f_%2.2f.pdf", name, min, max));
        c->Clear();
        
        h_dphi_deta_Mixed->Scale(1.0/scale);
        
        
        c2->cd(2);
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);
        h_dphi_deta_Mixed->SetTitle("Mixed");
        h_dphi_deta_Mixed->Draw("SURF1 FB BB");
        myText(0.4,0.92,kBlack, "Mixed Events");
        h_dphi_deta_Mixed->GetXaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetYaxis()->SetNdivisions(4);
        h_dphi_deta_Mixed->GetZaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetZaxis()->SetRangeUser(0,1.0);
        h_dphi_deta_Mixed->SetMinimum(0);
        h_dphi_deta_Mixed->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        //c2->cd(2);
        c2->SaveAs(Form("PDFOUTPUT/%s_2D_d%2.2f_%2.2f.png", name, min, max));
        //c2->Clear();
        c2->Clear();;
    
    
        c->cd();    
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0,200);
        h_deta->Scale(1.0/h_dphi_deta_Mixed->GetNbinsY());
        h_deta->Draw();
        c->SaveAs(Form("PDFOUTPUT/%s_DetaNorm_%2.2f_%2.2f.pdf", name, min, max));
        c->Clear();
    
    
    
    
        //Divide 
        c->cd();
        h_dphi_deta->Divide(h_dphi_deta_Mixed);

    
        h_dphi_deta->Draw("SURF1 FB BB");
        myText(0.4,0.92,kBlack, "Corrected correlation");
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->SetMinimum(0);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        c->SaveAs(Form("PDFOUTPUT/%s_2DCorr_%2.2f_%2.2f.png", name, min, max));
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
        myText(.63,.85,kBlack, Form("%2.1f < p_{T}^{#pi0} < %2.0f GeV",8.0,20.0));
        
        TString label;
        if(name=="pt"){
        label = Form(" %2.1f < p_{T}^{h} < %2.1f GeV", min, max);
        }
        else if(name=="zt"){
        label = Form(" %2.1f < Z_{T} < %2.1f", min, max);
        }
        else if(name=="xi"){
        label = Form(" %2.1f < #xi < %2.1f", min, max);
        }
        myText(.63,.77,kBlack, label);
        myText(.63,.69,kAzure-3, "|#Delta#eta| <0.6");
        myText(.63,.61,kRed, "0.8< |#Delta#eta| <1.4");
        //c->SaveAs(Form("PDFOUTPUT/%s_hdphi_%2.2f_%2.2f.png", name, min, max));
        //h_dphi_near->SetMinimum(0);
        //Call fit function:
        TF1 *func = new TF1("fit",fitf,-0.5,1.5,6);
        func->SetLineColor(kAzure);
        func->SetLineWidth(1.0);
      //set the parameters to the mean and RMS of the histogram
        func->SetParameters(100,  0.2, 100,  0.2, 30);
         
        h_dphi_near->Fit(func);
        func->Draw("same");
        
        c->SaveAs(Form("PDFOUTPUT/%s_hdphi_withFit_%2.2f_%2.2f.png", name, min, max));
        //undo the cuts on track pt for next iteration
        h_PionTrack->GetAxis(axisNumber)->SetRange(0,100); //(track_ptmin, track_ptmax);
        h_PionTrack_Mixed->GetAxis(axisNumber)->SetRange(0,100);//track_ptmin, track_ptmax);
    }//end loop

return;
}

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
    ///////////////////////////Mass cut//////////////////////////////////////////
    h_PionTrack->GetAxis(axis_corr_mass)->SetRange(0.100/Width_Mass, 0.200/Width_Mass);
    h_PionTrack_Mixed->GetAxis(axis_corr_mass)->SetRange(0.100/Width_Mass, 0.200/Width_Mass);
    ///////////////////////////Pi0 pT cut////////////////////////////////////////
    h_PionTrack->GetAxis(axis_corr_pionpT)->SetRange(8.0/Width_PT, 20.0/Width_PT);
    h_PionTrack_Mixed->GetAxis(axis_corr_pionpT)->SetRange(8.0/Width_PT, 20.0/Width_PT);
    /////////////////////////////////////////////////////////////////////////////
   
    auto h_TrackpT = h_PionTrack->Projection(axis_corr_trackpT);
    auto h_Mass =  h_PionTrack->Projection(axis_corr_mass);
    auto h_Pionpt = h_PionTrack->Projection(axis_corr_pionpT);
    
    auto h_TrackpT_Mixed = h_PionTrack_Mixed->Projection(axis_corr_trackpT);
    auto h_Mass_Mixed    =  h_PionTrack_Mixed->Projection(axis_corr_mass);
    auto h_Pionpt_Mixed = h_PionTrack_Mixed->Projection(axis_corr_pionpT);
    
    auto h_xi = h_PionTrack->Projection(axis_corr_xi);
    auto h_zt = h_PionTrack->Projection(axis_corr_zt);
    
    h_Mass->Sumw2();
    h_Pionpt->Sumw2();
    
    auto c = new TCanvas();
    
    gPad->SetLogy();
    gPad->SetLogx();
    h_TrackpT->Draw();
    h_TrackpT->SetTitle("; p_{T} [GeV]; entries");
    c->SaveAs("Track_pt.png");
    c->Clear();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    h_Mass->Draw();
    h_Mass->SetTitle("; m_{#gamma#gamma} [GeV]; entries");
    c->SaveAs("Mass.png");
    
    c->Clear();
    h_Pionpt->Draw();
    h_Pionpt->SetTitle("; #pi^{0} p_{T} [GeV]; entries");
    c->SaveAs("PionpT.png");
    
    c->Clear(); 
    
    gPad->SetLogy();
    gPad->SetLogx();
    h_TrackpT_Mixed->Draw();
    c->SaveAs("Track_ptMixed.png");
    c->Clear();
    
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    h_Mass_Mixed->Draw();
    c->SaveAs("MassMixed.png");
    c->Clear();
    h_Pionpt_Mixed->Draw();
    c->SaveAs("PionpTMixed.png");
    c->Clear(); 
    
    h_xi->Draw();
    h_xi->SetTitle("; #xi; entries");
    c->SaveAs("xi.png");
    c->Clear();
    
    h_zt->Draw();
    h_zt->SetTitle(";Z_{T}; entries");
    gPad->SetLogy();
    c->SaveAs("zt.png");
    c->Clear();
    
    std::vector<double> bins_pt;
    bins_pt.push_back(0.5);
    bins_pt.push_back(2.0);
    bins_pt.push_back(5.0);
    //double bins[3] = {0.5, 2.0, 10.0}; //bins of pt hadron.
    std::vector<double> bins_xi;
    bins_xi.push_back(0.0);
    bins_xi.push_back(2.0);
    bins_xi.push_back(4.0);
    std::vector<double> bins_zt;
    bins_zt.push_back(0.0);
    bins_zt.push_back(0.5);
    bins_zt.push_back(1.0);
    bins_zt.push_back(2.0);
    
    
    //Looping(h_PionTrack, h_PionTrack_Mixed, axis_corr_trackpT, bins_pt, "pt");
    Looping(h_PionTrack, h_PionTrack_Mixed, axis_corr_xi     , bins_xi, "xi");
    //Looping(h_PionTrack, h_PionTrack_Mixed, axis_corr_zt     , bins_zt, "zt");
    
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
    c->SaveAs("2DMassPT.png");
    c->Clear();
    
    auto h_2D_TrackEtaPhi = h_Pion->Projection( axis_pionTrackEta ,  axis_pionTrackPhi );
    h_2D_TrackEtaPhi->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DTrackEtaPhi.png");
    c->Clear();
    
    auto h_2D_M02 = h_Pion->Projection( axis_pion_PhM02_1 ,  axis_pion_PhM02_2 );
    h_2D_M02->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DM02.png");
    c->Clear();
    
    auto h_2D_Pt = h_Pion->Projection( axis_pion_PhPt_1  , axis_pion_PhPt_2);
    h_2D_Pt->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("2DPt.png");
    c->Clear();
    
    auto h_pT = h_Pion->Projection(axis_pionPt);
    gPad->SetLogy(1);
    h_pT->Draw();
    c->SaveAs("pionpt.png");
    gPad->SetLogy(0);
    

    for(int n = 1; n<1; n++){
    h_Pion->GetAxis(axis_pionPt)->SetRange(5*n,5*(n+1) );
    auto h_Mass = h_Pion->Projection(axis_pionMass);
    auto h_pTtemp = h_Pion->Projection(axis_pionPt);
    h_Mass->Draw();
    c->SaveAs(Form("mass_%i.png",n));
    c->Clear();
    h_pTtemp->Draw();
    c->SaveAs(Form("pt_%i.png",n));
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
    
    h_Pion->Print(); 
     
    PlotCorrelation();
    //PionMass(h_Pion);
    //Pion histograms:
   
     
    //Plot(0,50);
    //Plot(14,18);
    //Plot(18,22);
    //Plot(22,26);
    //Plot(26,30);
}