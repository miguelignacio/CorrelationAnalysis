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

//variables of hPion
const int axis_pion_Cen  = 0;
const int axis_pion_Zvtx = 1; 
const int axis_pionMass  = 2;
const int axis_pionPt    = 3;
const int axis_pionEta   = 4;
const int axis_pionPhi   = 5;


const int axis_pion_PhE_1 = 7;
const int axis_pion_PhE_2 = 8;
const int axis_pion_asymmetry = 9;
const int axis_pion_PhPt_1 = 10;
const int axis_pion_PhPt_2 = 11;
const int axis_pionTrackEta = 12;
const int axis_pionTrackPhi = 14;
const int axis_pion_PhM02_1 = 16;
const int axis_pion_PhM02_2 = 17;

//variables of hPionTrack
const int axis_corr_dphi    = 10;
const int axis_corr_deta    = 11;
const int axis_corr_dabseta = 12;
const int axis_corr_pionpT  = 3; 
const int axis_corr_mass    = 2;
const int axis_corr_trackpT = 7;
const int axis_corr_zt      = 13;
const int axis_corr_xi      = 14;


//fit function
double fitf(Double_t *x,Double_t *par) {
  double arg1 = 0;
  if (par[1]!=0) arg1 = (x[0] - 0.0)/par[1];
  double arg2 = 0;
  if (par[3]!=0) arg2 = (x[0] - 1.0)/par[3];
  
  double normfactor_1 = 1.0/(sqrt(2*TMath::Pi())*par[1]);
  double normfactor_2 = 1.0/(sqrt(2*TMath::Pi())*par[3]);
  
  double fitval = par[0]*normfactor_1*TMath::Exp(-0.5*arg1*arg1)+par[2]*normfactor_2*TMath::Exp(-0.5*arg2*arg2)+par[4];
  return fitval;
}

double GaussP0(Double_t *x, Double_t *par){
    double arg1 = 0;
    if (par[1]!=0) arg1 = (x[0] - 0.0)/par[1];
    double fitval = par[0]*TMath::Exp(-0.5*arg1*arg1) + par[2];
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
    
    auto c = new TCanvas("c","c",800,800);
    double Width   = h_PionTrack->GetAxis(axisNumber)->GetBinWidth(1);
    
    std::vector<double> Yield_near;
    std::vector<double> Yield_near_err;
    std::vector<double> Yield_far;
    std::vector<double> Yield_far_err;
    
    for(int n=0; n<bins.size()-1; n++){    //loop over bins of hadron pT/zt or xi
        
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
        //c->SaveAs(Form("PDFOUTPUT/TEMP%s_%2.2f_%2.2f.png", name, min, max));
        c->Clear();
                
        //h_dphi_deta = (TH2F*)h_PionTrack->Projection(axis_corr_dphi,  axis_corr_dabseta);
        h_dphi_deta = (TH2F*)h_PionTrack->Projection(axis_corr_dphi,  axis_corr_deta);
        h_dphi_deta->Sumw2();
        
        auto c2 = new TCanvas("c2","c2",1200,600);
        c2->Divide(2);
        c2->cd(1);
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);

        c2->cd(1);
        h_dphi_deta->Draw("SURF1 FB BB");
        h_dphi_deta->SetLineWidth(1.0);        
        myText(0.4,0.82,kBlack, "Same Events");
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        
        h_dphi_deta_Mixed = (TH2F*)h_PionTrack_Mixed->Projection(axis_corr_dphi,  axis_corr_deta);
        h_dphi_deta_Mixed->Sumw2();
        
        c->cd();    
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0.0,200.0);
        
        h_deta->Draw();
        double scale = h_deta->GetMaximum()/h_dphi_deta_Mixed->GetNbinsY();
        std::cout << "SCALE " << h_deta->GetMaximum() << " scale:" << scale << std::endl;
        c->SaveAs(Form("PDFOUTPUT/%s_Deta_%2.2f_%2.2f.png", name, min, max));
        c->Clear();
    
        h_dphi_deta_Mixed->Scale(1.0/scale);
            
        c2->cd(2);
        gPad->SetPhi(-55);
        gPad->SetTheta(+20);
        h_dphi_deta_Mixed->SetTitle("Mixed");
        h_dphi_deta_Mixed->Draw("SURF1 FB BB");
        h_dphi_deta_Mixed->SetLineWidth(1.0);
        //h_dphi_deta_Mixed->SetLineColor(kGray);
        myText(0.4,0.82,kBlack, "Mixed Events");
        h_dphi_deta_Mixed->GetXaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetYaxis()->SetNdivisions(4);
        h_dphi_deta_Mixed->GetZaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetZaxis()->SetRangeUser(0,1.0);
        h_dphi_deta_Mixed->SetMinimum(0);
        h_dphi_deta_Mixed->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");

        c2->SaveAs(Form("PDFOUTPUT/%s_2D_d%2.2f_%2.2f.png", name, min, max));
        c2->Clear();
        c2->Close();
        
        c->cd();    
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0,200);
        h_deta->Scale(1.0/h_dphi_deta_Mixed->GetNbinsY());
        h_deta->Draw();
        //c->SaveAs(Form("PDFOUTPUT/%s_DetaNorm_%2.2f_%2.2f.png", name, min, max));
        c->Clear();

        //Divide 
        c->cd();
        h_dphi_deta->Divide(h_dphi_deta_Mixed);
        h_dphi_deta->Draw("SURF1 FB BB");
        myText(0.4,0.82,kBlack, "Corrected correlation");
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->SetMinimum(0);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetLineWidth(1.0);
        
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        c->SaveAs(Form("PDFOUTPUT/%s_2DCorr_%2.2f_%2.2f.png", name, min, max));
        c->Clear();
        //Projecting dphi
        
        for(int i=1; i< h_dphi_deta->GetXaxis()->GetNbins()+1 ; i++){
        double width =  h_dphi_deta->GetXaxis()->GetBinWidth(i) ;
        std::cout << " GetBinEdge" << h_dphi_deta->GetXaxis()->GetBinLowEdge(i) << " " <<  h_dphi_deta->GetXaxis()->GetBinLowEdge(i) + width << " i " << i << std::endl;
        }
        h_dphi_near = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_near",7, 14);
        
        auto h_dphi_far1 = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_far1",1, 4);
        h_dphi_far1->Sumw2();
        h_dphi_far = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_far",17, 20);
        h_dphi_far->Sumw2();
        
        h_dphi_far->Add(h_dphi_far1);
        
        h_dphi_near->SetLineColor(kAzure-3);
        h_dphi_far->SetLineColor(kRed);
        h_dphi_near->SetMarkerColor(kAzure-3);
        h_dphi_far->SetMarkerColor(kRed);
        
        h_dphi_near->Sumw2();
        h_dphi_far->Sumw2();
        h_dphi_near->Draw();
        h_dphi_near->SetMinimum(0);
        h_dphi_near->GetYaxis()->SetNdivisions(5);
        h_dphi_far->SetLineColor(2);
        h_dphi_far->SetMarkerColor(2);
        h_dphi_far->SetMarkerStyle(24);
        h_dphi_far->SetLineStyle(kDashed);
        h_dphi_far->Draw("same");
        h_dphi_near->GetYaxis()->SetTitle("arb units");
       
        
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
        ///myText(.58,.85,kBlack, Form("%2.1f < p_{T}^{#pi^{0}} < %2.0f GeV",8.0,20.0));
        myText(.58,.85,kBlack, label);
        myText(.60,.79,kAzure-3, " |#Delta#eta| <0.6");
        myText(.60,.73,kRed,     "0.8< |#Delta#eta| <1.4");

        //Call fit function:
        TF1 *func = new TF1("fit",fitf,-0.5,1.5,6);
        func->SetLineColor(kAzure);
        func->SetLineWidth(1.0);
        func->SetParameters(100,  0.2, 100,  0.5, 3000);
        func->SetParLimits(1, 0.05, 1.0); // widths
        func->SetParLimits(3, 0.05, 1.0); // widths 
        func->SetParLimits(0, 1, 10000);  //yield
        func->SetParLimits(2, 1, 10000);  //yield
         
        h_dphi_near->Fit(func);
        func->Draw("same");
        h_dphi_near->GetYaxis()->CenterTitle(kTRUE);
        c->SaveAs(Form("PDFOUTPUT/%s_hdphi_withFit_%2.2f_%2.2f.png", name, min, max));
        
        //fill vectors with fit results
        Yield_near.push_back(func->GetParameter(0));
        Yield_near_err.push_back(func->GetParError(0));
        Yield_far.push_back(func->GetParameter(2));
        Yield_far_err.push_back(func->GetParError(2));
        ///Do |eta| projection
        
        c->Clear();
        h_deta =  (TH1F*)h_dphi_deta->ProjectionX("h_deta",1,10); // project the deta bu only in the near side
        h_deta->Draw();
        h_deta->GetXaxis()->SetNdivisions(9);
        h_deta->GetYaxis()->SetNdivisions(4);
        myText(.65,.85,kBlack, "|#Delta#phi| < #pi/2");
        for(int i=1; i< h_dphi_deta->GetYaxis()->GetNbins()+1 ; i++){
            double width =  h_dphi_deta->GetYaxis()->GetBinWidth(i) ;
            std::cout << " GetBinEdge" << h_dphi_deta->GetYaxis()->GetBinLowEdge(i) << " " <<  h_dphi_deta->GetYaxis()->GetBinLowEdge(i) + width << " i " << i << std::endl;
        }
        
        std::cout << " FIT GAUSS + p0" << std::endl;
        TF1 *f = new TF1("f",GaussP0,-1.5,1.5,3);
        f->SetLineColor(kRed);
        f->SetLineWidth(1.0);
        f->SetParameters(100, 0.5, 3000);
        f->SetParLimits(0, 0.0, 1e5);
        f->SetParLimits(1, 0.05, 2.0 );
        
        h_deta->Fit(f);
        myText(.65,.79,kRed, Form("#sigma = %2.2f", f->GetParameter(1)));
        myText(.65,.73,kBlack, label); 
        c->SaveAs(Form("PDFOUTPUT/%s_DetaCorr_%2.2f_%2.2f.png", name, min, max));
        
        //undo the cuts on track pt for next iteration
        h_PionTrack->GetAxis(axisNumber)->SetRange(0,100); //(track_ptmin, track_ptmax);
        h_PionTrack_Mixed->GetAxis(axisNumber)->SetRange(0,100);//track_ptmin, track_ptmax);
    }//end loop
    
    
    TGraphErrors* g_yieldnear = new TGraphErrors();
    TGraphErrors* g_yieldfar = new TGraphErrors(); 
    double midpoint;
    double widthbin;
    double error;
    for(int n=0; n < Yield_near.size(); n++){
        std::cout<< Form("%2.2f--%2.2f = %2.f ", bins.at(n), bins.at(n+1), Yield_near.at(n)) << std::endl;
        midpoint = (bins.at(n)+bins.at(n+1))/2.0;
        widthbin = (bins.at(n+1)-bins.at(n));
        error = widthbin/2.0;
        g_yieldnear->SetPoint(n, midpoint, Yield_near.at(n)/widthbin);
        g_yieldnear->SetPointError(n, error, Yield_near_err.at(n)/widthbin);
        
        g_yieldfar->SetPoint(n, midpoint, Yield_far.at(n)/widthbin);
        g_yieldfar->SetPointError(n, error, Yield_far_err.at(n)/widthbin);
    }
    
    g_yieldfar->SetMarkerColor(kRed);
    g_yieldfar->SetLineColor(kRed);
    g_yieldnear->SetMarkerColor(kAzure-3);
    g_yieldnear->SetLineColor(kAzure-3);
    
    
    TMultiGraph* multi = new TMultiGraph();
    multi->Add(g_yieldnear);
    multi->Add(g_yieldfar);

    multi->SetMinimum(0);
    
    TString title;
    if(name=="xi"){
         title = "#xi";
    }
    
    c->Clear();
    multi->Draw("AP");
    multi->GetXaxis()->SetNdivisions(6);
    multi->GetYaxis()->SetNdivisions(5);
    multi->GetXaxis()->SetTitle("#xi");
    multi->GetYaxis()->SetTitle("dN/d#xi");
    
    //multi->SetTitle("; " + title + " ; Yield");
    c->SaveAs("PDFOUTPUT/Yields.png");
    
    
 c->Close();
 delete c;
 
return;
}

void PlotCorrelation(){
    auto fIn = new TFile("Ntuple.root","READ");
    THnSparse *h_PionTrack=0;
    fIn->GetObject("h_PionTrack", h_PionTrack);
    THnSparse *  h_PionTrack_Mixed = 0;
    fIn->GetObject("h_PionTrack_Mixed", h_PionTrack_Mixed);

    double Width_Mass =     h_PionTrack->GetAxis(axis_corr_mass)->GetBinWidth(1);
    double Width_PT   =     h_PionTrack->GetAxis(axis_corr_pionpT)->GetBinWidth(1);
    double Width_TrackpT  = h_PionTrack->GetAxis(axis_corr_trackpT)->GetBinWidth(1);
    double Width_dabsEta  = h_PionTrack->GetAxis(axis_corr_dabseta)->GetBinWidth(1);
    double Width_dEta     = h_PionTrack->GetAxis(axis_corr_deta)->GetBinWidth(1);  
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
    c->SaveAs("PDFOUTPUT/Track_pt.png");
    
    c->Clear();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    h_Mass->Draw();
    h_Mass->SetTitle("; m_{#gamma#gamma} [GeV]; entries");
    c->SaveAs("PDFOUTPUT/Mass.png");
    
    c->Clear();
    h_Pionpt->Draw();
    h_Pionpt->SetTitle("; #pi^{0} p_{T} [GeV]; entries");
    c->SaveAs("PDFOUTPUT/PionpT.png");
    
    c->Clear(); 
    gPad->SetLogy();
    gPad->SetLogx();
    h_TrackpT_Mixed->Draw();
    c->SaveAs("PDFOUTPUT/Track_ptMixed.png");
    c->Clear();
    
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    h_Mass_Mixed->Draw();
    c->SaveAs("PDFOUTPUT/MassMixed.png");
    c->Clear();
    h_Pionpt_Mixed->Draw();
    c->SaveAs("PDFOUTPUT/PionpTMixed.png");
    c->Clear(); 
    
    h_xi->Draw();
    h_xi->SetTitle("; #xi; entries");
    c->SaveAs("PDFOUTPUT/xi.png");
    c->Clear();
    
    h_xi->Scale(1.0/h_xi->Integral());
    auto h_cumulative = h_xi->GetCumulative();
    h_cumulative->Draw();
    h_cumulative->SetTitle("; #xi; cumulative prob");
    c->SaveAs("PDFOUTPUT/xi_cumulative.png");
    
    h_zt->Draw();
    h_zt->SetTitle(";Z_{T}; entries");
    gPad->SetLogy();
    c->SaveAs("PDFOUTPUT/zt.png");
    c->Clear();
    
    gPad->SetLogy(0);
    h_zt->Scale(1.0/h_zt->Integral());
    h_cumulative = h_zt->GetCumulative();
    h_cumulative->Draw();
    h_cumulative->SetTitle("; z_{T}; cumulative prob");
    c->SaveAs("PDFOUTPUT/zt_cumulative.png");
    c->Clear();
    
    gPad->SetLogx(kTRUE);
    h_TrackpT->Scale(1.0/h_TrackpT->Integral());
    h_cumulative = h_TrackpT->GetCumulative();
    h_cumulative->Draw();
    h_cumulative->GetXaxis()->SetNdivisions(7);
    h_cumulative->SetTitle("; p_{T} hadron [GeV] cumulative prob");
    c->SaveAs("PDFOUTPUT/TrackpT_cumulative.png");
    c->Clear();
    
    
    
    
    
    std::vector<double> bins_pt;
    bins_pt.push_back(0.5);
    bins_pt.push_back(1.0);
    bins_pt.push_back(2.0);
    bins_pt.push_back(5.0);

    std::vector<double> bins_xi;
    bins_xi.push_back(0.0);
    bins_xi.push_back(1.75);
    bins_xi.push_back(2.5);
    bins_xi.push_back(3.5);
    //bins_xi.push_back(10.0);

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
    c->SaveAs("PDFOUTPUT/h_pT.png");
    c->Clear();
    h_mass->Draw();
    h_mass->SetTitle("; m_{#gamma#gamma} [GeV]; Entries");
    c->SaveAs(Form("PDFOUTPUT/h_mass_%2.f_%2.f.png",minpT,maxpT));
   
   

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
    c->SaveAs("PDFOUTPUT/h_dphi_deta.png");
    
    c->Clear();
    h_dphi_deta_Mixed->Draw("SURF1 Z FB BB");
    h_dphi_deta_Mixed->GetXaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetYaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetZaxis()->SetNdivisions(4);
    h_dphi_deta_Mixed->GetXaxis()->CenterTitle(kTRUE);
    h_dphi_deta_Mixed->GetYaxis()->CenterTitle(kTRUE);
    c->SaveAs("PDFOUTPUT/h_dphi_deta_Mixed.png");
    
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
    c->SaveAs("PDFOUTPUT/h_dphi_deta_Corrected.png");

    c->Close();  
    }


void PionMass(THnSparse* h_Pion){
       ////
        
    auto c = new TCanvas();  
    h_Pion->GetAxis(axis_pionMass)->SetTitle("m_{#gamma#gamma} [GeV]");
    h_Pion->GetAxis(axis_pionPt)->SetTitle("p^{#gamma#gamma}_{T} [GeV]");
       
    auto h_MassAll = h_Pion->Projection(axis_pionMass);
    h_MassAll->Draw();
    c->SaveAs("PDFOUTPUT/Mass_All.png");
    c->Clear();
    
    //////////////////Invariant Mass Selection////////////////////////////////////////////
    double width_Mass = h_Pion->GetAxis(axis_pionMass)->GetBinWidth(1);
    h_Pion->GetAxis(axis_pionMass)->SetRange(0.100/width_Mass, 0.200/width_Mass);
    //////////////////////////////////////////////////////////////////////////////////////
    
    auto h_pT = h_Pion->Projection(axis_pionPt);
    gPad->SetLogy(1);
    h_pT->Draw();
    c->SaveAs("PDFOUTPUT/Pionpt_All.png");
    gPad->SetLogy(0);
    c->Clear();
    
    /////////////////////// pion pT selection /////////////////////////////////////////////
    double width_Pt = h_Pion->GetAxis(axis_pionPt)->GetBinWidth(1);
    h_Pion->GetAxis(axis_pionPt)->SetRange(8.0/width_Pt,20.0/width_Pt);
    //////////////////////////////////////////////////////////////////////////////////////

    h_MassAll = h_Pion->Projection(axis_pionMass);
    h_MassAll->Draw();
    c->SaveAs("PDFOUTPUT/Mass_withpTCut.png");
    c->Clear();


    auto h_2D = h_Pion->Projection(axis_pionMass,axis_pionPt);
    h_2D->Draw("colz");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    gPad->SetLogz(1);
    auto h_proj = h_2D->ProfileX("h_proj",0,100); 
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    c->SaveAs("PDFOUTPUT/2D_MassPT.png");
    c->Clear();
    
    h_2D = h_Pion->Projection(axis_pionMass,axis_pion_Cen);
    h_2D->Draw("colz");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_massproj",0,100); 
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    c->SaveAs("PDFOUTPUT/2D_MassCen.png");
        
    h_2D = h_Pion->Projection(axis_pionMass,    axis_pion_Zvtx);
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_2D->Draw("colz");
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_proj",0,100); 
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    c->SaveAs("PDFOUTPUT/2D_MassZvtx.png");

    //Mass--asymetry
    h_Pion->GetAxis(axis_pionPt)->SetRange(12.0/width_Pt,20.0/width_Pt); //setting range of pT
    
    h_2D = h_Pion->Projection(axis_pionMass, axis_pion_asymmetry);
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_2D->Draw("colz");
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_proj",0,100); 
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    c->SaveAs("PDFOUTPUT/2D_MassAsym.png");
   
    h_Pion->GetAxis(axis_pionPt)->SetRange(8.0/width_Pt,20.0/width_Pt);


    h_2D =  h_Pion->Projection(axis_pion_asymmetry,axis_pionPt);
    h_2D->Draw("COLZ");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(5);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_proj = h_2D->ProfileX("h_proj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same"); 
    c->SaveAs("PDFOUTPUT/2D_AsymPT.png");
    c->Clear();
        
    auto h_2D_TrackEtaPhi = h_Pion->Projection( axis_pionTrackEta ,  axis_pionTrackPhi );
    h_2D_TrackEtaPhi->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("PDFOUTPUT/2D_TrackEtaPhi.png");
    c->Clear();
    
    auto h_2D_M02 = h_Pion->Projection( axis_pion_PhM02_1 ,  axis_pion_PhM02_2 );
    h_2D_M02->Draw("colz");
    gPad->SetLogz(0);
    c->SaveAs("PDFOUTPUT/2D_PhotonM02.png");
    c->Clear();
    
    auto h_2D_Pt = h_Pion->Projection( axis_pion_PhE_1   , axis_pion_PhE_2);
    
    h_2D_Pt->Draw("colz");
    h_2D_Pt->SetTitle("; #gamma_{1} Energy [GeV]; #gamma_{2} Energy [GeV]");
    h_2D_Pt->GetXaxis()->SetRangeUser(0,15.0);
    h_2D_Pt->GetYaxis()->SetRangeUser(0,15.0);
    gPad->SetLogz(0);
    c->SaveAs("PDFOUTPUT/2D_PhotonPt.png");
    c->Clear();
    
    
    auto h_1D = h_Pion->Projection(axis_pionEta);
    h_1D->Draw("PL");
    h_1D->Sumw2();
    c->SaveAs("PDFOUTPUT/PionEta.png");
    c->Clear();
    
    h_1D = h_Pion->Projection(axis_pionPhi);
    h_1D->Draw();
    h_1D->Sumw2();
    c->SaveAs("PDFOUTPUT/PionPhi.png");
    c->Clear();
    
    h_2D =  h_Pion->Projection(axis_pionEta,axis_pionPhi);
    h_2D->Draw("COLZ");
    c->SaveAs("PDFOUTPUT/2D_PionEtaPhi.png");
    c->Clear();
    
    c->Close();
    delete c;
}


void Plotting(){
    SetAtlasStyle();

    TFile* fIn = new TFile("Ntuple.root","READ");
    fIn->Print();
    fIn->ls();
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion);
    
    PlotCorrelation();
    PionMass(h_Pion);
     
    //Plot(0,50);
    //Plot(14,18);
    //Plot(18,22);
    //Plot(22,26);
    //Plot(26,30);
}