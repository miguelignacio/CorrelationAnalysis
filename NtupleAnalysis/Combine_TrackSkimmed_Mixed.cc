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
//    TFile *corr = TFile::Open((TString)argv[1]);

  int nTrackSkims = 3;
  float* trackPtSkims;
  trackPtSkims = new float[nTrackSkims+1];
  trackPtSkims[0] = 0.0; trackPtSkims[1] = 4.0; trackPtSkims[2] = 6.0; trackPtSkims[3] = 30.0;
  //trackPtSkims[2] = 30.0;

  TFile* MixFile[nTrackSkims];
  for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++){
    MixFile[iSkim] = TFile::Open(Form("Mix_Correlation_%1.0fGeVTracks_13defv1_ALL.root",trackPtSkims[iSkim]));
    std::cout<<Form("Mix_Correlation_%1.0fGeVTracks_13defv1_ALL.root",trackPtSkims[iSkim])<<std::endl;
    if (MixFile[iSkim] == NULL) {
      std::cout<<MixFile[iSkim]<<std::endl;
      std::cout << " Mixed file "<<iSkim<<" failed" << std::endl;
      exit(EXIT_FAILURE);}
  }
  //FIXME add Corr_config.yaml support
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.05; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
  //FIXME: Why was it written like this? why not zt[] = {a,b,c,...}

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  TH2F* Mix_DNN1_Corr[nztbins*nptbins];
  TH2F* Mix_DNN2_Corr[nztbins*nptbins];

  TFile *MyFile = new TFile("Combined_TrackSkim_Mix_Correlation_v1.root","RECREATE");

  for(int ipt = 0; ipt < nptbins; ipt++){
    for (int izt = 0; izt<nztbins; izt++){

      //Implement due to Track pT Skimming min bias
      for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++){
	
	if ((ptbins[ipt]*ztbins[izt] >= trackPtSkims[iSkim]) && ((ptbins[ipt]*ztbins[izt]) <= trackPtSkims[iSkim+1]-0.5)){
	  //a bit hacky, since ptbin*ztbin do not always correspond to integers...

	  fprintf(stderr,"pt min max: %1.2f %1.2f, zt min max %1.2f %1.2f\n",ptbins[ipt],ptbins[ipt+1],ztbins[izt],ztbins[izt+1]);
	  fprintf(stderr,"Track pT between %1.2f and %1.2f\n",ptbins[ipt]*ztbins[izt],ptbins[ipt+1]*ztbins[izt+1]);
	  
	  Mix_DNN1_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
          1,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));   

	  if (Mix_DNN1_Corr[izt+ipt*nztbins] == NULL) {
	    std::cout << " mix TH2D DNN1 read fail "<<"Mix file "<<iSkim<< std::endl;
	    exit(EXIT_FAILURE);}
      
	  Mix_DNN2_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	  2,ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));

	  if (Mix_DNN2_Corr[izt+ipt*nztbins] == NULL) {
	    std::cout << " mix TH2D DNN2 read fail "<<"Mix file "<<iSkim<<std::endl;
	    exit(EXIT_FAILURE);}     
	}

      }
      
      Mix_DNN1_Corr[izt+ipt*nztbins]->Write();
    }//zt bin
    
    for (int izt = 0; izt < nztbins; izt++){
      Mix_DNN2_Corr[izt+ipt*nztbins]->Write();
    } //second zt loop for structured writing of rootfile

  }//ipt

  TCanvas *canvas = new TCanvas("canv","canv", 1200,1000);
  for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++) MixFile[iSkim]->Close();
  MyFile->Close();
  std::cout<<"Success"<<std::endl;

  return 0;
}
    
