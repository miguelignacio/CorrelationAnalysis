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

const int MAX_INPUT_LENGTH = 200;

//int main(int argc, char *argv[])
int main()
{
//    if (argc < 3){
//     std::cout<<"Syntax is [Command] [Same_Event_root_file] [Mixed_Event_root_file]"<<std::endl;
//     exit(EXIT_FAILURE);
//   }
//    TFile *corr = TFile::Open((TString)argv[1]);

  TFile *corr = TFile::Open("Same_Event_Correlation_13defv1.root");
  if (corr == NULL) {
    std::cout << "file 1 fail" << std::endl;
    exit(EXIT_FAILURE);
    }

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


  FILE* config = fopen("Corr_config.yaml", "r");

  int n_eta_bins = 0;
  int n_phi_bins = 0;

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


  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
    if (line[0] == '#') continue;

    char key[MAX_INPUT_LENGTH];
    char dummy[MAX_INPUT_LENGTH];
    char value[MAX_INPUT_LENGTH];

    // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays       
    key[0] = '\0';
    value[0] = '\0';
    sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
    
    if (strcmp(key, "Zt_bins") == 0) {
      nztbins = -1;
      for (const char *v = value; *v != ']';) {
	while (*v != ']' && !isdigit(*v)) v++;
	nztbins++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      ztbins = new float[nztbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
	while (*v != ']' && !isdigit(*v)) v++;
	ztbins[i] = atof(v);
	i++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
      for (int i = 0; i <= nztbins; i++)
	std::cout << ztbins[i] << ", ";
      std::cout << "}\n";
    }
    
    else if (strcmp(key, "Pt_bins") == 0) {
      nptbins = -1;
      for (const char *v = value; *v != ']';) {
	while (*v != ']' && !isdigit(*v)) v++;
	nptbins++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      ptbins = new float[nptbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
	while (*v != ']' && !isdigit(*v))  v++;
	ptbins[i] = atof(v);
	i++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
      for (int i = 0; i <= nptbins; i++)
	std::cout << ptbins[i] << ", ";
      std::cout << "}\n";
    }
  }

  fclose(config);

  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;

  //Config File Done -------------//


  TH2F* Same_DNN1_Corr[nztbins*nptbins]; 
  TH2F* Same_DNN2_Corr[nztbins*nptbins];

  TH2F* Mix_DNN1_Corr[nztbins*nptbins];
  TH2F* Mix_DNN2_Corr[nztbins*nptbins];

  TH1D* N_Iso_Triggers[nptbins];
  TH1D* N_BKGD_Triggers[nptbins];

  TFile *MyFile = new TFile("Same_Mix_Ratio_v1.root","RECREATE");

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
    //Normalize MixOAing to 1 at ∆eta∆phi = 0,0
    TAxis *mix_xaxis = Mix_DNN1_Corr[izt+ipt*nztbins]->GetXaxis();
    TAxis *mix_yaxis = Mix_DNN1_Corr[izt+ipt*nztbins]->GetYaxis();

    Double_t mix_DNN1_intgrl = Mix_DNN1_Corr[izt+ipt*nztbins]->GetBinContent
      (mix_xaxis->FindBin(0.0),mix_yaxis->FindBin(0.0));
    std::cout<<"Mix Event DNN1 Value at (0,0) = "<<mix_DNN1_intgrl<<std::endl;

    Double_t mix_DNN2_intgrl = Mix_DNN2_Corr[izt+ipt*nztbins]->GetBinContent
      (mix_xaxis->FindBin(0.0),mix_yaxis->FindBin(0.0)); //binning is the same

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
  for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++) MixFile[iSkim]->Close();
  MyFile->Close();
  //  theApp.Run();
  std::cout<<"Success"<<std::endl;

  return 0;
}
    
