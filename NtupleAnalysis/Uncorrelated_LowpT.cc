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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  
  // Read configuration file for the DNN_min and DNN_max variables
  FILE* config = fopen("Corr_config.yaml", "r");
  
  // Default values of various variables used in the file (actual values are to be determined by the configuration file)
  // Cut variables
  double DNN_min = 0.55;
  double DNN_max = 0.85;
  double pT_min = 8;
  double pT_max = 16;
  double Eta_max = 0.6;
  double Cluster_min = 2;
  double EcrossoverE_min = 0.03;
  int selection_number = 3;
    
  // The bounds for the events to fal into the isolation and nonisolation areas
  double iso_max = 2.0;
  double noniso_min = 5.0;
  double noniso_max = 15.0;
  
  // Delta eta
  double deta_max = 0.6;
  
  // Zt bins
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
   
  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  //ptbins[0] = 10.0; ptbins[1] = 11.5; ptbins[2] = 13; ptbins[3] = 14.5; ptbins[4] = 16;
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  // Number of bins in correlation functions
  int n_correlationbins = 18;
    
  // Which branch should be used to determine whether a cluster should fall into iso, noniso, or neither
  //isolationDet determiner = CLUSTER_ISO_TPC_04;
  isolationDet determiner = CLUSTER_ISO_ITS_04;

  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
      if (line[0] == '#') {
          continue;
      }
      
      // Declare char arrays needed to read the line
      char key[MAX_INPUT_LENGTH];
      char dummy[MAX_INPUT_LENGTH];
      char value[MAX_INPUT_LENGTH];
      
      // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
      key[0] = '\0';
      value[0] = '\0';
      sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
      
      // Use if statements to detect, based on key, which variable the line's content should be used to fill and fill that variable
      if (strcmp(key, "DNN_min") == 0) {
          // Assign DNN_min to the double-converted version of value
          DNN_min = atof(value);
          std::cout << "DNN_min is " << DNN_min << std::endl;
      }
      else if (strcmp(key, "DNN_max") == 0) {
          // Assign DNN_max to the double-converted version of value
          DNN_max = atof(value);
          std::cout << "DNN_max is " << DNN_max << std::endl;
      }
      else if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min is " << pT_min << std::endl;
      }
      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max is " << pT_max << std::endl;
      }
      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max is " << Eta_max << std::endl;
      }
      else if (strcmp(key, "Cluster_min") == 0) {
          Cluster_min = atof(value);
          std::cout << "Cluster_min is " << Cluster_min << std::endl;
      }
      else if (strcmp(key, "EcrossoverE_min") == 0) {
          EcrossoverE_min = atof(value);
          std::cout << "EcrossoverE_min is " << EcrossoverE_min << std::endl;
      }
      else if (strcmp(key, "iso_max") == 0) {
          iso_max = atof(value);
          std::cout << "iso_max is " << iso_max << std::endl;
      }
      else if (strcmp(key, "noniso_min") == 0) {
          noniso_min = atof(value);
          std::cout << "noniso_min is " << noniso_min << std::endl;
      }
      else if (strcmp(key, "noniso_max") == 0) {
          noniso_max = atof(value);
          std::cout << "noniso_max is " << noniso_max << std::endl;
      }
      else if (strcmp(key, "deta_max") == 0) {
          deta_max = atof(value);
          std::cout << "deta_max is " << deta_max << std::endl;
      }
      else if (strcmp(key, "Correlation_func_bins") == 0) {
          n_correlationbins = atoi(value);
          std::cout << "Bins in a correlation function: " << n_correlationbins << std::endl;
      }
      else if (strcmp(key, "track_selection_num") == 0) {
          selection_number = atoi(value);
          std::cout << "Number of Selection that tracks must pass: " << selection_number << std::endl;
      }
      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins = -1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) {
                  v++;
              }
	      nztbins++;
              
              while (*v != ']' && (isdigit(*v) || *v == '.')) {
                  v++;
              }
          }
          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) {
                  v++;
              }              
              ztbins[i] = atof(v);
              i++;              
              while (*v != ']' && (isdigit(*v) || *v == '.')) {
                  v++;
              }
          }
          std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
          for (int i = 0; i <= nztbins; i++)
              std::cout << ztbins[i] << ", ";
          std::cout << "}\n";
          }

      else if (strcmp(key, "Pt_bins") == 0) {
        nptbins = -1;
        for (const char *v = value; *v != ']';) {
          while (*v != ']' && !isdigit(*v)) {
            v++;
          }
          nptbins++;

          while (*v != ']' && (isdigit(*v) || *v == '.')) {
            v++;
          }
        }
        ptbins = new float[nptbins + 1];
        int i = 0;
        for (const char *v = value; *v != ']' ;) {
          while (*v != ']' && !isdigit(*v)) {
            v++;
          }
          ptbins[i] = atof(v);
          i++;
          while (*v != ']' && (isdigit(*v) || *v == '.')) {
            v++;
          }
        }
	std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
        for (int i = 0; i <= nptbins; i++)
	  std::cout << ptbins[i] << ", ";
	std::cout << "}\n";
      }

      else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
          if (strcmp(value, "cluster_iso_tpc_04") == 0){
              determiner = CLUSTER_ISO_TPC_04;
              std::cout << "cluster_iso_tpc_04 will determine the isolation and non-isolation placement" << std::endl;
          }
          else if (strcmp(value, "cluster_iso_its_04") == 0){
              determiner = CLUSTER_ISO_ITS_04;
              std::cout << "cluster_iso_its_04 will determine the isolation and non-isolation placement" << std::endl;
          }
          else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_TPC_04_02;
              std::cout << "cluster_frixione_tpc_04_02 will determine the isolation and non-isolation placement" << std::endl;
          }
          else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_ITS_04_02;
              std::cout << "cluster_frixione_its_04_02 will determine the isolation and non-isolation placement" << std::endl;
          }
          else {
              std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
              exit(EXIT_FAILURE);
          }
      }
      else {
          std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
      }
  }
  fclose(config);
  
  for (int i = 0; i <= nztbins; i++){
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  }

  for (int i = 0; i <= nptbins; i++){
    std::cout << "pt bound: " << ptbins[i] << std::endl;
  }

  TCanvas canvas("canvas", "");
  TH1D histogram0("histogram0", "", 16, 8.0, 16.0);
  TH1D histogram3("histogram3", "", 18, -0.5,1.5);
  TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);

  TH2D* PtIsoDist = new TH2D("PtIsoDist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",24,10,16,5,-0.5,2);
  TH2D* LowDNN_PtIsoDist = new TH2D("LowDNN_PtIsoDist","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",24,10,16,5,-0.5,2);

  // Function declarations of h_dPhi_iso and h_dPhi_noniso

  TH2D* Corr[nztbins*nptbins];
  TH2D* IsoCorr[nztbins*nptbins];
  TH2D* AntiIsoCorr[nztbins*nptbins];


  TH1D* h_Iso_triggers[nptbins];
  TH1D* h_AntiIso_triggers[nptbins];
  //TH1D h_ntrig_low_pT("N_Trig_LowpT_Uncorrelated","Number of low p_{T}^{Clus.} Triggers",2,-0.5,1);

  float ntriggers_iso = 0;
  float ntriggers_Antiiso = 0;
  float n_lowpt_triggers = 0;
  
  int n_eta_bins = 14;
  int n_phi_bins = 24;
  
  TH1D* Uncorrelated_lowpT[nztbins];

    for (int izt = 0; izt<nztbins; izt++){
      Uncorrelated_lowpT[izt] = new TH1D(
      Form("Uncorrelated_LowpT__zT%1.0f_%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
      Form("Uncorrelated Low p_{T}^{Clus.} 1-2 GeV z_{T} %1.2f-%1.2f",ztbins[izt],ztbins[izt+1]), 
      n_phi_bins,-M_PI/2,3*M_PI/2);
    }

  for (int ipt = 0; ipt <nptbins; ipt++) {
    h_Iso_triggers[ipt] = new TH1D(
    Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]),
    "Number of Isolated Photon Triggers", 2, -0.5,1.0);

    h_AntiIso_triggers[ipt] = new TH1D(
    Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),
    "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

    for (int izt = 0; izt<nztbins; izt++){

      Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      10*ztbins[izt],10*ztbins[izt+1]),"#gamma-H [all] Correlation", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);

      Corr[izt+ipt*nztbins]->Sumw2();
      Corr[izt+ipt*nztbins]->SetMinimum(0.);

      IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1],
      10*ztbins[izt],10*ztbins[izt+1]),"#gamma-H [Iso] Correlation", n_phi_bins,-M_PI/2,3*M_PI/2,n_eta_bins, -1.4, 1.4);

      IsoCorr[izt+ipt*nztbins]->Sumw2();
      IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);

      AntiIsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      10*ztbins[izt],10*ztbins[izt+1]),"#gamma-H [AntiIso] Correlation", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);

      AntiIsoCorr[izt+ipt*nztbins]->Sumw2();
      AntiIsoCorr[izt+ipt*nztbins]->SetMinimum(0.);

    }//zt bins
  }//pt bins                           
  
  //histogram2.Sumw2();
  //histogram3.Sumw2();

  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    
    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }  
    //_tree_event->Print();

    //variables
    Double_t primary_vertex[3];
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
    UChar_t track_its_ncluster[NTRACK_MAX];
    Float_t track_its_chi_square[NTRACK_MAX];
    Float_t track_dca_xy[NTRACK_MAX];

    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];

    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
     
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);
  
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);
    _tree_event->SetBranchAddress("track_its_ncluster", &track_its_ncluster);
    _tree_event->SetBranchAddress("track_its_chi_square", &track_its_chi_square);
    _tree_event->SetBranchAddress("track_dca_xy", &track_dca_xy);

    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);

    Long64_t nentries = _tree_event->GetEntries();         
    std::cout << " Total Number of entries in TTree: " << nentries << std::endl;

    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
    //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
 
      for (ULong64_t n = 0; n < ncluster; n++) {
	//if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue; //select pt of photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue; //cut edges of detector
	if( not(cluster_ncell[n]>Cluster_min)) continue;   //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters

	float isolation;
	if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
	else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
	else isolation = cluster_frixione_its_04_02[n];
	if (isolation>iso_max) continue;
	
	//if(isolation>iso_max) continue;    
	if (cluster_pt[n]>pT_min and cluster_pt[n]<pT_max){
	if(isolation<iso_max){
	  if ((cluster_s_nphoton[n][1] > 1) && (cluster_s_nphoton[n][1]<2)){
	    histogram0.Fill(cluster_pt[n]); //isolated deep-photon pt spectra
	    ntriggers_iso += 1;
	    //look at pt Distribution in Isolation bins
	    for (double Iso_bin = -0.5; Iso_bin < iso_max; Iso_bin += 0.5){
	      if (isolation > Iso_bin && isolation < Iso_bin+0.5) PtIsoDist->Fill(cluster_pt[n],isolation);
	    }
	    for (int ipt = 0; ipt < nptbins; ipt++){
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
	      h_Iso_triggers[ipt]->Fill(0);
	      }
	    }
	  }
	}
	//if(isolation>noniso_min && isolation<noniso_max){
	
	if(isolation<iso_max){
	  if ((cluster_s_nphoton[n][1]>0.0) && (cluster_s_nphoton[n][1]<0.3)){
	    h_ntrig.Fill(0.5);
	    ntriggers_Antiiso += 1;
	    //look at pt Distrobution in Isolation bins
	    for (double Iso_bin = -0.5; Iso_bin < iso_max; Iso_bin += 0.5){
	      if (isolation > Iso_bin && isolation < Iso_bin+0.5) LowDNN_PtIsoDist->Fill(cluster_pt[n],isolation);
	    }
	    for (int ipt = 0; ipt < nptbins; ipt++){
              if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
		h_AntiIso_triggers[ipt]->Fill(0);
	      } 
	    }  
	  }
	}
	}//cluster pT if
	if (cluster_pt[n]>1 && cluster_pt[n]<2){
	  if (isolation<iso_max) n_lowpt_triggers +=1;
	}
	const int TrackCutBit =16;
	for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {            
	  if(track_pt[itrack] < 0.5) continue; //1GeV Tracks
	  if(track_pt[itrack] > 30) continue;
	  if((track_quality[itrack]&TrackCutBit)==0) continue; //select only tracks that pass selection 3
	  if(abs(track_eta[itrack]) > 0.8) continue;
	  if( not(track_its_ncluster[itrack]>4)) continue;
	  if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	  if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;

	  Double_t zt = track_pt[itrack]/cluster_pt[n];
	  Float_t deta =  cluster_eta[n]-track_eta[itrack];
	  Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack])/TMath::Pi();
	  //if(!(TMath::Abs(deta)<deta_max)) continue; // delta eta cut
	  if(dphi<-0.5) dphi +=2;
        
	  Float_t DeltaPhi = cluster_phi[n] - track_phi[itrack];
	  if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi                                                       
	  if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
	  Float_t DeltaEta = cluster_eta[n] - track_eta[itrack];
	  //if ((TMath::Abs(DeltaPhi) < 0.1) && (TMath::Abs(DeltaEta) < 0.1)) continue;
	  if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue;
	  //if (sqrt(DeltaPhi*DeltaPhi + DeltaEta*DeltaEta) < 0.4) continue;

	  if (cluster_pt[n]>pT_min and cluster_pt[n]<pT_max){
	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
	      for(int izt = 0; izt<nztbins ; izt++){
		if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
		  if(isolation< iso_max){
		    if (cluster_s_nphoton[n][1] > DNN_min && cluster_s_nphoton[n][1] < DNN_max){
		      IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		    }
		  }
		  if(isolation<iso_max){
		    if (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < 0.3){ //sel deep photons                                       
		      AntiIsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		    }
		  }
		  Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		}//if in zt bin                                                                                               
	      } // end loop over zt bins                                                                         
	    }//end if in pt bin                                                                                                  
	  }//end pt loop bin   
	  }//if pT cluster
	  if (cluster_pt[n]>1 && cluster_pt[n]<2 && isolation<iso_max){
	    if (TMath::Abs(DeltaEta) <= 0.6) //this is ok because careful binning makes eventual subtraction ok
	      for (int izt = 0; izt < nztbins; izt++){
		if (zt>ztbins[izt] and zt<ztbins[izt+1]){
		  Uncorrelated_lowpT[izt]->Fill(DeltaPhi);
		}
	      }
	  }       
	}//end loop over tracks
      }//end loop on clusters. 
    } //end loop over events  
  }//end loop over samples

    // Write to fout

    //TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");
    TFile* fout = new TFile("fout_lowDNN.root","RECREATE");
    histogram0.Write("DeepPhotonSpectra");
    h_ntrig.Write("ntriggers");
    std::cout<<"Clusters Passed Iosalation: "<<ntriggers_iso<<std::endl;

    PtIsoDist->Write();
    LowDNN_PtIsoDist->Write();
    for (int izt = 0; izt < nztbins; izt++){
      Uncorrelated_lowpT[izt]->Scale(1.0/1.2);//scale by eta
      Uncorrelated_lowpT[izt]->Scale(1.0/n_lowpt_triggers); //per trigger yield
      Uncorrelated_lowpT[izt]->Write();
    }
    std::cout<<n_lowpt_triggers<<std::endl;

    for (int ipt = 0; ipt<nptbins; ipt++){
      h_Iso_triggers[ipt]->Write();
    }
    for (int ipt = 0; ipt<nptbins; ipt++){
      h_AntiIso_triggers[ipt]->Write();
    }
    
    for (int ipt = 0; ipt<nptbins; ipt++){
      for (int izt = 0; izt<nztbins; izt++){
        Corr[izt+ipt*nztbins]->Write();
      }
      for (int izt = 0; izt<nztbins; izt++){
        IsoCorr[izt+ipt*nztbins]->Write();
      }
      for (int izt = 0; izt<nztbins; izt++){
        AntiIsoCorr[izt+ipt*nztbins]->Write();
      }
    }
    fout->Close();     

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
