/**
   This program produces energy response plots from Monte-Carlo simulations
*/

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
#include <TGraphAsymmErrors.h>

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];
    
  dummyv[0] = strdup("main");
  TApplication application("", &dummyc, dummyv);    
  std::cout <<" Number of arguments " << argc << std::endl; 
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);
        
    if (file == NULL) {
      std::cout << " fail; could not open file" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
        
    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = NULL;
    std::cout << " About to try getting the ttree" << std::endl;
    _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout << "First try did not got trying again" << std::endl;
      _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail; could not find _tree_event " << std::endl;
	exit(EXIT_FAILURE);
      }
    }
        
    if (_tree_event == NULL) {
      std::cout << " fail; the _tree_event is NULL " << std::endl;
      exit(EXIT_FAILURE);
    }
        
    //
 
    const int xjbins = 20;
    const int phibins = 20;  

    TH1D h_cutflow("h_cutflow","cut flow for photons", 10, -0.5,9.5);

   
    TH1D h_reco("h_reco", "reco photons filled with pt reco", 50, 0, 50);    
    TH1D h_reco_truthpt("h_reco_truthpt", "reco photons filled with truthpt reco", 50, 0, 50); 
    TH1D h_truth("h_truth", "truth photons", 50, 0, 50);


    TH1D hBR_reco("hBR_reco", "Isolated cluster, bkg region", 30, 0, 30);
    TH1D hSR_reco("hSR_reco", "Isolated cluster, signal region", 30, 0, 30);

    TH1D hSR_Xj("hSR_Xj", "Xj distribution, Signal region", xjbins, 0.0,2.0);
    TH1D hBR_Xj("hBR_Xj", "Xj distribution, BKG region", xjbins, 0.0,2.0);

    TH1D hSR_Xj_truth("hSR_Xj_truth", "True Xj distribution, Signal region", xjbins, 0.0,2.0);
    TH1D hBR_Xj_truth("hBR_Xj_truth", "True Xj distribution, BKG region",xjbins, 0.0,2.0);

    TH2D h_Xj_Matrix("h_Xj_Matrix", "Truth/Reco matrix", xjbins, 0.0,2.0, xjbins, 0.0,2.0);
    TH1D h_Xj_truth("h_Xj_truth", "Xj truth distribution", xjbins, 0.0,2.0);    
    TH1D hSR_dPhi("hSR_dPhi", "delta phi gamma-jet signal region", phibins, 0, TMath::Pi());
    TH1D hBR_dPhi("hBR_dPhi", "delta phi gamma-jet background region", phibins, 0, TMath::Pi());
    TH1D h_dPhi_truth("h_dPhi_truth", "delta phi gamma-jet truth MC", phibins, 0, TMath::Pi());

    hSR_Xj.Sumw2();
    hBR_Xj.Sumw2();
    hSR_Xj_truth.Sumw2();
    hBR_Xj_truth.Sumw2();
    h_Xj_truth.Sumw2();     
    hSR_dPhi.Sumw2();
    hBR_dPhi.Sumw2();
    h_dPhi_truth.Sumw2();
    hSR_reco.Sumw2();
    hBR_reco.Sumw2();

    h_Xj_Matrix.Sumw2();

    hSR_Xj.SetTitle("; X_{j} ; 1/N_{#gamma} dN_{J#gamma}/dX_{j}");
    hBR_Xj.SetTitle("; X_{j} ; 1/N_{#gamma} dN_{J#gamma}/dX_{j}");   
    h_Xj_truth.SetTitle("; X_{j}^{true} ; counts");

    
    TCanvas* canvas = new TCanvas();
        
    //you define variables
    Double_t primary_vertex[3];
    Bool_t is_pileup_from_spd_5_08;
    Bool_t is_pileup_from_spd_3_08;
   
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
        
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
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];   
 
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];
    Float_t cell_e[17664];

    //Jets
    UInt_t njet_ak04its;
    Float_t jet_ak04its_pt_raw[NTRACK_MAX];
    Float_t jet_ak04its_eta_raw[NTRACK_MAX];
    Float_t jet_ak04its_phi[NTRACK_MAX];
    
    Float_t jet_ak04its_pt_truth[NTRACK_MAX];
    Float_t jet_ak04its_eta_truth[NTRACK_MAX];
    Float_t jet_ak04its_phi_truth[NTRACK_MAX];

    
    //Truth Jets
    UInt_t njet_truth_ak04;
    Float_t jet_truth_ak04_pt[NTRACK_MAX];
    Float_t jet_truth_ak04_eta[NTRACK_MAX];
    Float_t jet_truth_ak04_phi[NTRACK_MAX];
   
    //Int_t eg_ntrial;

    Float_t eg_cross_section;
    Int_t   eg_ntrial;
    
    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
        
    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    
    
        
    // Set the branch addresses of the branches in the TTrees
    //    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    _tree_event->SetBranchAddress("is_pileup_from_spd_3_08", &is_pileup_from_spd_3_08);

    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);
        
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
        
    _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);

    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
    

    //jets
    _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
    _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
    _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi);
    _tree_event->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth);
    _tree_event->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth);
    _tree_event->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth);


    //truth jets
    _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
    _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
    _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
    _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta);    
    // Loop over events

    Bool_t isRealData = true;
    _tree_event->GetEntry(1);
    if(nmc_truth>0) isRealData= false;
    Int_t N_SR = 0;
    Int_t N_BR = 0;
    Int_t N_eventpassed = 0;
    Int_t N_truth = 0;

    

    for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){
    //for(Long64_t ievent = 0; ievent < 500000 ; ievent++){
      _tree_event->GetEntry(ievent);
      //Eevent Selection: 
      if(not( TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position    
      if(is_pileup_from_spd_5_08) continue; //removes pileup     
      N_eventpassed +=1;


      double weight = 1.0;
      if(eg_ntrial>0) weight = eg_cross_section/(double)eg_ntrial;

      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
        //Photon Selection
	h_cutflow.Fill(0);
        if( not(cluster_pt[n]>10.0)) continue; //select pt of photons
        if( not(cluster_pt[n]<16.0)) continue;
        h_cutflow.Fill(1);
        if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells
        h_cutflow.Fill(2);
	if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; //removes "spiky" clusters
        h_cutflow.Fill(3);
        if( not(cluster_nlocal_maxima[n]<= 2)) continue; //require to have at most 2 local maxima.
        h_cutflow.Fill(4);
	if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue;
        h_cutflow.Fill(5);
        if( not(cluster_iso_its_04[n] < 1.0)) continue;
        h_cutflow.Fill(6);
        //if( not(cluster_lambda_square[n][0]<0.27)) continue; //single-photon selection (as opposed to merged photon).
        //
	Bool_t isTruePhoton = false;
        Float_t truth_pt = -999.0;
	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   
          if(isTruePhoton) break;
          if(index==65535) continue;
          if(mc_truth_pdg_code[index]!=22) continue;
          if(mc_truth_first_parent_pdg_code[index]!=22) continue;
          if( not (mc_truth_status[index] >0)) continue;        
          isTruePhoton = true;
          truth_pt     = mc_truth_pt[index];
	}//end loop over indices
	
	if( not(isTruePhoton or isRealData)) continue; 

	//start jet loop 

        Bool_t inSignalRegion = false;               
	if( cluster_lambda_square[n][0]<0.30){
            N_SR +=1;
            inSignalRegion = true;
	}
        else{
            N_BR +=1;
            inSignalRegion = false;
	}

	for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) { //start loop over jets
          if(not (jet_ak04its_pt_raw[ijet]>5)) continue; 
          if(not (TMath::Abs(jet_ak04its_eta_raw[ijet]<0.5))) continue;
          Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jet_ak04its_phi[ijet] - cluster_phi[n]));
	  

          if(inSignalRegion){
	    hSR_dPhi.Fill(dphi,weight);
	  }
	  else{
	    hBR_dPhi.Fill(dphi,weight);
	  }
 	  if(not (dphi>0.4)) continue; 
	  Float_t xj = jet_ak04its_pt_raw[ijet]/cluster_pt[n];
	  //std::cout <<"truthptjet: " << jet_ak04its_pt_truth[ijet] << "reco pt jet" << jet_ak04its_pt_raw[ijet] << std::endl;           
          Float_t xj_truth = jet_ak04its_pt_truth[ijet]/truth_pt; 
           
	  if( inSignalRegion){
            hSR_Xj.Fill(xj, weight);
            hSR_Xj_truth.Fill(xj_truth, weight);
            h_Xj_Matrix.Fill(xj_truth, xj, weight);
	    //std::cout<<" xj " << xj << " " << " xj_truth "<< xj_truth << std::endl;

	  }
          else{
	    hBR_Xj.Fill(xj,weight);
            hBR_Xj_truth.Fill(xj_truth,weight);
	  }



	}//end loop over jets

        if(inSignalRegion){
	  hSR_reco.Fill(cluster_pt[n],weight);
	}
        else{
          hBR_reco.Fill(cluster_pt[n],weight);
	}
	//fill in this histogram only photons that can be traced to a generated non-decay photon.	
        h_reco_truthpt.Fill(truth_pt,weight);
        h_reco.Fill(cluster_pt[n],weight); 

      }//end loop on clusters
      //loop over truth particles
     
      //std::cout << " about to loop over mc truth particles " << std::endl;
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
        if(not(mc_truth_pt[nmc]>10.0)) continue;
	if(not(mc_truth_pt[nmc]<16.0)) continue;
	if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22){
	  //std::cout << mc_truth_pt[nmc] << "phi " << mc_truth_phi[nmc] << " eta " << mc_truth_eta[nmc] << 
	  //  " code: " << mc_truth_pdg_code[nmc] << " status " << int(mc_truth_status[nmc]) << " parentpdg " << mc_truth_first_parent_pdg_code[nmc] << std::endl;    
	  h_truth.Fill(mc_truth_pt[nmc],weight);
          N_truth +=1;
	  //std::cout << " number of truth jets " << njet_truth_ak04 << std::endl;
	  for (ULong64_t ijet = 0; ijet < njet_truth_ak04; ijet++) {
            if(not(jet_truth_ak04_pt[ijet]>5.0)) continue;
	    if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<0.5)) continue;
	    Float_t dphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(jet_truth_ak04_phi[ijet] - mc_truth_phi[nmc]));
	    //std::cout<< dphi_truth << std::endl;
            
            if( not(dphi_truth>0.4)) continue;
	    h_dPhi_truth.Fill(dphi_truth,weight);
	    Float_t xj_truth = jet_truth_ak04_pt[ijet]/mc_truth_pt[nmc];
            h_Xj_truth.Fill(xj_truth,weight);
	  }//end loop over truth jets
	}
      }//end loop over mc particles

     // Create the file label, to be used within the filenames, to represent the source file
     

      h_truth.SetLineColor(2);
      THStack* hs = new THStack("hs","stack histo for plotting");
      hs->Add(&h_reco);
      hs->Add(&h_truth);

      if (ievent % 10000 == 0) {
	//SR_Xj.Draw("e1x0nostack");
        hSR_dPhi.Draw("e1x0nostack");
	std::cout << ievent << " " << _tree_event->GetEntries() << std::endl;
        canvas->Update();
      } 

     }//end over events

    std::cout << " Numbers of events passing selection " << N_eventpassed << std::endl;
    std::cout << " Number of clusters in signal region " << N_SR << std::endl;
    std::cout << " Number of clusters in background region " << N_BR << std::endl;
    std::cout << " Number of truth photons " << N_truth << std::endl;

    std::string opened_files = "";
    //  for (int iarg = 1; iarg < argc; iarg++) {
    std::string filepath = argv[iarg];    
    opened_files = "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    
    TFile* fout = new TFile(Form("GammaJet_config%s.root", opened_files.c_str()),"RECREATE");
    fout->Print();
    hSR_Xj.Scale(1.0/N_SR);
    hBR_Xj.Scale(1.0/N_BR);
    hSR_Xj_truth.Scale(1.0/N_SR);
    hBR_Xj_truth.Scale(1.0/N_BR);

    hSR_dPhi.Scale(1.0/N_SR);
    hBR_dPhi.Scale(1.0/N_BR);

    h_dPhi_truth.Scale(1.0/N_truth);
    h_Xj_truth.Scale(1.0/N_truth);

    hSR_Xj.Write("hSR_Xj");
    hBR_Xj.Write("hBR_Xj");
    h_cutflow.Write("ClusterCutFlow");
    h_Xj_truth.Write("h_Xj_truth");
    hSR_dPhi.Write("hSR_dPhi");
    hBR_dPhi.Write("hBR_dPhi");
    h_dPhi_truth.Write("h_dPhi_truth");
    hSR_reco.Write("hSR_recopt");
    hBR_reco.Write("hBR_recopt");
    h_Xj_Matrix.Write("xj_matrix");     
    hSR_Xj_truth.Write("hSR_Xj_truth");
    hBR_Xj_truth.Write("hBR_Xj_truth");

  std::cout << " ending " << std::endl;
  }//end of arguments
  return EXIT_SUCCESS;
 
}
