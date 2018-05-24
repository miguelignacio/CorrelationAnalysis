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
        
    //
    TH1D h_reco("h_reco", "reco photons filled with pt reco", 50, 0, 50);    
    TH1D h_reco_truthpt("h_reco_truthpt", "reco photons filled with truthpt reco", 50, 0, 50); 
    TH1D h_truth("h_truth", "truth photons", 50, 0, 50);
    
    TApplication application("", &dummyc, dummyv);
    TCanvas* canvas = new TCanvas();
        
    //you define variables
    Double_t primary_vertex[3];
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
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
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


    
    // Loop over events
    for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){
    // for(Long64_t ievent = 0; ievent < 10000 ; ievent++){
      _tree_event->GetEntry(ievent);
            
      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
        //Photon Selection
        //if( not(cluster_pt[n]>8)) {continue;} //select pt of photons
        if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; //removes "spiky" clusters
        if( not(cluster_nlocal_maxima[n]<= 2)) continue; //require to have at most 2 local maxima.
	if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue;

        //Isolation and shower shape selection:
        if( not(cluster_iso_its_04[n] < 1.0)) continue;
        if( not(cluster_lambda_square[n][0]<0.27)) continue; //single-photon selection (as opposed to merged photon).
	// Access the corresonding mc_truth particle; skip if index is 65535, which is invalid, 
        //or the truth particle pT is less than 10, or the mc_truth_pdg_code is not 22 (it's not a photon)
	//std::cout << " cluster pt " << cluster_pt[n] << " phi" << cluster_phi[n] << " eta " << cluster_eta[n] << std::endl;
	Bool_t isTruePhoton = false;
        Float_t truth_pt = -999.0;
	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   

          if(isTruePhoton) break;
          if(index==65535) continue;
	  //std::cout<<"truth, pt: " << mc_truth_pt[index] << "phi " << mc_truth_phi[index] << " eta " << mc_truth_eta[index] << 
	  //  " code: " << mc_truth_pdg_code[index] << " status " << int(mc_truth_status[index]) << " parentpdg " << mc_truth_first_parent_pdg_code[index] << std::endl; 
          if(mc_truth_pdg_code[index]!=22) continue;
          if(mc_truth_first_parent_pdg_code[index]!=22) continue;
          if( not (mc_truth_status[index] >0)) continue;        
          isTruePhoton = true;
          truth_pt     = mc_truth_pt[index];
	}//end loop over indices
	
	if(isTruePhoton){
	  //fill in this histogram only photons that can be traced to a generated non-decay photon.	
            h_reco_truthpt.Fill(truth_pt);
            h_reco.Fill(cluster_pt[n]); 
	} 
      }//end loop on clusters
       //loop over truth particles
     
      //std::cout << " about to loop over mc truth particles " << std::endl;
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
        //if(mc_truth_pt[nmc]<5.0) continue;
	if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22){
	  //std::cout << mc_truth_pt[nmc] << "phi " << mc_truth_phi[nmc] << " eta " << mc_truth_eta[nmc] << 
	  //  " code: " << mc_truth_pdg_code[nmc] << " status " << int(mc_truth_status[nmc]) << " parentpdg " << mc_truth_first_parent_pdg_code[nmc] << std::endl;    
	  h_truth.Fill(mc_truth_pt[nmc]);
	}
      } //end loop over mc truth particles
        
      // std::cout<<" ----------- "<< std::endl;
     // Create the file label, to be used within the filenames, to represent the source file
      std::string opened_files = "";
      for (int iarg = 1; iarg < argc; iarg++) {
      std::string filepath = argv[iarg];
            
      opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
      }

      h_truth.SetLineColor(2);
      THStack* hs = new THStack("hs","stack histo for plotting");
      hs->Add(&h_reco);
      hs->Add(&h_truth);

      if (ievent % 10000 == 0) {
	hs->Draw("e1x0nostack");
	std::cout << ievent << " " << _tree_event->GetEntries() << std::endl;
        canvas->Update();
      } 

    }//end over events

    TFile* fout = new TFile("PhotonEfficiency.root","RECREATE");
    TGraphAsymmErrors* eff = new TGraphAsymmErrors(&h_reco_truthpt, &h_truth);
    eff->Write("efficiency");
    h_reco.Write("Reco");
    h_truth.Write("Generated");
    h_reco_truthpt.Write("Reco_truthpt");

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
  }
}
