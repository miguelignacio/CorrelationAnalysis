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

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

float c_eta_array[17664];
float c_phi_array[17664];


unsigned int GetSuperModule(const unsigned int n){
   unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;


  return sm;
}

void cell_5_5(unsigned int n_5_5[], const unsigned int n,
              const unsigned int ld = 5)
{
  const unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
        const unsigned int nphi =
	  sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        n_5_5[0 * ld + 0] = n - 2 * nphi - 4;
        n_5_5[0 * ld + 1] = n - 2 * nphi - 2;
        n_5_5[0 * ld + 2] = n - 2 * nphi;
        n_5_5[0 * ld + 3] = n - 2 * nphi + 2;
        n_5_5[0 * ld + 4] = n - 2 * nphi + 4;
        if (n % 2 == 0) {
	  n_5_5[1 * ld + 0] = n - 3;
          n_5_5[1 * ld + 1] = n - 1;
          n_5_5[1 * ld + 2] = n + 1;
          n_5_5[1 * ld + 3] = n + 3;
          n_5_5[1 * ld + 4] = n + 5;
        }
        else {
          n_5_5[1 * ld + 0] = n - 2 * nphi - 5;
          n_5_5[1 * ld + 1] = n - 2 * nphi - 3;
          n_5_5[1 * ld + 2] = n - 2 * nphi - 1;
          n_5_5[1 * ld + 3] = n - 2 * nphi + 1;
          n_5_5[1 * ld + 4] = n - 2 * nphi + 3;
        }
        n_5_5[2 * ld + 0] = n - 4;
        n_5_5[2 * ld + 1] = n - 2;
        n_5_5[2 * ld + 2] = n;
        n_5_5[2 * ld + 3] = n + 2;
        n_5_5[2 * ld + 4] = n + 4;
        if (n % 2 == 0) {
	  n_5_5[3 * ld + 0] = n + 2 * nphi - 3;
          n_5_5[3 * ld + 1] = n + 2 * nphi - 1;
          n_5_5[3 * ld + 2] = n + 2 * nphi + 1;
          n_5_5[3 * ld + 3] = n + 2 * nphi + 3;
          n_5_5[3 * ld + 4] = n + 2 * nphi + 5;
        }
        else {
          n_5_5[3 * ld + 0] = n - 5;
          n_5_5[3 * ld + 1] = n - 3;
          n_5_5[3 * ld + 2] = n - 1;
          n_5_5[3 * ld + 3] = n + 1;
          n_5_5[3 * ld + 4] = n + 3;
        }
        n_5_5[4 * ld + 0] = n + 2 * nphi - 4;
        n_5_5[4 * ld + 1] = n + 2 * nphi - 2;
        n_5_5[4 * ld + 2] = n + 2 * nphi;
        n_5_5[4 * ld + 3] = n + 2 * nphi + 2;
        n_5_5[4 * ld + 4] = n + 2 * nphi + 4;
}


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
	std::cout<<"About to get TTree" << std::endl;

        TTree *_tree_event = NULL;
        _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
	  _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
	  if (_tree_event == NULL) {
	      std::cout << " fail " << std::endl;
	      exit(EXIT_FAILURE);
	  }
        }  
        //_tree_event->Print();


   	UInt_t ncluster;
        UInt_t cluster_nmc_truth[NTRACK_MAX];
        Float_t cluster_e[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX];
          
	UShort_t  cluster_cell_id_max[NTRACK_MAX];

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

	//new branches:
        Float_t cluster_NN1[NTRACK_MAX];
        Float_t cluster_NN2[NTRACK_MAX];
        Float_t cluster_Lambda[NTRACK_MAX];
        Float_t cluster_Lambda_angle[NTRACK_MAX];
        Float_t cluster_isHardPhoton[NTRACK_MAX];
        Float_t cluster_minMass[NTRACK_MAX];
        Float_t cluster_mindR_trackbit16[NTRACK_MAX];
        Float_t cluster_mindR_trackbit3[NTRACK_MAX];
	Float_t cluster_mindR_trackbit16_trackpt[NTRACK_MAX];
        Float_t cluster_mindR_trackbit3_trackpt[NTRACK_MAX];
        Float_t cluster_isoetaband_its[NTRACK_MAX];
	Float_t cluster_isoetaband_tpc[NTRACK_MAX];

        Int_t cluster_pi0tagged[NTRACK_MAX];
        ////
	Double_t primary_vertex[3];

        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];

        Float_t cluster_lambda_square[NTRACK_MAX][2];   
        Float_t cluster_lambda_square_angle[NTRACK_MAX][2];

	UChar_t track_quality[NTRACK_MAX];
	UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
	Float_t track_eta_emcal[NTRACK_MAX];
        Float_t track_phi_emcal[NTRACK_MAX];
	UChar_t track_its_ncluster[NTRACK_MAX];
	Float_t track_dca_xy[NTRACK_MAX];
	Float_t track_dca_z[NTRACK_MAX];
	Float_t track_its_chi_square[NTRACK_MAX];

	_tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
	_tree_event->SetBranchAddress("track_eta_emcal", track_eta_emcal);
        _tree_event->SetBranchAddress("track_phi_emcal", track_phi_emcal);
        _tree_event->SetBranchAddress("track_quality", track_quality);
	_tree_event->SetBranchAddress("track_dca_xy", track_dca_xy);
	_tree_event->SetBranchAddress("track_dca_z", track_dca_z);
	_tree_event->SetBranchAddress("track_its_chi_square", track_its_chi_square);
	_tree_event->SetBranchAddress("track_its_ncluster", track_its_ncluster);

	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
	_tree_event->SetBranchAddress("cluster_eta", cluster_eta);
	_tree_event->SetBranchAddress("cluster_phi", cluster_phi);
	_tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_lambda_square",  cluster_lambda_square);
	//_tree_event->SetBranchAddress("cluster_lambda_square_angle",  cluster_lambda_square_angle);
        _tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);
	_tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);

	_tree_event->SetBranchAddress("nmc_truth",&nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_pt",mc_truth_pt);
        _tree_event->SetBranchAddress("mc_truth_eta",mc_truth_eta);
        _tree_event->SetBranchAddress("mc_truth_phi",mc_truth_phi);
        _tree_event->SetBranchAddress("mc_truth_charge",mc_truth_charge);
        _tree_event->SetBranchAddress("mc_truth_pdg_code",mc_truth_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code", mc_truth_first_parent_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_e", mc_truth_first_parent_e);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pt", mc_truth_first_parent_pt);
        _tree_event->SetBranchAddress("mc_truth_first_parent_eta", mc_truth_first_parent_eta);
        _tree_event->SetBranchAddress("mc_truth_first_parent_phi", mc_truth_first_parent_phi);
	_tree_event->SetBranchAddress("mc_truth_status",  mc_truth_status);

        _tree_event->SetBranchStatus("*",0);
        _tree_event->SetBranchStatus("*cluster*",1);
	_tree_event->SetBranchStatus("*cell*",1);
	_tree_event->SetBranchStatus("track*",1);
        _tree_event->SetBranchStatus("*mc*", 1);
        _tree_event->SetBranchStatus("*vertex*",1);
        _tree_event->SetBranchStatus("*trigger*",1);
        _tree_event->SetBranchStatus("*period*",1);
        _tree_event->SetBranchStatus("*bunch*",1);
        _tree_event->SetBranchStatus("*number*",1);
        _tree_event->SetBranchStatus("*time*",1);
        _tree_event->SetBranchStatus("*pileup*",1);
	_tree_event->SetBranchStatus("*event*",1);

        _tree_event->SetBranchStatus("*centrality*",1);
	_tree_event->SetBranchStatus("*multiplicity*",1);
        _tree_event->SetBranchStatus("run_number",1);
        //_tree_event->SetBranchStatus("*Mix*",1);
	_tree_event->SetBranchStatus("*jet*",1);
	_tree_event->SetBranchStatus("*muon*",0);


 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

	TFile *newfile = new TFile("Small.root","recreate");
	TTree *newtree = _tree_event->CloneTree(0);

        //const Long64_t nevents = _tree_event->GetEntries();
	const Long64_t nevents = 1000;  
	for(Long64_t ievent = 0; ievent < nevents ; ievent++){     
	  _tree_event->GetEntry(ievent);
	  fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent+1, nevents);
          bool KeepEvent = false;

	  for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
	    if((track_quality[itrack]&16)==0) continue; //select only tracks that pass selection 16
	    if( not(track_pt[itrack]>4.0))    continue;
	    if( not(track_pt[itrack]<30.0))   continue;
	    if(abs(track_eta[itrack]) > 0.8)  continue;
	    if( not(track_its_ncluster[itrack]>4)) continue;
	    if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	    if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;
	    //cluster for loop for electron veto
	    Float_t dR = 1.0;
	    double dRmin = 0.1;
	    for (ULong64_t n = 0; n < ncluster; n++) {
	      Float_t deta =  cluster_eta[n]-track_eta_emcal[itrack];
	      if (std::isnan(track_phi_emcal[itrack])|| std::isnan(cluster_phi[n])) continue;
	      Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi_emcal[itrack])/TMath::Pi();
	      dR = TMath::Sqrt(deta*deta + dphi*dphi);
	      if (dR < dRmin) break;
	    }
	    if (dR < dRmin) continue;
	    KeepEvent = true; //passes cuts and electron veto

	  }//End track loop
	  
	  if(ievent==0){
	    newtree->Fill();
	    //std::cout <<" cell e, eta , phi " << cell_e[0] << " " << cell_eta[0] << std::endl;
	  }
	  //if(TMath::Abs(primary_vertex[2])>10) continue;
	  if(KeepEvent){
            newtree->Fill();
          }
	}
	//newtree->Print();
	std::cout << std::endl<<"#events passing skimming: " << newtree->GetEntries() << std::endl;
	newtree->AutoSave();
	delete file;
	delete newfile;
	  
    }//end loop over samples

    return EXIT_SUCCESS;
}
