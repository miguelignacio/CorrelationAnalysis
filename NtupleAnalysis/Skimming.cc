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

        TTree *_tree_event = dynamic_cast<TTree *>
            (dynamic_cast<TDirectoryFile *>
             (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << " fail " << std::endl;
	  exit(EXIT_FAILURE);
        }  


   	UInt_t ncluster;
        UInt_t cluster_nmc_truth[NTRACK_MAX];
        Float_t cluster_e[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX]; 
	Float_t cell_e[NTRACK_MAX];
        Float_t cell_eta[NTRACK_MAX];
        Float_t cell_phi[NTRACK_MAX];

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
        Float_t cluster_isHardPhoton[NTRACK_MAX];
        Float_t cluster_minMass[NTRACK_MAX];
        Int_t cluster_pi0tagged[NTRACK_MAX];
        ////
	Double_t primary_vertex[3];

        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];

        Float_t cluster_lambda_square[NTRACK_MAX][2];   
    
	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
	_tree_event->SetBranchAddress("cluster_eta", cluster_eta);
	_tree_event->SetBranchAddress("cluster_phi", cluster_phi);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_lambda_square",  cluster_lambda_square);
        _tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);

	_tree_event->SetBranchAddress("cell_e", cell_e);
	_tree_event->SetBranchAddress("cell_eta", cell_eta);
	_tree_event->SetBranchAddress("cell_phi", cell_phi);

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
        _tree_event->SetBranchStatus("*centrality*",1);
	_tree_event->SetBranchStatus("*multiplicity*",1);
        _tree_event->SetBranchStatus("run_number",1);
	_tree_event->SetBranchStatus("*jet*",0);
	_tree_event->SetBranchStatus("*muon*",0);


 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

	TFile *newfile = new TFile("small.root","recreate");
	TTree *newtree = _tree_event->CloneTree(0);
        newtree->Branch("cluster_NN1", cluster_NN1, "cluster_NN1[ncluster]/F");
        newtree->Branch("cluster_NN2", cluster_NN2, "cluster_NN2[ncluster]/F");

	newtree->Branch("cluster_Lambda", cluster_Lambda, "cluster_Lambda[ncluster]/F");
        newtree->Branch("cluster_isHardPhoton", cluster_isHardPhoton, "cluster_isHardPhoton[ncluster]/F");
        newtree->Branch("cluster_minMass", cluster_minMass, "cluster_minMass[ncluster]/F");
        newtree->Branch("cluster_pi0tagged",cluster_pi0tagged,"cluster_pi0tagged[ncluster]/I");
       
        const Long64_t nevents = _tree_event->GetEntries();
        //const Long64_t nevents = 3000;  
	for(Long64_t ievent = 0; ievent < nevents ; ievent++){     
	  _tree_event->GetEntry(ievent);
          

	  
          bool KeepEvent = false;

	  for (ULong64_t n = 0; n < ncluster; n++) {
            Float_t IsTrueHardPhoton = -666.0;
	    cluster_NN1[n] = cluster_s_nphoton[n][1]; 
            cluster_NN2[n] = cluster_s_nphoton[n][2];
            cluster_Lambda[n] = cluster_lambda_square[n][0];
	    // std::cout << " cluster truth " << cluster_nmc_truth[n] << std::endl;
            //Looping over truth-particles matched to  cluster
	    for(UInt_t counter = 0 ; counter<cluster_nmc_truth[n]; counter++)
	      {
                unsigned short index = cluster_mc_truth_index[n][counter];
                if(index!=65535){
		  //if(mc_truth_first_parent_pdg_code[index]==22) continue; 
		  //std::cout << "reco pt" << cluster_pt[n] << " pdgcode:" << mc_truth_pdg_code[index] << "  parent pdg:" << mc_truth_first_parent_pdg_code[index] << " mc_truth_pt "<< mc_truth_pt[index] << "status" << 
		  //  (Int_t)mc_truth_status[index] << std::endl;

		  if(mc_truth_pdg_code[index]==22 && mc_truth_first_parent_pdg_code[index]==22) IsTrueHardPhoton = 1.0;
                }
	      }//end loop over indices
	    //if(IsTrueHardPhoton!=1.0) std::cout << "is Hard Photon? " << IsTrueHardPhoton << std::endl;
            cluster_isHardPhoton[n] = IsTrueHardPhoton;


            //pi0 veto: cluster_minMass
            float minmass = 999;
            Int_t pi0tag = 0;

            TLorentzVector v1;
	    v1.SetPtEtaPhiM(cluster_pt[n], cluster_eta[n], cluster_phi[n], 0.0);


	    for (ULong64_t m = 0; m < ncluster; m++) {
	      if(cluster_pt[m]<0.7) continue;
	      if(m==n) continue;


	      TLorentzVector v2;
	      v2.SetPtEtaPhiM(cluster_pt[m], cluster_eta[m],   cluster_phi[m], 0.0);
	      TLorentzVector pi0 = v1+v2;
              if(pi0.M()<minmass) minmass = pi0.M();
              if(pi0.M()<0.160 and pi0.M()>0.100) pi0tag=1;
	    }
	    cluster_pi0tagged[n] =pi0tag; 
	    cluster_minMass[n] = minmass;

            //event selection
	    if(cluster_pt[n]>8.0){
	      KeepEvent = true;
	      //std::cout << cluster_pt[n] << " " << cluster_NN1[n] << " " << cluster_NN2[n] << " " << cluster_Lambda[n] << std::endl;
              break;
	    }



	  }//end looop clusters 
	  if (ievent % 100000 == 0) { std::cout << " event " << ievent << std::endl;}
	  if(ievent==0){
	    newtree->Fill();
	    //std::cout <<" cell e, eta , phi " << cell_e[0] << " " << cell_eta[0] << std::endl;
	  }
	  if(TMath::Abs(primary_vertex[2])>10) continue;
	  if(KeepEvent){
            newtree->Fill();
          }
	}

	//newtree->Print();
	std::cout << "#events passing skimming: " << newtree->GetEntries() << std::endl;
	newtree->AutoSave();
	delete file;
	delete newfile;
	  
    }//end loop over samples

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
