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


      std::cout << "Opening: " << (TString)argv[1] << std::endl;
        TFile *file = TFile::Open((TString)argv[1]);

	std::string filepath = argv[1];
	filepath = filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
	std::cout << filepath << std::endl;
      
	Float_t ptmin = 12.0;
        if(argc>2){
	  ptmin = strtol(argv[2],NULL,0);
        }
        Int_t nevents = 0;
        if(argc>3){
            nevents = strtol(argv[3],NULL,0);
	}

        std::cout << " Number of events requested " << nevents << std::endl; 
	std::cout << " Minimum pt for skimming " << ptmin << std::endl;
	
        if (file == NULL) {
	  std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
	std::cout<<"About to get TTree" << std::endl;

        TTree *_tree_event = NULL;
	_tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
	  _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
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

	Float_t cluster_iso_its_04[NTRACK_MAX];
    
	Float_t cluster_iso_its_04_ue[NTRACK_MAX];
            
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
        Float_t cluster_iso_its_new[NTRACK_MAX];

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
	UChar_t track_its_ncluster[NTRACK_MAX];
	Float_t track_dca_xy[NTRACK_MAX];
	Float_t track_dca_z[NTRACK_MAX];
	Float_t track_its_chi_square[NTRACK_MAX];

	Float_t ue_estimate_its_const;
        _tree_event->SetBranchAddress("ue_estimate_its_const",&ue_estimate_its_const);
	_tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
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
	_tree_event->SetBranchAddress("cluster_lambda_square_angle",  cluster_lambda_square_angle);
        _tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);
	_tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
	_tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
	_tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);



	Float_t cell_e[17664];
        Float_t cell_eta[17664];
        Float_t cell_phi[17664];

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
        _tree_event->SetBranchStatus("*trigger*",1);
        _tree_event->SetBranchStatus("*period*",1);
        _tree_event->SetBranchStatus("*bunch*",1);
        _tree_event->SetBranchStatus("*number*",1);
        _tree_event->SetBranchStatus("*time*",1);
        _tree_event->SetBranchStatus("*pileup*",1);
	_tree_event->SetBranchStatus("*event*",1);
	_tree_event->SetBranchStatus("*eg*",1);

        _tree_event->SetBranchStatus("*centrality*",1);
	_tree_event->SetBranchStatus("*multiplicity*",1);
        _tree_event->SetBranchStatus("run_number",1);
        _tree_event->SetBranchStatus("*Mix*",1);
	_tree_event->SetBranchStatus("*jet*",1);
        _tree_event->SetBranchStatus("*ue*", 1);
	_tree_event->SetBranchStatus("*muon*",0);
        //_tree_event->SetBranchStatus("*tpc*",0);


 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

	if( not(nevents>0)){
	  nevents = _tree_event->GetEntries();
        }
	TFile *newfile = new TFile(Form("Skimmed_%s_ptmin%2.1f_Nevent_%5.0i.root",filepath.c_str(),ptmin,nevents),"recreate");
	TTree *newtree = _tree_event->CloneTree(0);
        newtree->Branch("cluster_NN1", cluster_NN1, "cluster_NN1[ncluster]/F");
        newtree->Branch("cluster_NN2", cluster_NN2, "cluster_NN2[ncluster]/F");

	newtree->Branch("cluster_Lambda", cluster_Lambda, "cluster_Lambda[ncluster]/F");
        newtree->Branch("cluster_Lambda_angle", cluster_Lambda_angle, "cluster_Lambda_angle[ncluster]/F");

        newtree->Branch("cluster_isHardPhoton", cluster_isHardPhoton, "cluster_isHardPhoton[ncluster]/F");
        
        newtree->Branch("cluster_mindR_trackbit16",cluster_mindR_trackbit16,"cluster_mindR_trackbit16[ncluster]/F");
	
	newtree->Branch("cluster_mindR_trackbit16_trackpt",cluster_mindR_trackbit16_trackpt,"cluster_mindR_trackbit16_trackpt[ncluster]/F");
        
        newtree->Branch("cluster_isoetaband_its", cluster_isoetaband_its, "cluster_isoetaband_its[ncluster]/F");

        newtree->Branch("cluster_iso_its_new", cluster_iso_its_new, "cluster_iso_its_new[ncluster]/F");
 
	//newtree->Branch("cluster_isoetaband_tpc", cluster_isoetaband_tpc, "cluster_isoetaband_tpc[ncluster]/F");

        newtree->Branch("cluster_pi0tagged",cluster_pi0tagged,"cluster_pi0tagged[ncluster]/I");
     

	Float_t cluster_b5x5[NTRACK_MAX];
	Float_t cluster_b5x5_lin[NTRACK_MAX];
	unsigned int cluster_SuperModule[NTRACK_MAX];
	//newtree->Branch("cluster_b5x5", cluster_b5x5, "cluster_b5x5[ncluster]/F");
        //newtree->Branch("cluster_b5x5_lin", cluster_b5x5_lin, "cluster_b5x5_lin[ncluster]/F");
        newtree->Branch("cluster_SuperModule", cluster_SuperModule, "cluster_SuperModule[ncluster]/i");
       
	_tree_event->GetEntry(0);
        for(int i=0; i<17664; i++){
          c_eta_array[i] = cell_eta[i];
          c_phi_array[i] = cell_phi[i]; 
        }
        
       	//const Long64_t nevents = 100000;  
	for(Long64_t ievent = 0; ievent < nevents ; ievent++){     
	  _tree_event->GetEntry(ievent);
          
          bool KeepEvent = true;
          if(ptmin>0.0) KeepEvent = false;

	  if(TMath::Abs(primary_vertex[2])>10) continue;
          if(primary_vertex[2]==0.00) continue;

	  for (ULong64_t n = 0; n < ncluster; n++) {
            Float_t IsTrueHardPhoton = -666.0;
	    cluster_NN1[n] = cluster_s_nphoton[n][1]; 
            cluster_NN2[n] = cluster_s_nphoton[n][2];
            cluster_Lambda[n] = cluster_lambda_square[n][0];
	    cluster_Lambda_angle[n] = cluster_lambda_square_angle[n][0];

	    //std::cout << " pt " << cluster_pt[n] << " " << cluster_nmc_truth[n] << std::endl;

            Float_t truept = -999.9;
	    for(UInt_t counter = 0 ; counter<cluster_nmc_truth[n]; counter++)
	      {
                unsigned short index = cluster_mc_truth_index[n][counter];
	
                if(index!=65535){
		  // std::cout << mc_truth_pdg_code[index] << " " << mc_truth_first_parent_pdg_code[index] << std::endl;
		  if(mc_truth_pdg_code[index]==22 && mc_truth_first_parent_pdg_code[index]==22){
		    IsTrueHardPhoton = 1.0;
		    //std::cout <<cluster_pt[n] << " " << cluster_pt << " true: " << mc_truth_pt[index] << std::endl; 
		  }
                }
	      }//end loop over indices
            cluster_isHardPhoton[n] = IsTrueHardPhoton;
             //pi0 veto: cluster_minMass
            float minmass = 1.0;
            Int_t pi0tag = 0;

	    cluster_pi0tagged[n] =pi0tag; 
	    cluster_minMass[n] = minmass;
            cluster_SuperModule[n] = GetSuperModule(cluster_cell_id_max[n]);

	    if(cluster_pt[n]>ptmin){
              KeepEvent = true;
            }


            //Save Variables to veto high-pt tracks matching to clusters
            double dRmin = 1.0;
            double dRmin_trackpt = 0;
	    for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
	      if((track_quality[itrack]&16)==0) continue; //select only tracks that pass selection 16
	      if( not(track_pt[itrack]>5.0))    continue;
	      if( not(track_pt[itrack]<30.0))   continue;
	      Float_t deta =  cluster_eta[n]-track_eta[itrack];
	      Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack]);
	      Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
              if(dR<dRmin){
                  dRmin = dR;
                  dRmin_trackpt = track_pt[itrack];
	      }
	    }
            cluster_mindR_trackbit16_trackpt[n] = dRmin_trackpt;
            cluster_mindR_trackbit16[n] = dRmin;
            dRmin =1.0;
	    dRmin_trackpt = 999.0;
	    for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
              if((track_quality[itrack]&3)==0) continue; //select only tracks that pass selection 3 = ITS+TPC
              if( not(track_pt[itrack]>5.0)) continue;
	      if( not(track_pt[itrack]<30.0)) continue;
              Float_t deta =  cluster_eta[n]-track_eta[itrack];
              Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack]);
              Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
              if(dR<dRmin){
                  dRmin = dR;
	          dRmin_trackpt = track_pt[itrack];
	      }
            }
	    cluster_mindR_trackbit3_trackpt[n] = dRmin_trackpt;
	    cluster_mindR_trackbit3[n] = dRmin;

            //Calculate UE using "eta band method, that is used by Erwann et al"
            double iso_etaband =0.0;

            for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
                if((track_quality[itrack]&16)==0) continue; //select only tracks that pass selection 16 = ITS onluy
	        Float_t deta =  cluster_eta[n]-track_eta[itrack];
		Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack]);
		Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                if(TMath::Abs(dphi)>0.4) continue; // select tracks only within +/- 0.4 in phi (eta band)
                if(dR<0.4) continue; //rejects events inside cone
                iso_etaband = iso_etaband + track_pt[itrack];
	    }
            cluster_isoetaband_its[n] = iso_etaband;

            iso_etaband = 0.0;
            for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
	      if((track_quality[itrack]&3)==0) continue; //select only tracks that pass selection 3
	      Float_t deta =  cluster_eta[n]-track_eta[itrack];
	      Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack]);
	      Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
	      if(TMath::Abs(dphi)>0.4) continue; // select tracks only within +/- 0.4 in phi (eta band)
	      if(dR<0.4) continue; //rejects events inside cone
	      iso_etaband = iso_etaband + track_pt[itrack];
            }
	    cluster_isoetaband_tpc[n] = iso_etaband;
            cluster_iso_its_new[n]    = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n] - 0.4*0.4*TMath::Pi()*ue_estimate_its_const;

	    unsigned int cell_id_5_5[25];
	    cell_5_5(cell_id_5_5, cluster_cell_id_max[n]);

            
	    double eta_center = 0;
	    double phi_center = 0;
	    double sumw        = 0;
	    for (size_t i = 0; i < 25; i++) {
	      auto c_id = cell_id_5_5[i];
	      if( not(c_id <17664)) continue;
	      if( not(cell_e[c_id] >0.100)) continue;
	      eta_center =+ cell_e[c_id]* c_eta_array[c_id];
	      phi_center =+ cell_e[c_id]* c_phi_array[c_id];
	      sumw       =+ cell_e[c_id];
	    }
	    eta_center = eta_center/sumw;
	    phi_center = phi_center/sumw;
	    float ce  = 0;
	    float cep = 0;
	    float cp  = 0;
	    float m   = 0;
	    for(size_t i=0; i<25; i++){
	      auto c_id = cell_id_5_5[i];
              if( not(c_id <17664)) continue;
	      if( not(cell_e[c_id] >0.100)) continue;
              float w = TMath::Log(cell_e[c_id]/cluster_e[n]);
	      ce  += (c_eta_array[c_id] - eta_center) * (c_eta_array[c_id]- eta_center) * w;
	      cep += (c_eta_array[c_id] - eta_center) * (c_phi_array[c_id] - phi_center) * w;
	      cp  += (c_phi_array[c_id] - phi_center) * (c_phi_array[c_id] - phi_center) * w;
	      m   += w;
	    }
	    ce /= m;
	    cep /= m;
	    cp /= m;
	    cluster_b5x5[n] = 1000*0.5 * (ce + cp + TMath::Sqrt( TMath::Power(ce - cp, 2.0) + 4*TMath::Power(cep,2.0)));



	    ce  = 0;
            cep = 0;
            cp  = 0;
            m   = 0;
            for(size_t i=0; i<25; i++){
              auto c_id = cell_id_5_5[i];
              if( not(c_id <17664)) continue;
              if( not(cell_e[c_id] >0.100)) continue;
              float w = cell_e[c_id];
              ce  += (c_eta_array[c_id] - eta_center) * (c_eta_array[c_id]- eta_center) * w;
              cep += (c_eta_array[c_id] - eta_center) * (c_phi_array[c_id] - phi_center) * w;
              cp  += (c_phi_array[c_id] - phi_center) * (c_phi_array[c_id] - phi_center) * w;
              m   += w;
            }
            ce /= m;
            cep /= m;
            cp /= m;
            cluster_b5x5_lin[n] = 1000*0.5 * (ce + cp + TMath::Sqrt( TMath::Power(ce - cp, 2.0) + 4*TMath::Power(cep,2.0)));
	  
	  }//end looop clusters 
	  if (ievent % 25000 == 0) { std::cout << " event " << ievent << std::endl;}
	  if(ievent==0){
	    newtree->Fill();
	    //std::cout <<" cell e, eta , phi " << cell_e[0] << " " << cell_eta[0] << std::endl;
	  }

	  if(KeepEvent){
            newtree->Fill();
          }
	}

	//newtree->Print();
	std::cout << "#events passing skimming: " << newtree->GetEntries() << std::endl;
	newtree->AutoSave();
	delete file;
	delete newfile;
	  


    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
