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
#include <stdio.h>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>


void to_sm_nphi(unsigned int &sm, unsigned int &nphi,
		unsigned int n)
{
  sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
  nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
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

void cell_3_3(unsigned int n_3_3[], const unsigned int n,
	      const unsigned int ld = 3)
{
  unsigned int sm;
  unsigned int nphi;

  to_sm_nphi(sm, nphi, n);

  if (n % 2 == 0) {
    n_3_3[0 * ld + 0] = n - 1;
    n_3_3[0 * ld + 1] = n + 1;
    n_3_3[0 * ld + 2] = n + 3;
  }
  else {
    n_3_3[0 * ld + 0] = n - 2 * nphi - 3;
    n_3_3[0 * ld + 1] = n - 2 * nphi - 1;
    n_3_3[0 * ld + 2] = n - 2 * nphi + 1;
  }
  n_3_3[1 * ld + 0] = n - 2;
  n_3_3[1 * ld + 1] = n;
  n_3_3[1 * ld + 2] = n + 2;
  if (n % 2 == 0) {
    n_3_3[2 * ld + 0] = n + 2 * nphi - 1;
    n_3_3[2 * ld + 1] = n + 2 * nphi + 1;
    n_3_3[2 * ld + 2] = n + 2 * nphi + 3;
  }
  else {
    n_3_3[2 * ld + 0] = n - 3;
    n_3_3[2 * ld + 1] = n - 1;
    n_3_3[2 * ld + 2] = n + 1;
  }
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

        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

        //TTree *_tree_event = dynamic_cast<TTree *>
        //    (dynamic_cast<TDirectoryFile *>
        //     (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << " fail " << std::endl;
	  exit(EXIT_FAILURE);
        }  
        //_tree_event->Print();

	TApplication application("", &dummyc, dummyv);
	TCanvas canvas("canvas", "");//, 960 + 4, 720 + 28);
	//canvas.SetLogy();

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

        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];   
        Float_t cluster_b5x5_lin[NTRACK_MAX];


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
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
	_tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
        _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
 	_tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max); 
        _tree_event->SetBranchAddress("cluster_b5x5_lin",cluster_b5x5_lin);
        _tree_event->SetBranchAddress("cell_e", cell_e);

	 
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
 	
 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
    
	int counter =0;

        //define output files:
        FILE * pFile[10];
        pFile[0] = fopen ("output_All.txt","w");
        pFile[1] = fopen ("output_1.txt","w");
	pFile[2] = fopen ("output_2.txt","w");
	pFile[3] = fopen ("output_3.txt","w");
	pFile[4] = fopen ("output_4.txt","w");
	pFile[5] = fopen ("output_5.txt","w");
	pFile[6] = fopen ("output_LambdaPeak.txt","w");
	pFile[7] = fopen ("output_MC.txt","w");
        pFile[8] = fopen ("output_MC_DNNPeak.txt","w");
	pFile[9] = fopen ("output_MC_DNNTail.txt","w");
 
   

	//for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){     
        for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
             _tree_event->GetEntry(ievent);
	     bool stop = false;
	      for (ULong64_t n = 0; n < ncluster; n++) {
		  
	          if(TMath::Abs(primary_vertex[2])>10) continue;
	          if(cluster_pt[n]<11.0) continue;
                  //if(cluster_pt[n]>16.0) continue;
                  //if(TMath::Abs(cluster_eta[n])>0.6) continue;
                  if(cluster_ncell[n]<3) continue;  //ncell
                  if(cluster_e_cross[n]/cluster_e[n]<0.05) continue; //exotic clusters
                  //if(cluster_s_nphoton[n][1]>0.86) continue; //eliminate weirdos. 

	          unsigned int cell_id_5_5[25];
	          cell_5_5(cell_id_5_5, cluster_cell_id_max[n]);

                  int nbin = 0; //0
                  //if(cluster_s_nphoton[n][1] > 0.75 and cluster_s_nphoton[n][1] <0.85)      nbin=8;
		  //if(cluster_s_nphoton[n][1] > 0.55 and cluster_s_nphoton[n][1] <0.70)      nbin=9;

		  //else if(cluster_s_nphoton[n][1] > 0.18 and cluster_s_nphoton[n][1] <0.19)      nbin=2;
                  //else if(cluster_s_nphoton[n][1] > 0.19 and cluster_s_nphoton[n][1] <0.20) nbin=3;
                  //else if(cluster_s_nphoton[n][1] > 0.21 and cluster_s_nphoton[n][1] <0.22) nbin=4;
                  //else if(cluster_s_nphoton[n][1] > 0.22 and cluster_s_nphoton[n][1] <0.23) nbin=5;
		  //else if(cluster_lambda_square[n][0] > 0.249 and cluster_lambda_square[n][0] <0.251) nbin=6;
                  //else nbin=0;                  

                  //if(!(cluster_lambda_square[n][0]<0.252 and cluster_lambda_square[n][0]>0.248)) continue;
		  //if(cluster_s_nphoton[n][1] >0.85) continue; 
                  //if(cluster_s_nphoton[n][1] <0.75) continue;

                  fprintf (pFile[nbin], "%2.1f, %2.2f, %2.2f", cluster_pt[n], cluster_s_nphoton[n][1], cluster_b5x5_lin[n]);
		  //std::cout << Form("%2.1f, %2.3f, %2.3f", cluster_pt[n], cluster_s_nphoton[n][1], cluster_lambda_square[n][0]);
	          
                  for (size_t i = 0; i < 25; i++) {
                    auto cellenergy = cell_e[cell_id_5_5[i]];
                    double weight = TMath::Log(cellenergy/cluster_e[n]);
		    //std::cout << cellenergy/cell_e[n] << std::endl;
                    if(cellenergy<0.100) cellenergy = std::numeric_limits<double>::quiet_NaN();
		     //if(weight <-4.5){
		    //cellenergy = std::numeric_limits<double>::quiet_NaN();
		     //}
                    fprintf (pFile[nbin], ",%2.2f", 100*cellenergy/cluster_e[n]);
		  }
                  
                  fprintf (pFile[nbin],"\n");
		  ///		  std::cout<<std::endl;
	    }//end loop on clusters. 
	    if(stop) break;

	} //end loop over events

    }//end loop over samples

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
