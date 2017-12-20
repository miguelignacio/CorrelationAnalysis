#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TProfile.h>

#include <iostream>
#include <fstream>

#define NTRACK_MAX (1U << 14)

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }

    int dummyc = 1;
    char **dummyv = new char *[1];

    dummyv[0] = strdup("main");



    TApplication application("", &dummyc, dummyv);
    TCanvas canvas("canvas", "", 960 + 4, 720 + 28);
    TH1D histogram0("histogram0", "", 40, 10.0, 50.0);
    TH1D histogram1("histogram1", "", 100, 0.0, 2.0);
   

    //canvas.SetLogy();
    histogram0.Sumw2();
    histogram0.SetMarkerStyle(20);
    histogram0.SetMarkerColor(EColor::kBlack);
    histogram0.SetLineColor(EColor::kBlack);
    histogram0.SetXTitle("track pt (GeV)");


    for (int iarg = 1; iarg < argc; iarg++) {
	std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);

        if (file == NULL) {
	  std::cout << " fail" << std::endl;
	  exit(EXIT_FAILURE);
        }
        file->Print();

        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << " fail " << std::endl;
          exit(EXIT_FAILURE);
        }

        UInt_t ntrack;
	// ULong64_t ntrack;
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
        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];


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
        _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
        _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        
      
        
	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
        const int TrackCutBit=3;
        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
             _tree_event->GetEntry(i);

	    for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
              if(track_pt[itrack]<10.0) continue;
	      if((track_quality[itrack]&TrackCutBit)==0) continue; //select only tracks that pass selection 16
              float dRmin = 999;
              float eop   = 0;
	      for (ULong64_t iclus = 0; iclus < ncluster; iclus++) {
                  if( not(TMath::Abs(cluster_eta[iclus])<0.6)) continue; //cut edges of detector
                  if( not(cluster_ncell[iclus]>2)) continue;   //removes clusters with 1 or 2 cells
                  if( not(cluster_e_cross[iclus]/cluster_e[iclus]>0.05)) continue; //removes "spiky" clusters
                  float deta =  cluster_eta[iclus]-track_eta[itrack];
	          float dphi =  TVector2::Phi_mpi_pi(cluster_phi[iclus]-track_phi[itrack]);
	          float dR = TMath::Sqrt(deta*deta + dphi*dphi);
                  if(dR<dRmin){
                      dRmin = dR; 
                      eop   =  cluster_e[iclus]/track_pt[itrack];
		  }
	      }
              if(not(dRmin<999)) continue;
	      histogram0.Fill(track_pt[itrack]);
              if(eop>2.0) eop = 1.99;
              histogram1.Fill(eop);

	    }//end cluster loop
	
            if (i % 1000 == 0) {
                histogram1.Draw("e1x0");
                canvas.Update();
		std::cout << "Event # " << i << " / " << _tree_event->GetEntries() << std::endl;
            }
        }
    }

    histogram0.Draw("e1x0");
    histogram1.Draw("e1x0same");
    canvas.Update();

    application.Run();

    return EXIT_SUCCESS;
}
