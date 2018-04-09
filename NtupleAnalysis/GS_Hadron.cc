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
#include "H5Cpp.h"

#define HI_TREE "AliAnalysisTaskNTGJ/_tree_event"
#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

using namespace H5;

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout<<"Argument syntax is [Command] [root_file] [hdf5_file]"<<std::endl;
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

    std::cout << "Opening: " << (TString)argv[1] << std::endl;
    TFile *file = TFile::Open((TString)argv[1]);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
    //TTree *_tree_event = dynamic_cast<TTree *>(file->Get(HI_TREE));

    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }  
    //_tree_event->Print();

    TApplication application("", &dummyc, dummyv);

    //Tree Stuff    
    UInt_t ntrack;
    UInt_t ncluster;
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    Float_t track_eta_emcal[NTRACK_MAX];
    Float_t track_phi_emcal[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 
    //Long64_t Mix_Events[10];
    Long64_t LimitUse_Mixed_Events[10];

    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta",track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_eta_emcal",track_eta_emcal);
    _tree_event->SetBranchAddress("track_phi_emcal",track_phi_emcal);
    _tree_event->SetBranchAddress("track_quality", track_quality);
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    //_tree_event->SetBranchAddress("Mix_Events", Mix_Events);
    _tree_event->SetBranchAddress("LimitUse_Mixed_Events", LimitUse_Mixed_Events);

    //TH2D *Corr = new TH2D("Correlation", "Gale Shapley Mixed H-H Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
    TH2D Corr = TH2D("Correlation", "GS Mixed H-H Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
    TCanvas *Corr2D = new TCanvas("Corr2D", "GS_Mix_HH", 800,800);
    Corr.Sumw2();
    Corr.SetMinimum(0.); //set Y-Axis minimum to 0 when drawing 

    Long64_t nentries = _tree_event->GetEntries();

    //open hdf5
    const H5std_string hdf5_file_name(argv[2]);

    const H5std_string track_ds_name( "track" );
    H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY );
    DataSet track_dataset = h5_file.openDataSet( track_ds_name );
    DataSpace track_dataspace = track_dataset.getSpace();

    const H5std_string cluster_ds_name( "cluster" );
    DataSet cluster_dataset = h5_file.openDataSet( cluster_ds_name );
    DataSpace cluster_dataspace = cluster_dataset.getSpace();

    UInt_t ntrack_max = 0; //set to zero and find below
    UInt_t ncluster_max = 0;
       
      for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
      _tree_event->GetEntry(i);
      ntrack_max = std::max(ntrack_max, ntrack);
      ncluster_max = std::max(ncluster_max, ncluster);
      fprintf(stderr, "\r%s:%d: %llu", __FILE__, __LINE__, i);
      }

      fprintf(stderr, "\r%s:%d: %iu %iu", __FILE__, __LINE__, ntrack_max, ncluster_max);   

    //Define array hyperslab will be read into    
    float track_data_out[1][ntrack_max][7];
    float cluster_data_out[1][ncluster_max][5];

    //Define hyperslab size and offset in the FILE;
    hsize_t track_offset[3] = {0, 0, 0};           
    hsize_t track_count[3] = {1, ntrack_max, 7};   
    hsize_t cluster_offset[3] = {0, 0, 0};           
    hsize_t cluster_count[3] = {1, ncluster_max, 5};   

    track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
    cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

    //Define the memory dataspace to place hyperslab
    const int RANK_OUT = 3;
    hsize_t track_dimsm[3] = {1, ntrack_max, 7};      //memory space dimensions
    DataSpace track_memspace( RANK_OUT, track_dimsm );
    hsize_t cluster_dimsm[3] = {1, ncluster_max, 5};      //memory space dimensions
    DataSpace cluster_memspace( RANK_OUT, cluster_dimsm );

    //Define memory offset for hypreslab->array
    hsize_t track_offset_out[3] = {0};
    hsize_t cluster_offset_out[3] = {0};

    //define 2D memory hyperslab
    hsize_t track_count_out[3] = {1, ntrack_max, 7};    // size of the hyperslab in memory
    hsize_t cluster_count_out[3] = {1, ncluster_max, 5};

    //define space in memory for hyperslab, then write from file to memory
    track_memspace.selectHyperslab( H5S_SELECT_SET, track_count_out, track_offset_out );
    track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "track dataset read into array: OK");

    cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
    cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "cluster dataset read into array: OK");

    const int TrackCutBit =16;
    for (Long64_t ievent = 0; ievent < nentries; ievent++){
      //for (Long64_t ievent = 0; ievent < 1000; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
      for (Long64_t imix = 0; imix < 10; imix++){
	//Long64_t mix_event = Mix_Events[imix];
	Long64_t mix_event = LimitUse_Mixed_Events[imix];
	if (mix_event == ievent) continue;

	for (unsigned int n = 0; n < ntrack; n++){
	  //if(track_pt[n] < 4) continue;
	  if((track_pt[n] < 8) || (track_pt[n] > 16)) continue;

	  if((track_quality[n]&TrackCutBit)==0) continue;
	  bool HasMatch = false;

	  for (unsigned int c = 0; c < ncluster; c++) {
	    if (TMath::Abs(cluster_eta[c] - track_eta_emcal[n]) < 0.015) {
	      HasMatch = true;
	      break; }
	    if (TMath::Abs(cluster_phi[c] - track_phi_emcal[n]) < 0.015) {
	      HasMatch = true;
	      break; }
	  }
	  if (HasMatch) continue;

	  //read mixed event into hyperslab Place in track loop if high pt cut.
	  track_offset[0]=mix_event;
	  track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
	  track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
	  
	  cluster_offset[0]=mix_event;
	  cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
	  cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
	  
	  for (unsigned int i = 0; i < ntrack_max; i++){
	    if(std::isnan(track_data_out[0][i][0])) break;
	    if (track_data_out[0][i][1] < 0.15) continue;
	    //if (track_data_out[0][i][1] < 1) continue;
	    if (HasMatch) continue;
	    if ((int(track_data_out[0][i][4]+0.5)&TrackCutBit)==0) continue;

	    bool Mix_HasMatch = false;
	    for (unsigned int l = 0; l < ncluster_max; l++){
	      if (std::isnan(cluster_data_out[0][l][0])) break;
	      if (TMath::Abs(cluster_data_out[0][l][2] - track_data_out[0][i][5]) < 0.015) {
		Mix_HasMatch = true;
		break; }
	      if (TMath::Abs(cluster_data_out[0][l][3] - track_data_out[0][l][6]) < 0.015) {
		Mix_HasMatch = true;
		break; }
	    }
	    if (Mix_HasMatch) continue;	    

	    //FIXME: loop through clusters, see if matches
	    //don't forget to skip if emcal_eta = -999, or track_emcal_phi = 0.026....
	    

	    Float_t DeltaPhi = track_phi[n] - track_data_out[0][i][3];
	      if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi
	      if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
	    Float_t DeltaEta = track_eta[n] - track_data_out[0][i][2];
	    if ((TMath::Abs(DeltaPhi) < 0.01) && (TMath::Abs(DeltaEta) < 0.01)) continue;

	    Corr.Fill(DeltaPhi,DeltaEta);

	  }//Mixed tracks

	  if (ievent%10000==0){
	    std::cout<<std::endl;
	    for (unsigned int i = 0; i < ncluster_max; i++){
	      if (std::isnan(cluster_data_out[0][i][0])) break;
	      std::cout<<"Mixed Event: "<<imix<<":  ";
	      for(unsigned int j = 0; j < 5; j++)
		std::cout<<cluster_data_out[0][i][j]<<" ";
	      std::cout<<std::endl;
	    }
	    std::cout<<std::endl;
	  }//will repeat if offset has not been adjusted due to pT cut

	}//Main tracks
      }//imix
      if(ievent % 1000 == 0){
      Corr.Draw("SURF1");                                                                                                                         
      Corr2D->Update();                                                                                                                            
      }
    }//ievent

    Corr.Draw("SURF2");
    // TString NewFileName = argv[1];
//     NewFileName.ReplaceAll(".root", "");
//     NewFileName += "_Mixed.root";
    TString NewFileName = "LBin_Corr_8to16GeV.root";
    //TString NewFileName = "LBin_Corr_150MeV.root";
    TFile *MyFile = new TFile(NewFileName,"RECREATE");
    Corr.Write();
    MyFile->Print(); //write to file so I don't have to run this again
}
