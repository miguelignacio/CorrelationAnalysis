#include "TList.h"
#include "TFile.h"
#include "THnSparse.h"
#include <iostream>


void Reformat(){
    TFile* fIn = new TFile("ROOTOUTPUT/Results_060717.root","READ");
    fIn->Print();
    //TList* list = (TList*)fIn->Get("AliAnalysisTask_histos");
    AliEmcalList* list = (AliEmcalList*)fIn->Get("AliAnalysisTask_histos");
    THnSparse* h_Cluster = (THnSparse*)list->FindObject("h_Cluster");
    THnSparse* h_ClusterTrack = (THnSparse*)list->FindObject("h_ClusterTrack");
    THnSparse* h_ClusterTrack_Mixed = (THnSparse*)list->FindObject("h_ClusterTrack_Mixed");
    
    THnSparse* h_Pion = (THnSparse*)list->FindObject("h_Pi0");
    THnSparse* h_PionTrack = (THnSparse*)list->FindObject("h_Pi0Track");
    THnSparse* h_PionTrack_Mixed = (THnSparse*)list->FindObject("h_Pi0Track_Mixed");
    
    THnSparse* h_Track = (THnSparse*)list->FindObject("h_Track");
 
    TFile* fOut = new TFile("THnSparses.root","RECREATE");
    h_Cluster->Write("h_Cluster");
    h_ClusterTrack->Write("h_ClusterTrack");
    h_ClusterTrack_Mixed->Write("h_ClusterTrack_Mixed");
    
    h_Pion->Write("h_Pion");
    h_PionTrack->Write("h_PionTrack");
    h_PionTrack_Mixed->Write("h_PionTrack_Mixed");
    
    h_Track->Write("h_Track");
    fOut->Close();
    fIn->Close();
}