#ifndef PIONHADRON_H
#define PIONHADRON_H

// $Id$
#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include <THn.h>
#include "AliStaObjects.h"
#include "AliEventCuts.h"

class TH1;
class TH2;
class TH3;
class THnSparse;
class THnSparse;
class AliVVZERO;
class AliEvtPoolManager;

using std::vector;

// class to store EMC hits    <<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
class EmcHitPi0 {
  TLorentzVector thishit;
  Short_t hittype; // 100 for bg event/real event, 1 for added pi0, 2 for added eta, -1 for not needed; from primary pi0: 101, from secondary pi0: 102, from K0: 103, from material: 104
  Byte_t smno;
  Int_t imo; // index of original mother in monte carlo stack
  Int_t pid; // particle ID
  Float_t weight; // weight from mother particle
  Bool_t bclean; // clean if only one contributor
  Int_t NCells;
  std::vector<int> CellRay;
public:
  //virtual ~EmcHit();
  EmcHitPi0();
  friend class EmcEventPi0;
  friend class PionHadron;
  void Print(){Printf("E=%.2f, type=%d, MoID=%d, PID=%d, w=%.3f",thishit.E(),hittype,imo,pid,weight);   }
};

class EmcEventPi0 {
  //    Int_t fCenPercent;
  //    Int_t fVtx;
	Float_t TrigPhi; // phi of highest pT hit on EMCal
  Float_t TrigTheta; // eta of highest pT hit ...
  const static int nMaxHit = 1000;
  int nHits;
  EmcHitPi0 hit[nMaxHit];
public:
  EmcEventPi0();
  EmcEventPi0(const EmcEventPi0 &obj);
  //virtual ~EmcEvent();
  //    void SetGlobalInfo(const Int_t&, const Int_t&, const Int_t&, const Int_t&, const Double_t&, const Double_t&);
  void SetGlobalInfo(const Int_t&, const Float_t&, const Float_t&);
  int evsize() {return nHits;}
  void Reset();
  void Print();
  friend class PionHadron;
};


class PionHadron : public AliAnalysisTaskEmcal {
 public:
	PionHadron();
	PionHadron(Bool_t InputSameEventAnalysis);
virtual ~PionHadron();

  void SetEffHistGamma(THnF *h)                              { fHistEffGamma    = h      ; }
  void SetEffHistHadron(THnF *h)                             { fHistEffHadron   = h      ; }
  void SetSavePool(Bool_t input)                             { fSavePool        = input  ; }
  void SetEvtTriggerType(UInt_t input)                       { fTriggerType     = input  ; }
  void SetEvtMixType(UInt_t input)                           { fMixingEventType = input  ; }
  void SetUseManualEvtCuts(Bool_t input)                     { fUseManualEventCuts = input;}

  //Functions for mixed event purposes
  void SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  // Set which pools will be saved
  void AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);
  
  private:
  AliEventCuts fEventCuts;                   ///< event selection utility
  
  void                        InitArrays()                                                 ;
  // overwritten EMCal framework functions
  Bool_t                      Run()                             	                          ;
  void                        ExecOnce()         									      ;
  Bool_t                      IsEventSelected()											  ;
  void                        UserCreateOutputObjects()        		                      ;
  
  //Functions for mixed event purposes
  void                        InitEventMixer()											  ;
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;

  Bool_t                      FillHistograms()                                              ;
  
  void                        FillClusterHisto(AliVCluster* cluster, THnSparse* histo);
  void                        FillPionHisto(AliVCluster* cluster1, AliVCluster* cluster2, THnSparse* histo);
  void                        FillCorrelation(AliVCluster* cluster1, AliVCluster* cluster2, AliVParticle* track, THnSparse* histo);
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Bool_t                      AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster);
  
   //<<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
    TObjArray*                  CloneClustersTObjArray(AliClusterContainer* clusters)          ;
    void GetMulClassPi0(Int_t&);
    void GetZVtxClassPi0(Int_t&);
    void AddMixEventPi0(const Int_t, const Int_t, const Int_t, Int_t&, const Float_t&, const Float_t&);
  //<<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
  
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec) ;
  Double_t                    GetEff(AliTLorentzVector ParticleVec)                         ;

  Bool_t                      fDoMixing;                 ///< This option enables mixed events being used in the analysi
  Bool_t                      fMCorData;                 //<Are we looking at simulations or at the real thing
  Bool_t                      fDebug;			        ///< Can be set for debugging
  Bool_t                      fSavePool;                 ///< Defines whether to save output pools in a root file
  Bool_t                      fUseManualEventCuts;       ///< Use manual cuts if automatic setup is not available for the period
  
  
  
  //..Input histograms
  THnF                       *fHistEffGamma;             ///< ??input efficiency for trigger particles
  THnF                       *fHistEffHadron;            ///< ??input efficiency for associate particles

  //..Constants
  Double_t                    fRtoD;                     ///< conversion of rad to degree
  static const Int_t          kNIdentifier=3;            ///< number of different versions of the same histogram type, can later be used for centrality or mixed event eg.
  static const Int_t          kNvertBins=20;             ///< vertex bins in which the ME are mixed
  static const Int_t          kNcentBins=8;              ///< centrality bins in which the ME are mixed
  static const Int_t          kNoGammaBins=9;            ///< Bins in gamma pT
  static const Int_t          kNoZtBins=7;               ///< Bins in Zt
  static const Int_t          kNoXiBins=8;               ///< Bins in Xi
  Double_t                    fArray_G_Bins[10];         ///< 10=kNoGammaBins+1
  Double_t                    fArray_ZT_Bins[8];         ///< 8=kNoZtBins+1
  Double_t                    fArray_XI_Bins[9];         ///< 9=kNoXiBins+1
  Double_t                    fArrayNVertBins[21];       ///< 21=kNvertBins+1

  //..Event pool variables
  TAxis                      *fMixBCent;                 ///< Number of centrality bins for the mixed event
  TAxis                      *fMixBZvtx;                 ///< Number of vertex bins for the mixed event
  AliEventPoolManager        *fPoolMgr;                  ///< event pool manager
  Int_t                       fTrackDepth;               ///<  #tracks to fill pool
  Int_t                       fPoolSize;                 ///<  Maximum number of events
  vector<vector<Double_t> >   fEventPoolOutputList;      //!<! ???vector representing a list of pools (given by value range) that will be saved
  //..Event selection types
  UInt_t                      fTriggerType;              ///<  Event types that are used for the trigger (gamma or pi0)
  UInt_t                      fMixingEventType;          ///<  Event types that are used for the tracks in the mixed event
  UInt_t                      fCurrentEventTrigger;      //!<! Trigger of the current event

  // MC stuff
  Bool_t                      fParticleLevel;            ///< Set particle level analysis
  Bool_t                      fIsMC;                     ///< Trigger, MC analysis
  UInt_t                      fAODfilterBits[2];         ///< AOD track filter bit map
  
   // THnSparse
  THnSparse                 *h_Cluster;                 //!<!
  THnSparse                 *h_Pi0;                 //!<!
  THnSparse                 *h_Pi0Track;                 //!<! THnSparse with info on pi0 and track.
  THnSparse                 *h_Pi0Track_Mixed;                 //!<!
  // Other stuff
  TList                      *fEventCutList;           //!<! Output list for event cut histograms
  TList                      *OutputList;            //!<! Output list
 
  
  const static int nMulClass =   8;  // <<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
  const static int nZClass   =   6;
  const static int nPtClass = 1;
  int iEvt[nMulClass][nZClass][nPtClass];
  const static int nEvt      =   10;//30; // mixing "depth"

  EmcEventPi0 evt;
  EmcEventPi0 EmcEventList[nMulClass][nZClass][nPtClass][nEvt];

  EmcEventPi0 thisEvent;

 private:
  PionHadron(const PionHadron&);            // not implemented
  PionHadron &operator=(const PionHadron&); // not implemented

  ClassDef(PionHadron, 10) // Class to analyse gamma hadron correlations
};
#endif



