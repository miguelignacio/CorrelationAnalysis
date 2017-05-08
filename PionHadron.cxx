//
#include <Riostream.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "PionHadron.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"
#include "AliEventPoolManager.h"
#include "AliMCEvent.h"

#include <memory>
using std::cout;
using std::endl;

ClassImp(PionHadron)

////////////////////////////////////////////////////////////////////////////////////////
PionHadron::PionHadron():
AliAnalysisTaskEmcal("PionHadron", kTRUE),

fDoMixing(0),
fDebug(0),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fRtoD(0),
fMixBCent(0),
fMixBZvtx(),
fPoolMgr(0x0),
fTrackDepth(0),
fPoolSize(0),
fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7),
fMixingEventType(AliVEvent::kINT7),
fCurrentEventTrigger(0),
fParticleLevel(kFALSE),
fIsMC(kFALSE),
fEventCutList(0),
thisEvent(),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0)
{
	InitArrays();
	for(int i=0;i<nMulClass;i++){
	  for(int j=0;j<nZClass;j++){
	    for(int k=0;k<nPtClass;k++){
	      iEvt[i][j][k] = 0;
	      for(int l=0;l<nEvt;l++){
		EmcEventList[i][j][k][l].SetGlobalInfo(0,0.,0.);
	      }
	    }
	  }
	}
}

// -------------------------------------------------------------------------------------
// Constructor with inputs
PionHadron::PionHadron(Bool_t InputDoMixing):
AliAnalysisTaskEmcal("PionHadron", kTRUE),
fDoMixing(0),
fDebug(0),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fRtoD(0),
fMixBCent(0),
fMixBZvtx(),
fPoolMgr(0x0),
fTrackDepth(0),
fPoolSize(0),
fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), 
fMixingEventType(AliVEvent::kINT7),
fCurrentEventTrigger(0),
fParticleLevel(kFALSE),
fIsMC(kFALSE),
fEventCutList(0),
thisEvent(),
h_Cluster(0),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0)
{
	InitArrays();
	fDoMixing          =InputDoMixing;

	for(int i=0;i<nMulClass;i++){
	  for(int j=0;j<nZClass;j++){
	    for(int k=0;k<nPtClass;k++){
	      iEvt[i][j][k] = 0;
	      for(int l=0;l<nEvt;l++){
		EmcEventList[i][j][k][l].SetGlobalInfo(0,0.,0.);
	      }
	    }
	  }
	}
}//End constructor PiHadron that receives input


void PionHadron::InitArrays()
{
	fDoMixing          =0; //= 0 (do only same event analyis with correct triggers), =1 (do event mixing)
	fMCorData          =0; // 0->MC, 1->Data
	fDebug             =0; //set only 1 for debugging
	fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.
	fUseManualEventCuts=0; //=0 use automatic setting from AliEventCuts. =1 load manual cuts
    fRtoD=180.0/TMath::Pi();
    
    //Setting bins for the mixing of events.
    Double_t centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	fMixBCent = new TAxis(kNcentBins,centmix);

    Double_t zvtxmix[kNvertBins+1] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
	fMixBZvtx = new TAxis(kNvertBins,zvtxmix);
    fTrackDepth     = 50000;    //Raymonds/Megans value
    fPoolSize       = 1; 
    SetMakeGeneralHistograms(kTRUE); //from aliEMCAL TASK
} //end of function init arrays


PionHadron::~PionHadron()
{

//   delete fOutput;
}

void PionHadron::UserCreateOutputObjects()
{
    AliWarning("Entering UserCreateOutPutObjects");   
	AliAnalysisTaskEmcal::UserCreateOutputObjects();
    fEventCutList = new TList();
	fEventCutList ->SetOwner();
	fEventCutList ->SetName("EventCutOutput");
	fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
    
    if(fUseManualEventCuts==1)
	{  
    AliWarning("Entering UserCreateOutPutObjects");   
    fEventCuts.SetManualMode();
    
    fEventCuts.fCentralityFramework=1; //..only for Run1!!
	fEventCuts.fTriggerMask = fOffTrigger;
	fEventCuts.fMinVtz = fMinVz;
	fEventCuts.fMaxVtz = fMaxVz;
	fEventCuts.fRequireTrackVertex = true;
	fEventCuts.fMaxDeltaSpdTrackAbsolute=fZvertexDiff;
	fEventCuts.fTrackletBGcut = fTklVsClusSPDCut; //(false by default for 15o)
	fEventCuts.fMinCentrality = fMinCent;
	fEventCuts.fMaxCentrality = fMaxCent;
    }
    fEventCuts.AddQAplotsToList(fEventCutList);
	fOutput->Add(fEventCutList);
    
    //initializing event mixer:
    if(fDoMixing==1 || fPoolMgr) //do this for either a mixed event analysis or when an external pool is given
	{
        AliWarning("Initializing Event Mixer");
		InitEventMixer();
	}
    else{  AliWarning("Did not initialize Event Mixer"); }
    //Initializing the histograms to be saved. For the moment, only pT of clusters and Mgammagamma.
    int nbins_Mass = 20;
    int nbins_Pt   = 50;
    int nbins_E    = 50;
    int nbins_dphi     = 18;
    int nbins_deta     = 20;
    int nbins_phi     = 18;
    int nbins_eta     = 20;
    int nbins_zt      = 50;
    int nbins_xi      = 50;
    int nbins_M02     = 50;
    int nbins_Ncells  = 20;
    int nbins_Centrality = 5;
    
    double min_Mass = 0;
    double max_Mass = 0.5;
    double min_Pt = 0;
    double max_Pt = 10.0;
    double min_dphi = -0.5; // rads
    double max_dphi = 1.5; // rads 
    double min_phi = -1.0*TMath::Pi();
    double max_phi = TMath::Pi();

    double max_deta = 2.0;
    double max_eta =  1.0;
    double min_deta = -max_deta;
    double min_eta =  -max_eta;
    
    double min_zt = 0;
    double min_xi = 0;
    double max_zt = 2.0;
    double max_xi = 2.0;
    
    double min_E =0;
    double max_E =10;
    double min_M02 = 0.0;
    double max_M02 = 1.0;
    
    double min_Ncells = -0.5;
    double max_Ncells = 19.5;
     
    double min_Centrality = 0;
    double max_Centrality = 100;
    Int_t bins[22]    = {nbins_Centrality, nbins_Mass, nbins_Pt,  nbins_eta, nbins_phi, nbins_E, nbins_Pt,  nbins_eta, nbins_phi, nbins_dphi, nbins_deta, nbins_deta/2, nbins_zt, nbins_xi,
                        nbins_Pt, nbins_Pt, nbins_eta, nbins_eta, nbins_phi, nbins_phi, nbins_M02, nbins_M02}; //these are for the photons from pion decays
    Double_t xmin[22] = {min_Centrality, min_Mass,   min_Pt  , min_eta    , min_phi , min_E, min_Pt, min_eta, min_phi    , min_dphi   , min_deta   , 0.0, min_zt, min_xi,
                        min_Pt, min_Pt, min_eta, min_eta, min_phi, min_phi, min_M02, min_M02};
    Double_t xmax[22] = {max_Centrality, max_Mass,   max_Pt  , max_eta    , max_phi , max_E, max_Pt, max_eta,  max_phi   , max_dphi,   max_deta, max_deta, max_zt, max_xi,
                        max_Pt, max_Pt, max_eta, max_eta, max_phi, max_phi, max_M02, max_M02};

    TString axisNames = "Pion-Track THnSparse; Centrality; #pion Mass; #pionpT; #pion Eta; #pion phi; #pion E";
    axisNames = axisNames + "; track_pT; track Eta; track Phi; #Dphi ; #Deta; #|Deta|; Zt; Xi";
    axisNames = axisNames + "; ph1_pT; ph2_pT; ph1_eta; ph2_eta; ph1_phi; ph2_phi; ph1_M02; ph2_M02";

    h_Pi0Track = new THnSparseD("h_Pi0Track", axisNames, 22, bins, xmin,xmax);
    fOutput->Add(h_Pi0Track);
    h_Pi0Track_Mixed = new THnSparseD("h_Pi0Track_Mixed", axisNames, 21, bins, xmin,xmax);
    fOutput->Add(h_Pi0Track_Mixed);
    
    axisNames = "Pion THnSparse; Centrality; #pion Mass; #pionpT; #pion Eta; #pion phi; #pion E";
    axisNames = axisNames + "; ph1_pT; ph2_pT; ph1_eta; ph2_eta; ph1_phi; ph2_phi; ph1_M02; ph2_M02";
    Int_t binsPi0[14] = {nbins_Centrality, nbins_Mass, nbins_Pt, nbins_eta, nbins_phi, nbins_E, 
                         nbins_Pt, nbins_Pt, nbins_eta, nbins_eta, nbins_phi, nbins_phi, nbins_M02, nbins_M02};
    
    Double_t xminPi0[14] = {min_Centrality, min_Mass, min_Pt, min_eta, min_phi , min_E, min_Pt, min_Pt, min_eta, min_eta, min_phi, min_phi, min_M02, min_M02};
    Double_t xmaxPi0[14] = {max_Centrality, max_Mass, max_Pt, max_eta , max_phi , max_E,max_Pt, max_Pt, max_eta, max_eta, max_phi, max_phi, max_M02, max_M02};
    h_Pi0= new THnSparseD("h_Pi0", axisNames, 14, binsPi0, xminPi0, xmaxPi0);
    fOutput->Add(h_Pi0);
    
    axisNames = "Cluster THnSparse; Centrality; Cluster E; Cluster Pt; Cluster Eta; Cluster Phi; Cluster M02; Cluster NCells";
    Int_t binsCluster[7] = {nbins_Centrality, nbins_E, nbins_Pt, nbins_eta, nbins_phi, nbins_M02, nbins_Ncells};
    Double_t xminCluster[7] = {min_Centrality, min_E, min_Pt, min_eta, min_phi, min_M02, min_Ncells};
    Double_t xmaxCluster[7] = {max_Centrality, max_E, max_Pt, max_eta, max_phi, max_M02, max_Ncells};
    h_Cluster = new THnSparseD("h_Cluster", axisNames, 7, binsCluster, xminCluster, xmaxCluster);
    fOutput->Add(h_Cluster);
    
	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}



void PionHadron::InitEventMixer()
{
	Int_t nCentBins=fMixBCent->GetNbins();
	Double_t centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(Int_t i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}
	//..Z-vertex pools
	Int_t nZvtxBins=fMixBZvtx->GetNbins();
	Double_t zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}
	//..in case no external pool is provided create one here
	if(!fPoolMgr)
	{
		cout<<"....  Pool Manager Created ...."<<endl;
		fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin);
		fPoolMgr->SetTargetValues(fTrackDepth, 0.1, 5);  //pool is ready at 0.1*fTrackDepth = 5000 or events =5
		//save this pool by default
	}
	else
	{
		fPoolMgr->ClearPools();
		cout<<"....  Pool Manager Provided From File ...."<<endl;
	}
	//..Check binning of pool manager (basic dimensional check for the time being) to see whether external pool fits the here desired one??
	if( (fPoolMgr->GetNumberOfMultBins() != nCentBins) || (fPoolMgr->GetNumberOfZVtxBins() != nZvtxBins) )
	{
		AliFatal("Binning of given pool manager not compatible with binning of correlation task!");
	}

	if(fSavePool==1)
	{
	fPoolMgr->SetSaveFlag(-1, 10000, -10000, 100000, 0, 0, -1, 10000000);
	fOutput->Add(fPoolMgr);
	}
	//..Basic checks and printing of pool properties
	fPoolMgr->Validate();
}



void PionHadron::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt)
{
	//..This allows you to add only specific pools and not the full pool manager to the output file
	std::vector<Double_t> binVec;
	binVec.push_back(minCent);
	binVec.push_back(maxCent);
	binVec.push_back(minZvtx);
	binVec.push_back(maxZvtx);
	binVec.push_back(minPt);
	binVec.push_back(maxPt);
	fEventPoolOutputList.push_back(binVec);
}


void PionHadron::ExecOnce()
{
    AliAnalysisTaskEmcal::ExecOnce();
}


Bool_t PionHadron::IsEventSelected()
{
    //std::cout << " Is Event Selected function " << std::endl;
	//..checks: rejects DAQ incompletes, magnetic field selection, offline trigger
	//..vertex selection, SPD pile-up (if enabled), centrality cuts
	if (!fEventCuts.AcceptEvent(InputEvent()))
	{
		PostData(1, fOutput);
        //AliWarning("Did not pass selection ");
		return kFALSE;
	}
    //std::cout << " Centrality " << fEventCuts.GetCentrality() << std::endl;
    
    TString Trigger;
    Trigger = fInputEvent->GetFiredTriggerClasses();
    bool PassedGammaTrigger = kFALSE;
    bool PassedMinBiasTrigger = kFALSE;
    
    if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
    if(Trigger.Contains("INT7")) PassedMinBiasTrigger = kTRUE;
    
    if(PassedGammaTrigger) AliWarning("Passed Gamma Trigger ");
    if(PassedMinBiasTrigger) AliWarning("Passed MinBias Trigger ");
    //std::cout << Trigger << std::endl;
    if(!PassedGammaTrigger && !PassedMinBiasTrigger) return kFALSE;
//..should be DELETED!!! later
   

	AliAnalysisTaskEmcal::IsEventSelected();
	return kTRUE;
}


Bool_t PionHadron::Run()
{
    //std::cout<< " Run function " << std::endl;
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	//..Determine the trigger for the current event
	//fCurrentEventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	//..check here some basic properties of the event
	//if (fDoMixing==0 && !fCaloClusters)                         return kFALSE;
	//if (fDoMixing==0 && !(fCurrentEventTrigger & fTriggerType)) return kFALSE;
	return kTRUE;
}

Bool_t PionHadron::FillHistograms()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	// 1. First get an event pool corresponding in mult (cent) and
	//    zvertex to the current event. Once initialized, the pool
	//    should contain nMix (reduced) events. This routine does not
	//    pre-scan the chain. The first several events of every chain
	//    will be skipped until the needed pools are filled to the
	//    specified depth. If the pool categories are not too rare, this
	//    should not be a problem. If they are rare, you could lose
	//    statistics.
	// 2. Collect the whole pool's content of tracks into one TObjArray
	//    (bgTracks), which is effectively a single background super-event.
	// 3. The reduced and bgTracks arrays must both be passed into
	//    FillCorrelations(). Also nMix should be passed in, so a weight
	//    of 1./nMix can be applied.

    //..Get pool containing tracks from other events like this one
    
    //AliWarning("Filling Histograms");
    TString Trigger;
    Trigger = fInputEvent->GetFiredTriggerClasses();
    bool PassedGammaTrigger = kFALSE;
    bool PassedMinBiasTrigger = kFALSE;
    
    if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
    if(Trigger.Contains("INT7")) PassedMinBiasTrigger = kTRUE;
    

	Double_t zVertex = fVertex[2];
    //AliWarning(Form("Centrality %f, ZVertex %f, FEventCuts centrality %f",fCent, zVertex,fEventCuts.GetCentrality() ));
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);

    fCent = fEventCuts.GetCentrality();
    AliWarning(Form("Centrality %f, ZVertex %f, FEventCuts centrality %f",fCent, zVertex,fEventCuts.GetCentrality() ));
	AliEventPool* pool = 0x0;
	pool = fPoolMgr->GetEventPool(fCent, zVertex);
	if (!pool)
	{
    	AliWarning(Form("No pool found. Centrality %f, ZVertex %f",fCent, zVertex));
		return kFALSE;
	}
	if(pool->IsReady() && PassedGammaTrigger)
	{
		AliWarning(Form("Pool ready. Centrality %f, ZVertex %f",fCent, zVertex));
		Int_t nMix = pool->GetCurrentNEvents();
		cout<<"number of events in pool: "<<nMix<<endl;
		for(Int_t jMix=0; jMix<nMix; jMix++)
		{
			TObjArray* bgTracks=0x0;
			bgTracks = pool->GetEvent(jMix);
			if(!bgTracks)
			{
				cout<<"could not retrieve TObjArray from EventPool!"<<endl;
			}
	    	//..Loop over clusters and fill histograms
            std::cout << " Correlated Cluster and Track, MIXED EVENTS " << std::endl;
			CorrelateClusterAndTrack(0,bgTracks,0,1.0/nMix);//correlate with mixed event
		}
	}
    if(PassedMinBiasTrigger && !PassedGammaTrigger )//&& ((fCurrentEventTrigger & 1000000000000000)==0))
	{
        TObjArray* tracksClone=0x0;
		tracksClone = CloneToCreateTObjArray(tracks);
		//..if there is no track object or the pool is locked do not update
		if(tracksClone && !pool->GetLockFlag())
		{
          //AliWarning(Form("Updating pool. Centrality %f, ZVertex %f",fCent, zVertex));
		  //AliWarning("Updating Pool");
          pool->UpdatePool(tracksClone);
		}
	}
   
     if(PassedGammaTrigger)
	{
    //Now do same-event analysis:
        std::cout << " Correlated Cluster and Track " << std::endl;
        CorrelateClusterAndTrack(tracks,0,1,1);//correlate with same event
    }
}



TObjArray* PionHadron::CloneToCreateTObjArray(AliParticleContainer* tracks)
{

	//..clones a track list
	if(!tracks)                            return 0;
	if(tracks->GetNAcceptedParticles()==0) return 0;
	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);
	Int_t NoOfTracksInEvent =tracks->GetNParticles();
	AliVParticle* track=0;
	for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
	{
		track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
		if(!track)continue; 
		tracksClone->Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!! Major error!!!! "<<"Accepted tracks in event: "<<tracks->GetNAcceptedParticles()<<", Tracks in TObjArray: "<<tracksClone->GetEntries()<<endl;
	return tracksClone;
}



//This function is the core of the analysis
Int_t PionHadron::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
    std::cout << " Number of clusters " <<  NoOfClustersInEvent << std::endl;
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight=1;    //weight to normalize mixed and same event distributions individually
	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* track=0;
    Weight=InputWeight; //..for mixed events normalize per events in pool
	Int_t GammaCounter=0;
    
    //AliWarning(Form("Number of Clusters %d",NoOfClustersInEvent));

	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ ) // Loop over clusters
	{   
		std::cout <<NoCluster1;
        cluster=(AliVCluster*) clusters->GetCluster(NoCluster1); // //it was GetAcceptCluster->GetCluster(NoCluster1);
     	if(!cluster) continue;
        if(!AccClusterForAna(clusters,cluster))continue ; 
        FillClusterHisto(cluster, h_Cluster);
        
        for( Int_t NoCluster2 = NoCluster1+1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
        {
            // std::cout <<NoCluster2 << std::endl;
            cluster2=(AliVCluster*) clusters->GetCluster(NoCluster2);
			if(!cluster2) continue;
            if(!AccClusterForAna(clusters,cluster2))continue ; 
			
            FillPionHisto(cluster, cluster2, h_Pi0); //filling variables
            if(SameMix==0){
		        for(Int_t ibg=0; ibg<bgTracksArray->GetEntries(); ibg++){
	    	        AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
		            if(!track) continue;
                    FillCorrelation(cluster, cluster2, track, h_Pi0Track_Mixed);
                }
            }
            else{
                for(Int_t NoTrack = 0; NoTrack < tracks->GetNParticles(); NoTrack++){ //correlate pion with tracks
				    track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
                    if(!track) continue;
                    FillCorrelation(cluster, cluster2, track, h_Pi0Track);
                    }
            }
        } //end 2 loop over clusters
	} //end  1 loop over clusters
	return GammaCounter; //return number of accepted gammas
}


void  PionHadron::FillCorrelation(AliVCluster* cluster1, AliVCluster* cluster2, AliVParticle* track, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph_1, ph_2, pi0; 
    //std::cout << " getting momentum";
	clusters->GetMomentum(ph_1, cluster1);
    clusters->GetMomentum(ph_2, cluster2);
	pi0= ph_1+ph_2;
       
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi = TVector2::Phi_mpi_pi(pi0.Phi()- trackphi);
    dphi = dphi/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    double deta = pi0.Eta()-track->Eta();
    double  Zt  = track->Pt()/pi0.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    double entries[21] = {pi0.M(), pi0.Pt(), pi0.Eta(), pi0.Phi(), pi0.E(), track->Pt(), track->Eta(), trackphi, dphi, deta, abs(deta), Zt, Xi,
                          ph_1.Pt(), ph_2.Pt(), ph_1.Eta(), ph_2.Eta(), ph_1.Phi(), ph_2.Phi() , cluster1->GetM02(), cluster2->GetM02()};                
    //std::cout<<"About to fill histogram with entries";
    histo->Fill(entries);
    
    //delete clusters;
    return;
}

void  PionHadron::FillPionHisto(AliVCluster* cluster1, AliVCluster* cluster2, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph_1, ph_2, pi0; 
    //std::cout << " getting momentum";
	clusters->GetMomentum(ph_1, cluster1);
    clusters->GetMomentum(ph_2, cluster2);
	pi0= ph_1+ph_2;
       
    double entries[13] = {pi0.M(), pi0.Pt(), pi0.Eta(), pi0.Phi(), pi0.E(), 
                          ph_1.Pt(), ph_2.Pt(), ph_1.Eta(), ph_2.Eta(), ph_1.Phi(), ph_2.Phi() , cluster1->GetM02(), cluster2->GetM02()};                
    histo->Fill(entries);
    
    //delete clusters;
    return;
}

void PionHadron::FillClusterHisto(AliVCluster* cluster, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph;
    clusters->GetMomentum(ph, cluster);
    
    double entries[7] = {fCent, ph.E(), ph.Pt(), ph.Eta(), ph.Phi(), cluster->GetM02(), cluster->GetNCells()};                
    histo->Fill(entries);
    
    //delete clusters;
    return;
}

TObjArray* PionHadron::CloneClustersTObjArray(AliClusterContainer* clusters)
{
	//..clones a track list
	if(!clusters)                            return 0;
	if(clusters->GetNClusters()==0) return 0;
	TObjArray* clustersCloneI = new TObjArray;
	clustersCloneI->SetOwner(kTRUE);
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	AliVCluster* cluster = 0;
	for(Int_t NoClus = 0; NoClus < NoOfClustersInEvent; NoClus++)
	{
		cluster = (AliVCluster*) clusters->GetAcceptCluster(NoClus);
		if(!cluster)continue; //check if the Cluster is good
		clustersCloneI->Add((AliVCluster*)cluster);// Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(clustersCloneI->GetEntries()!=clusters->GetNAcceptedClusters())cout<<"!!!!!!! Major error!!!! "<<"Accepted clusters in event: "<<clusters->GetNAcceptedClusters()<<", Tracks in TObjArray: "<<clustersCloneI->GetEntries()<<endl;
	return clustersCloneI;
}

 void PionHadron::GetMulClassPi0(Int_t& imcl) //What exacly is the use of this??
 {
   Int_t nclus = 0;
   AliClusterContainer* clustersCont  = GetClusterContainer(0);
   if (!clustersCont) return;
   TObjArray* clustersClone=0x0;
   clustersClone = CloneClustersTObjArray(clustersCont);
   TObjArray *clusters =clustersClone;// fEsdClusters;
   // if (!clusters)   clusters = fAodClusters;
   if (clusters)
   nclus = clusters->GetEntries();
   const int MultCut[nMulClass] = {5, 12, 20 ,50, 100, 250, 500, 9999};
   imcl=0;
   for (imcl=0; imcl<nMulClass; imcl++) {
     if (nclus < MultCut[imcl]) break;
   }
 }
 
  void PionHadron::GetZVtxClassPi0(Int_t& ivcl)
 {
   Double_t ZvTx=fEventCuts.GetPrimaryVertex()->GetZ();
   if (!ZvTx) return;
   const Double_t VtxCut[nZClass] = {-20, -5, -1.5, 2.5, 6, 20};
   ivcl=0;
   for (ivcl=0; ivcl<nMulClass; ivcl++) {
     if (ZvTx < VtxCut[ivcl]) break;
   }
 }



void PionHadron::AddMixEventPi0(const Int_t MulClass, const Int_t vtxClass, const Int_t PtClass, Int_t& iEvent, const Float_t& phitrig, const Float_t& thetatrig)
{
  if(iEvent >= nEvt){
    iEvent = 0;
  }
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  AliClusterContainer* clustersCont  = GetClusterContainer(0);
  if (!clustersCont) return;
  TObjArray* clustersClone=0x0;
  clustersClone = CloneClustersTObjArray(clustersCont);
  TObjArray *clusters =clustersClone;

  if (!clusters) return;
  Int_t nclus = 0;
  nclus = clusters->GetEntries();
  nclus = thisEvent.nHits;

  thisEvent.SetGlobalInfo(nclus,phitrig,thetatrig);
  EmcEventList[MulClass][vtxClass][PtClass][iEvent] = thisEvent;
  iEvent++;
  return;
}

Bool_t PionHadron::AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster)
{
    if(!caloCluster->IsEMCAL()) return kFALSE;
    if(caloCluster->E()<1.0) return kFALSE;
	//..at least 2 cells in cluster
	if(caloCluster->GetNCells()<2) return kFALSE;
	if(caloCluster->GetNExMax() > 1) return kFALSE; //local maxima should be 0 or 1
	if(caloCluster->GetM02()<0.1) return kFALSE;
    //..remove clusters with a matched track
	//if(fRmvMTrack==1 && caloCluster->GetNTracksMatched(s)!=0) return kFALSE;
	return kTRUE;
}

Double_t PionHadron::DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec)
{
	Double_t Phi_g = ClusterVec.Phi_0_2pi();
	Double_t Phi_h = TrackVec->Phi();
	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();
	dPhi = Phi_g-Phi_h;
	//--shift the second peak over the fist peak: \----/   --> ---
	//--to create a histogram that starts at -pi/2 and ends at 3/2pi
	if (dPhi <= -TMath::Pi()/2)    dPhi += 2*pi;
	if (dPhi > 3.0*TMath::Pi()/2.0)dPhi -= 2*pi;
	//--change from rad to degree:
	dPhi*= fRtoD;
	return dPhi;
}


Double_t PionHadron::GetEff(AliTLorentzVector ClusterVec)
{
	return 1;
}


//__________________________________________________________________________________________________
EmcEventPi0::EmcEventPi0()
: nHits(0),
TrigPhi(0),
TrigTheta(0)
{
  nHits = 0;
  TrigPhi = 0;
  TrigTheta = 0;
}
//__________________________________________________________________________________________________
EmcEventPi0::EmcEventPi0(const EmcEventPi0 &obj)
: nHits(0),
TrigPhi(0),
TrigTheta(0)
{
  nHits = obj.nHits;
  TrigPhi = obj.TrigPhi;
  TrigTheta = obj.TrigTheta;
  // copy all hits
  for(int i=0;i<nHits;i++){
    hit[i].hittype = obj.hit[i].hittype;
    hit[i].imo = obj.hit[i].imo;
    hit[i].pid = obj.hit[i].pid;
    hit[i].weight = obj.hit[i].weight;
    hit[i].thishit = obj.hit[i].thishit;
    hit[i].NCells = obj.hit[i].NCells;
    hit[i].smno = obj.hit[i].smno;
    hit[i].CellRay = obj.hit[i].CellRay;
  }
}
//__________________________________________________________________________________________________
void EmcEventPi0::SetGlobalInfo(const Int_t& Size, const Float_t& phiTrig, const Float_t& thetaTrig)
{
  //    fCenPercent = centPer;
  //    fVtx = vtxPos;
  nHits = Size;
  TrigPhi = phiTrig;
  TrigTheta = thetaTrig;
}

void EmcEventPi0::Print()
{
  Printf("%d hits",nHits);
  for(int i=0;i<nHits;i++){
    hit[i].thishit.Print();
  }
}
//__________________________________________________________________________________________________
void EmcEventPi0::Reset()
{
  for(int i=0;i<nHits;i++){
    hit[i].hittype = 0;
    hit[i].imo = 0;
    hit[i].pid = 0;
    hit[i].weight = 1;
    hit[i].smno = 0;
    TLorentzVector lv(0,0,0,0);
    hit[i].thishit = lv;
    hit[i].NCells = 0;
    hit[i].CellRay.clear();
  }
  nHits = 0;
}
//__________________________________________________________________________________________________
EmcHitPi0::EmcHitPi0()
  : thishit(),
    hittype(0),
    imo(0),
    pid(0),
    weight(1),
    bclean(1),
    NCells(0),
    CellRay()
{
}
