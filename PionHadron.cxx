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

fGammaOrPi0(0),
fDoMixing(0),
fDebug(0),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),
fClShapeMax(10),
fMaxNLM(10),
fRmvMTrack(0),
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
fHistNoClusPt(0),
fHistPi0(0),
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
PionHadron::PionHadron(Bool_t InputGammaOrPi0,Bool_t InputDoMixing):
AliAnalysisTaskEmcal("PionHadron", kTRUE),
fGammaOrPi0(0),
fDoMixing(0),
fDebug(0),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),
fClShapeMax(10),
fMaxNLM(10),
fRmvMTrack(0),
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
fHistNoClusPt(0),
fHistPi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0)
{
	InitArrays();
	fGammaOrPi0        =InputGammaOrPi0;
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
	fGammaOrPi0        =0; //= 0 ( Gamma analysis ), 1 (pi0 analysis)
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
    //OutputList= new TList();
    //OutputList->SetOwner();
    //OutputList->SetName("NewOutput");

	AliAnalysisTaskEmcal::UserCreateOutputObjects();
    fEventCutList = new TList();
	fEventCutList ->SetOwner();
	fEventCutList ->SetName("EventCutOutput");
	fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
    
    if(fUseManualEventCuts==1)
	{  
    fEventCuts.SetManualMode();
    
    fEventCuts.fCentralityFramework=2; //..only for Run1!!
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
		InitEventMixer();
	}
    
   
    //Initializing the histograms to be saved. For the moment, only pT of clusters and Mgammagamma.
    
    fHistNoClusPt = new TH1F(Form("fHistNoClusPt_Id%0d",1),Form("fHistNoClusPt_Id%0d",1), 31,0, 31);
	fHistNoClusPt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClusPt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClusPt->GetBinWidth(0)));
	fOutput->Add(fHistNoClusPt);
    
    
    fHistPi0 = new TH1F(Form("fHistPi0_%0d",1),Form("fHistPi0_%0d",1), 500, 0, 0.5);
	fHistPi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistPi0->GetYaxis()->SetTitle("Entries");
    //fOutput->Add(fHistPi0);
	//fOutput->Add(fHistPi0);    
    int nbins_Mass = 500;
    int nbins_Pt   = 50;
    int nbins_E    = 50;
    int nbins_dphi     = 18;
    int nbins_deta     = 20;
    int nbins_phi     = 18;
    int nbins_eta     = 20;
    int nbins_zt      = 50;
    int nbins_xi      = 50;
    
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
     
    Int_t bins[13]    = {nbins_Mass, nbins_Pt,  nbins_eta, nbins_phi, nbins_E, nbins_Pt,  nbins_eta, nbins_phi, nbins_dphi, nbins_deta, nbins_deta/2, nbins_zt, nbins_xi};
    Double_t xmin[13] = {min_Mass,   min_Pt  , min_eta    , min_phi , min_E, min_Pt, min_eta, min_phi    , min_dphi   , min_deta   , 0.0, min_zt, min_xi};
    Double_t xmax[13] = {max_Mass,   max_Pt  , max_eta    , max_phi , max_E, max_Pt, max_eta,  max_phi   , max_dphi,   max_deta, max_deta, max_zt, max_xi};

    TString axisNames = "Pion-Track THnSparse; #pion Mass; #pionpT; #pion Eta; #pion phi; #pion E; track_pT; track Eta; track Phi; #Dphi ; #Deta; #|Deta|; Zt; fxi";

    h_Pi0Track = new THnSparseD("h_Pi0Track", axisNames, 13, bins, xmin,xmax);
    fOutput->Add(h_Pi0Track);
    h_Pi0Track_Mixed = new THnSparseD("h_Pi0Track_Mixed", axisNames, 13, bins, xmin,xmax);
    fOutput->Add(h_Pi0Track_Mixed);
    
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
		return kFALSE;
	}
    
    //..Start of copy part from AliAnalysisTaskEmcal
	if (!fTrigClass.IsNull())
	{
		TString fired;
		const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
		if (aev){ fired = aev->GetFiredTriggerClasses(); }
		else cout<<"Error analysis only for AODs"<<endl;
		if (!fired.Contains("-B-")) 
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
		std::unique_ptr<TObjArray> arr(fTrigClass.Tokenize("|"));
		if (!arr)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
		Bool_t match = 0;
		for (Int_t i=0;i<arr->GetEntriesFast();++i)
		{
			TObject *obj = arr->At(i);
			if (!obj)
				continue;
			//Check if requested trigger was fired
			TString objStr = obj->GetName();
			if(fEMCalTriggerMode == kOverlapWithLowThreshold &&
					(objStr.Contains("J1") || objStr.Contains("J2") || objStr.Contains("G1") || objStr.Contains("G2"))) {
				// This is relevant for EMCal triggers with 2 thresholds
				// If the kOverlapWithLowThreshold was requested than the overlap between the two triggers goes with the lower threshold trigger
				TString trigType1 = "J1";
				TString trigType2 = "J2";
				if(objStr.Contains("G"))
				{
					trigType1 = "G1";
					trigType2 = "G2";
				}
				if(objStr.Contains(trigType2) && fired.Contains(trigType2.Data()))
				{ //requesting low threshold + overlap
					match = 1;
					break;
				}
				else if(objStr.Contains(trigType1) && fired.Contains(trigType1.Data()) && !fired.Contains(trigType2.Data())) { //high threshold only
					match = 1;
					break;
				}
			}
			else
			{
				// If this is not an EMCal trigger, or no particular treatment of EMCal triggers was requested,
				// simply check that the trigger was fired
				if (fired.Contains(obj->GetName()))
				{
					match = 1;
					break;
				}
			}
		}
		if (!match)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
	}
	if (fTriggerTypeSel != kND)
	{
		if (!HasTriggerType(fTriggerTypeSel))
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigTypeSel",1);
			return kFALSE;
		}
	}
//..should be DELETED!!! later
	AliAnalysisTaskEmcal::IsEventSelected();
	return kTRUE;
}


Bool_t PionHadron::Run()
{
    //std::cout<< " Run function " << std::endl;
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	//..Determine the trigger for the current event
	fCurrentEventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	//..check here some basic properties of the event
	if (fDoMixing==0 && !fCaloClusters)                         return kFALSE;
	if (fDoMixing==0 && !(fCurrentEventTrigger & fTriggerType)) return kFALSE;
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
	Double_t zVertex = fVertex[2];
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);

    if(fDoMixing==1)
	{
		AliEventPool* pool = 0x0;
		pool = fPoolMgr->GetEventPool(fCent, zVertex);
		if (!pool)
		{
			AliWarning(Form("No pool found. Centrality %f, ZVertex %f",fCent, zVertex));
			return kFALSE;
		}

    //..Start combining triggers (from fTriggerType events)
		//..with a pool filled with tracks (from fMixingEventType)
		if(pool->IsReady() && (fCurrentEventTrigger & fTriggerType))
		{
			//..get number of current events in pool
			Int_t nMix = pool->GetCurrentNEvents();
			//cout<<"number of events in pool: "<<nMix<<endl;
			for(Int_t jMix=0; jMix<nMix; jMix++)
			{
				TObjArray* bgTracks=0x0;
				bgTracks = pool->GetEvent(jMix);
				if(!bgTracks)
				{
					cout<<"could not retrieve TObjArray from EventPool!"<<endl;
				}
				//..Loop over clusters and fill histograms
				CorrelateClusterAndTrack(0,bgTracks,0,1.0/nMix);//correlate with mixed event
			}
		}
        if((fCurrentEventTrigger & fMixingEventType) && ((fCurrentEventTrigger & 1000000000000000)==0))
		{
			TObjArray* tracksClone=0x0;
			tracksClone = CloneToCreateTObjArray(tracks);
			//..if there is no track object or the pool is locked do not update
			if(tracksClone && !pool->GetLockFlag())
			{
				pool->UpdatePool(tracksClone);
			}
		}
	} //end of if DoMixing
    
    //Now do same-event analysis:
    CorrelateClusterAndTrack(tracks,0,1,1);//correlate with same event
    
    //Get all EMCAL clusters in the event
    AliVCluster* cluster = 0;
    AliClusterContainer* clusters  = GetClusterContainer(0); 
    if (!clusters) return 0;
    Int_t NoOfClustersInEvent = clusters->GetNClusters();
    
    thisEvent.Reset();
	Int_t NAccClus=0;
	Float_t Max_Phi=0;
	Float_t Max_Eta=0;

    //Loop over clusters in the event and apply selection
    for(Int_t NoCluster0 = 0; NoCluster0 < NoOfClustersInEvent; NoCluster0++ )
	  {
	    cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster0); //->GetCluster(NoCluster1);
	    if(!cluster || !AccClusterForAna(clusters,cluster) || cluster->E()<=2)continue; //check if the cluster is a good cluster
	    NAccClus++;
	  }
	thisEvent.SetGlobalInfo(NAccClus,Max_Phi,Max_Eta);
    
    //Loop over clusters in the event again  (?) ,Get the number of accepted clusters and put it in thisEvent (???)
    Int_t AccClus=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	  {
	  cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
	  if(!cluster || !AccClusterForAna(clusters,cluster) || cluster->E()<=2)continue; //check if the cluster is a good cluster
	  AccClus++;
	    TLorentzVector CaloClusterVec;
	    clusters->GetMomentum(CaloClusterVec, cluster);
	    thisEvent.hit[AccClus-1].thishit=CaloClusterVec;
	    thisEvent.hit[AccClus-1].NCells=cluster->GetNCells();
        //Loop over the cells within a cluster
	    for(Int_t c = 0; c<cluster->GetNCells(); c++){
	      thisEvent.hit[AccClus-1].CellRay.push_back(cluster->GetCellsAbsId()[c]); 
	     }
       }

    Int_t nclusthis = thisEvent.nHits;
	Int_t vtxClass = 1;
	Int_t MulClass = 4;
	GetMulClassPi0(MulClass);
	GetZVtxClassPi0(vtxClass);
    
    Float_t phitrig = 0;
	Float_t thetatrig = 0;
	Double_t pt_max = 0;
	Int_t ptClass = 0;
    
    AddMixEventPi0(MulClass, vtxClass, ptClass, iEvt[MulClass][vtxClass][ptClass], Max_Phi, Max_Eta);
	return kTRUE;
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
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight=1;    //weight to normalize mixed and same event distributions individually
	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* track=0;
    Weight=InputWeight; //..for mixed events normalize per events in pool
	Int_t GammaCounter=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ ) // Loop over clusters
	{   
        //std::cout << " Photon loop " << std::endl;
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster)continue; //check if the cluster is a good cluster
		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec, cluster);
		AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(CaloClusterVec); //..can acess phi from
		EffWeight_Gamma=GetEff(aliCaloClusterVec);
    	if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
        
        fHistNoClusPt->Fill(CaloClusterVec.Pt());
        GammaCounter++;
        for( Int_t NoCluster2 = NoCluster1+1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
        {
            //std::cout << " Pion loop " << std::endl;
            cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
			if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster
			TLorentzVector CaloClusterVec2;
			TLorentzVector vpi0;
			clusters->GetMomentum(CaloClusterVec2, cluster2);
			vpi0=CaloClusterVec+CaloClusterVec2;
            
            fHistPi0->Fill(vpi0.M()); //filling diphoton invariant mass distribution
			
            if(SameMix==0){
		        for(Int_t ibg=0; ibg<bgTracksArray->GetEntries(); ibg++){
	    	        AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
		            if(!track) continue;
                    FillCorrelation(vpi0, track, h_Pi0Track_Mixed);
                }
            }
            else{
                for(Int_t NoTrack = 0; NoTrack < tracks->GetNParticles(); NoTrack++){ //correlate pion with tracks
				    track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
                    if(!track) continue;
                    FillCorrelation(vpi0, track, h_Pi0Track);
                    }
            }
        } //end 2 loop over clusters
	} //end  1 loop over clusters
	return GammaCounter; //return number of accepted gammas
}


void  PionHadron::FillCorrelation(TLorentzVector vpi0, AliVParticle* track, THnSparse* histo){
    
    //check if the track is a good track
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi = TVector2::Phi_mpi_pi(vpi0.Phi()- trackphi);
    dphi = dphi/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    double deta = vpi0.Eta()-track->Eta();
    double  Zt  = track->Pt()/vpi0.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    double entries[13] = {vpi0.M(), vpi0.Pt(), vpi0.Eta(), vpi0.Phi(), vpi0.E(), track->Pt(), track->Eta(), trackphi, dphi, deta, abs(deta), Zt, Xi};                
    histo->Fill(entries);
    
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
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	Double_t deltaPhi=2;   //..phi away from detector edges.
	Double_t deltaEta=0.0; //..eta away from detector edges.
    //!!!! eventually transform to AliTLorentzvector
	//..Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=1; //..By default accepted
	//-----------------------------
    //Require at least 2 GeV of energy:
    //if(caloCluster->E()<0.50) return 0;
	//..at least 2 cells in cluster
	if(caloCluster->GetNCells()<2)		return 0;
	//..number of local maxima should be 1 or 0 (for cluster splitting this has to be changed)
	if(caloCluster->GetNExMax()>fMaxNLM) return 0;
	//..cut on the cluster shape
	if(fClShapeMin>0 && fClShapeMax>0   && (caloCluster->GetM02()<fClShapeMin || caloCluster->GetM02()>fClShapeMax)) return 0;
	//..remove clusters with a matched track
	if(fRmvMTrack==1 && caloCluster->GetNTracksMatched()!=0) return 0;
	return Accepted;
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
