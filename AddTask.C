// $Id$
PionHadron* AddTask(

  Bool_t      InputDoMixing          = 1,                 //..same event=0 mixed event =1 (currenlty used to init the pool=1, throw out events without clusters=0)
  UInt_t      evtTriggerType         = AliVEvent::kEMCEGA, //AliVEvent::kAnyINT,// AliVEvent::kEMCEGA,//..use this type of events to combine gammas(trigger) with hadrons
  UInt_t      evtMixingType          = AliVEvent::kAnyINT,//..use only this type of events to fill your mixed event pool with tracks
  Double_t    trackptcut             = 0.0,              //..
  Double_t    clusptcut              = 0.30,              //..
  Bool_t      SavePool               = 0,                 //..saves a mixed event pool to the output event
  const char *trackName              = "usedefault",
  const char *clusName               = "usedefault",
  const char *taskname               = "AliAnalysisTask",
  const char *suffix                 = ""
)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTask", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTask", "This task requires an input event handler");
    return 0;
  }
 
  //..in case of AOD the default names are:
  if(trackName=="usedefault")trackName = "tracks";
  if(clusName =="usedefault")clusName  = "caloClusters";
  //-------------------------------------------------------
  // Built the name of the Task together
  //-------------------------------------------------------
 
  TString SameMixName;
  if(InputDoMixing == 0)
  {
	  SameMixName += "SE";
  }
  else
  {
	  SameMixName += "ME";
  }
  TString combinedName;
  combinedName.Form("%s_%s_%s_%s",taskname,(const char*)SameMixName,trackName,clusName);
  if(suffix!="")
  {
	  combinedName += "_";
	  combinedName += suffix;
  }
  cout<<"combinedName: "<<combinedName<<endl;
  TString contName(combinedName);
  contName += "_histos";
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  PionHadron* AnalysisTask = new PionHadron(InputDoMixing);
  //..Add the containers and set the names
  AnalysisTask->AddClusterContainer(clusName);
  if (trackName == "mcparticles")
  {
	  AliMCParticleContainer* mcpartCont = AnalysisTask->AddMCParticleContainer(trackName);
	  mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks")
  {
	  AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer(trackName);
	  trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
  }
  else  //implemented for testing correction framework
  {
	  AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer(trackName);
	  trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
  }
  //..check that condition!! maybe for mixed events its different!!!!!!
  if(!AnalysisTask->GetTrackContainer(trackName) || !AnalysisTask->GetClusterContainer(clusName))
  {
	 cout<<"Task can not run like this!"<<endl;
	 return 0;
  }
  //-------------------------------------------------------
  // Add some selection criteria
  //-------------------------------------------------------
  AnalysisTask->SetOffTrigger(evtTriggerType|evtMixingType); //..select only evets of type evtTriggerType and evtMixingType
  
    //..for Run1 pPb
  //AnalysisTask->SetUseManualEvtCuts(kTRUE);
  //AnalysisTask->SetUseAliAnaUtils(kTRUE);
  //AnalysisTask->SetVzRange(-10,10);
  //AnalysisTask->SetCentRange(0.0,100.0);
  //..new task for run2
  AnalysisTask->SetNCentBins(5);
  AnalysisTask->SetUseNewCentralityEstimation(kTRUE); //maybe this is what is required
  if(AnalysisTask->GetTrackContainer(trackName))
  {
	  AnalysisTask->GetTrackContainer(trackName)->SetParticlePtCut(trackptcut);
	 // AnalysisTask->GetTrackContainer(trackName)->SetParticleEtaLimits(-trackEta,trackEta); //..Eta limits (-0.8,0.8 as in Pi0-h publication)
  }
  Double_t phiToR=TMath::Pi()/180.0;
  if(AnalysisTask->GetClusterContainer(clusName))
  {
	  AnalysisTask->GetClusterContainer(clusName)->SetClusECut(0);                 //by default set to 0
	  AnalysisTask->GetClusterContainer(clusName)->SetClusPtCut(clusptcut);        //by default set to 0.15
	  AnalysisTask->GetClusterContainer(clusName)->SetClusUserDefEnergyCut(AliVCluster::kHadCorr,0);
	  AnalysisTask->GetClusterContainer(clusName)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
//    AnalysisTask->GetClusterContainer(clusName)->SetClusTimeCut(,);
	 // AnalysisTask->GetClusterContainer(clusName)->SetEtaLimits(-0.9,0.9);
//	  AnalysisTask->GetClusterContainer(clusName)->SetPhiLimits(68*phiToR,174*phiToR);
  }
  //..some additional input for the analysis
  AnalysisTask->SetNeedEmcalGeom(kTRUE);
  AnalysisTask->SetSavePool(SavePool);
  AnalysisTask->SetEvtTriggerType(evtTriggerType);   //..Trigger to be used for filling same event histograms
  AnalysisTask->SetEvtMixType(evtMixingType);        //..Trigger to be used to fill tracks into the pool (no GA trigger!!)
 
//	  AnalysisTask->SetRmvMatchedTrack(1);           //..Removes all clusters that have a matched track
 // }
  //for later AnalysisTask->SetEffHistGamma(THnF *h);
  //for later AnalysisTask->SetEffHistHadron(THnF *h);
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(AnalysisTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),TList::Class(),
		  	  	  	  	  	  	  	  	  	  	  	  	   AliAnalysisManager::kOutputContainer,
		  	  	  	  	  	  	  	  	  	  	  	  	   Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
  mgr->ConnectOutput (AnalysisTask, 1, coutput1 );
  return AnalysisTask;
}