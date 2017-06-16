class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalSetupTask;
class AliAnalysisGrid;

void runGrid()
{

  Bool_t isLocal = kFALSE;
  //Setting track cuts
  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());
  AliTrackContainer::SetDefTrackCutsPeriod("lhc13d"); 
  //since we will compile a class, tell root where to look for headers
  //std::cout <<"Proceesing lines to include $ROOTSYS/include, $ALICE_ROOT/include and $ALICE_PHYSICS/include" << std::endl;
  //gROOT->ProcessLine(".include $ROOTSYS/include");
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");
  //gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  
  ///Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager();
  std::cout << "---***--- Adding ESDHandler" << std::endl;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  AliESDInputHandler* pESDHandler = AddESDHandler();
    
  //////Physics selection
  std::cout << "---***--- Adding AliPhysicsSelection" << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();  
  ////Centrality (only works on ESD)
  std::cout << "---***--- Adding Centrality" << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
  pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  ////CDB ////////////////////////////////////////////////////////////
  std::cout << "---***--- Adding TaskCDBConnect" << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  const UInt_t  kPhysSel       = AliVEvent::kEMCEGA + AliVEvent::kAnyINT;
  std::cout <<" Loading AddTaskEmcalCorrectionTask.C" << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
  AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
  correctionTask->SelectCollisionCandidates(kPhysSel);
  correctionTask->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
  correctionTask->SetUserConfigurationFilename("EMCalConfig.yaml");
  correctionTask->Initialize();
  //////////////////////////////////////////////////////////////////
  
  std::cout << " About to load AddTaskEMCALPi0GammaCorr macro " << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/EMCALTasks/macros/AddTaskEMCALPi0GammaCorr.C");
  //create an instance of your analysis task
  AliAnalysisTaskEMCALPi0GammaCorr *ptr = AddTaskEMCALPi0GammaCorr();

  Printf("About to set Force Beam Type");
  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpA;
  bool bIsRun2 = false;
  
  TObjArray *pTopTasks = mgr->GetTasks();
  for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i)
	{
		AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
		if (!pTask) continue;
        if (pTask->InheritsFrom("AliAnalysisTaskEmcal"))
		{
			AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
			Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcal->GetName());
			pTaskEmcal->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
            if (bIsRun2) pTaskEmcal->SetUseNewCentralityEstimation(kTRUE);
            if (bIsRun2) pTaskEmcal->SetNCentBins(5);
		}
        if (pTask->InheritsFrom("AliEmcalCorrectionTask"))
        {
            AliEmcalCorrectionTask * pTaskEmcalCorrection = static_cast<AliEmcalCorrectionTask*>(pTask);
            Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcalCorrection->GetName());
            pTaskEmcalCorrection->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
            if (bIsRun2) pTaskEmcalCorrection->SetNCentBins(5);
        }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  Printf("About to initialize AliAnalysisManager::InitAnalysis()");
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  if(isLocal){ 
    TChain *chain = new TChain("esdTree");
    chain->Add("LHC13d_ESD1.root");
    chain->Add("LHC13d_ESD2.root");
    std::cout << "Running locally on tree, Chain has " << chain->GetEntries() << " events " << std::endl;
     mgr->StartAnalysis("local",chain);
     return;
  }

    
  AliAnalysisAlien *plugin = new AliAnalysisAlien()
  plugin->SetOverwriteMode();
 // plugin->SetAliPhysicsVersion("vAN-20170519-1"); //this is maybe causing ISSUES.
  plugin->SetGridDataDir("/alice/data/2013/LHC13d");
  plugin->SetDataPattern("*/pass4/*/*ESDs.root");
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(195783);

  plugin->SetGridWorkingDir("workdir");
  plugin->SetGridOutputDir("output_pi01");
  
  //plugin->AddIncludePath ("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");//noot needed.
  //plugin->SetAdditionalLibs("$ALICE_PHYSICS/PWGGA/EMCALTasks/AliAnalysisTaskEMCALPi0GammaCorr.cxx $ALICE_PHYSICS/PWGGA/EMCALTasks/AliAnalysisTaskEMCALPi0GammaCorr.h"); //commented because crashed when using aliphysics
  //plugin->SetAnalysisSource("$ALICE_PHYSICS/PWGGA/EMCALTasks/AliAnalysisTaskEMCALPi0GammaCorr.cxx");  //commented because crashed when using aliphysics
  
  plugin->SetSplitMaxInputFileNumber(40);
  plugin->SetExecutable("myTask.sh");
  plugin->SetTTL(10000);
  plugin->SetJDLName("AliAnalysisTaskEMCALPi0GammaCorr.jdl");  
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetKeepLogs(kTRUE);
  plugin->SetMaxMergeStages(1);
  plugin->SetMergeViaJDL(kTRUE);

  mgr->SetGridHandler(plugin);

  Bool_t isTesting = kTRUE;
  if(isTesting){
    std::cout<< "Running Grid Test " << std::endl;
    plugin->SetNtestFiles(1);
    plugin->SetRunMode("test"); // 
  }
  else{
    plugin->SetRunMode("full"); //full, terminate
  }  
  mgr->StartAnalysis("grid");
  return;
}

