class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalCorrectionTask;
class AliAnalysisGrid;


void LoadMacros();

void run2()
{
  Bool_t isLocal = kTRUE;
  //Setting track cuts
  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());
  AliTrackContainer::SetDefTrackCutsPeriod("lhc13d");

  LoadMacros();
  ///Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager();
  AliESDInputHandler* pESDHandler = AddESDHandler();

  //////Physics selection
  AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();
  ////Centrality (only works on ESD)
  AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
  pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  ////CDB ////////////////////////////////////////////////////////////
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  const UInt_t  kPhysSel       = AliVEvent::kEMCEGA + AliVEvent::kAnyINT;
  //std::cout <<" Loading AddTaskEmcalCorrectionTask.C" << std::endl;
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
  AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
  correctionTask->SelectCollisionCandidates(kPhysSel);
  correctionTask->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
  correctionTask->SetUserConfigurationFilename("EMCalConfig.yaml");
  correctionTask->Initialize();

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


  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetAliPhysicsVersion("vAN-20170617-1");
  // plugin->SetAliPhysicsVersion("vAN-20170519-1"); //this is maybe causing ISSUES.
  plugin->SetGridDataDir("/alice/data/2013/LHC13d");
  plugin->SetDataPattern("*/pass4/*/*ESDs.root");
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(195783);

 plugin->SetGridWorkingDir("workdir");
  plugin->SetGridOutputDir("output_chi01");

 plugin->SetSplitMaxInputFileNumber(40);
  plugin->SetExecutable("myTask.sh");
  plugin->SetTTL(5000);
  plugin->SetJDLName("AliAnalysisTaskEMCALPi0GammaCorr.jdl");
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetKeepLogs(kTRUE);
  plugin->SetMaxMergeStages(1);
  plugin->SetMergeViaJDL(kFALSE);

  mgr->SetGridHandler(plugin);

const char* alien_close_se = gSystem->Getenv("alien_CLOSE_SE");
  if(alien_close_se !=NULL){
    std::cout << " Setting Close SE " << std::endl;
    const char* file = mgr->GetCommonFileName();
    plugin->SetDefaultOutputs(kFALSE);
    plugin->SetOutputFiles(Form("%s@%s", file, alien_close_se));
    // plugin->SetOutputArchive(Form(
    //                            "log_archive.zip:stdout,stderr@%s"
    //                            "*.root@%s",
    //                            alien_close_se, alien_close_se));
  }
Bool_t isTesting = kFALSE;
  if(isTesting){
    std::cout<< "Running Grid Test " << std::endl;
    plugin->SetNtestFiles(1);
    plugin->SetRunMode("test"); //
  }
  else{
    plugin->SetRunMode("terminate"); //full, terminate
  }
  mgr->StartAnalysis("grid");
  return;
}

void LoadMacros()
{
  // Aliroot macros
  //  std::cout << "---***--- Adding ESDHandler" << std::endl;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/EMCALTasks/macros/AddTaskEMCALPi0GammaCorr.C");
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALPi0GammaCorr.C");
}
