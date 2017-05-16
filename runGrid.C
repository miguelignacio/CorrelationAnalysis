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

  //Setting track cuts
  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());
  AliTrackContainer::SetDefTrackCutsPeriod("lhc13d"); //this was not here
  //since we will compile a class, tell root where to look for headers
  std::cout << "command to ignore warnings" << std::endl;
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;") ;
  std::cout <<"Proceesing lines to include $ROOTSYS/include, $ALICE_ROOT/include and $ALICE_PHYSICS/include" << std::endl;
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  //create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager();
  //For AODs:  
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  std::cout << "---***--- Adding ESDHandler" << std::endl;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  //AliAODInputHandler *inputHandler = new AliAODInputHandler();
  //AliESDInputHandler *inputHandler = new AliESDInputHandler();
  AliESDInputHandler* pESDHandler = AddESDHandler();
  //  inputHandler->SetReadFriends(kFALSE);
  //mgr->SetInputEventHandler(inputHandler);
  
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
  
  
  /*////EmcalCorrection/////////////////////////////////////////////////
  std::cout <<" Loading AddTaskEmcalCorrectionTask.C" << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
  const UInt_t  kPhysSel       = AliVEvent::kEMCEGA + AliVEvent::kAnyINT;
  AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask();
  correctionTask->SelectCollisionCandidates(kPhysSel);
  correctionTask->SetRunPeriod(""); //it was "lhc13d"
  correctionTask->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
  correctionTask->SetUserConfigurationFilename("userGH_Analysis_LHC13.yaml");
  correctionTask->Initialize();
  //////////////////////////////////////////////////////////////////
  */
  //compile the class (locally)
  gROOT->LoadMacro("PionHadron.cxx++g");
  //load the addtask macro
  gROOT->LoadMacro("AddTask.C");
  //create an instance of your analysis task
  PionHadron *ptr = AddTask();




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
  //mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  //mgr->SetUserProgressBar(1,25);

  Bool_t isLocal = kFALSE;
  //if you want to run locally, we need to define some input
  if(isLocal){ 
    TChain *chain = new TChain("esdTree");
    chain->Add("LHC13d_ESD1.root");
    chain->Add("LHC13d_ESD2.root");
    //chain->Add("LHC13d.root");f
    std::cout << "Running locally on tree, Chain has " << chain->GetEntries() << " events " << std::endl;
     mgr->StartAnalysis("local",chain);
     return;
  }

    
  AliAnalysisAlien *plugin = new AliAnalysisAlien()
  plugin->SetOverwriteMode();
  plugin->SetAliPhysicsVersion("vAN-20170510-1");
  //plugin->SetAliPhysicsVersion("vAN-20160203-1");
  plugin->SetGridDataDir("/alice/data/2013/LHC13d");
  plugin->SetDataPattern("*/pass4/*/*ESDs.root");
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(195783);
  
  plugin->SetGridWorkingDir("workdir");
  plugin->SetGridOutputDir("output");
  
  //also specify the include (header) paths on grid
  plugin->AddIncludePath ("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  //make sure your source files get copied to the grid
  plugin->SetAdditionalLibs("PionHadron.cxx PionHadron.h");
  plugin->SetAnalysisSource("PionHadron.cxx");  
  
  ///alice/data/2013/LHC13d/000195682/pass4/AOD/
  //Run 2 pp simulation Dijet
  ///alice/sim/2016/LHC16h3/11/244617
  //plugin->SetGridDataDir("/alice/sim/2016/LHC16h3/1");

  plugin->SetSplitMaxInputFileNumber(40);
  plugin->SetExecutable("myTask.sh");
  //specify how many seconds your job may take
  plugin->SetTTL(10000);
  plugin->SetJDLName("PionHadron.jdl");
  
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetKeepLogs(kTRUE);
  plugin->SetMaxMergeStages(1);
  plugin->SetMergeViaJDL(kTRUE);

  mgr->SetGridHandler(plugin);

  const char *alien_close_se = gSystem->Getenv("alien_CLOSE_SE");

  if (alien_close_se != NULL) {
    const char *file = mgr->GetCommonFileName();

    plugin->SetDefaultOutputs(kFALSE);
    plugin->SetOutputFiles(Form(
				"%s@%s", file, alien_close_se));
    plugin->SetOutputArchive(Form(
"log_archive.zip:stdout,stderr@%s "
"root_archive.zip:%s,*.stat@%s",
alien_close_se, file, alien_close_se));
  }

  Bool_t isTesting = kFALSE;
  if(isTesting){
    std::cout<< "Running Grid Test " << std::endl;
    plugin->SetNtestFiles(1);
    plugin->SetRunMode("test"); // 
  }
  else{
    plugin->SetRunMode("full"); //full, merge, terminate
  }
  
  mgr->StartAnalysis("grid");
  return;
}

