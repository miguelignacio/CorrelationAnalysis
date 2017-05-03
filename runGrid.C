void runGrid()
{
  //since we will compile a class, tell root where to look for headers
  std::cout <<"Proceesing lines to include $ROOTSYS/include, $ALICE_ROOT/include and $ALICE_PHYSICS/include" << std::endl;
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  //create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager();
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  //  inputHandler->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(inputHandler);
  //compile the class (locally)
  gROOT->LoadMacro("PionHadron.cxx++g");
  //load the addtask macro
  gROOT->LoadMacro("AddTask.C");
  //create an instance of your analysis task
  PionHadron *ptr = AddTask();

  if(!mgr->InitAnalysis()) return;
  // mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  // mgr->SetUserProgressBar(1,25);

  Bool_t isLocal = kTRUE;
  //if you want to run locally, we need to define some input
 
  if(isLocal){ 
    TChain *chain = new TChain("aodTree");
    //chain->Add("AliAOD.root");
    chain->Add("MC_AliAOD_EMCAL.root");
    //chain->Add("MC_MinBias.root");
  //star the analysis locally
    std::cout << "Running locally on AOD tree, Chain has " << chain->GetEntries() << " events " << std::endl;
     mgr->StartAnalysis("local",chain);
     return;
  }

  //All the rest is just in case you need to run in the grid  
  //create and configure the plugin
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  //also specify the include (header) paths on grid
  plugin->AddIncludePath ("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/incl\
ude -I$ALICE_PHYSICS/include");
  //make sure your source files get copied to the grid
  plugin->SetAdditionalLibs("PionHadron.cxx PionHadron.h");
  plugin->SetAnalysisSource("PionHadron.cxx");  
  //select the aliphysics version
  //AliPhysics/vAN-20170425-1
  plugin->SetAliPhysicsVersion("v5-08-14-01-1");
  //select the Alien API version
  plugin->SetAPIVersion("V1.1x");

  //13 TeV pp in 2016 run
  plugin->SetGridDataDir("/alice/data/2015/LHC15o");
  plugin->SetDataPattern("*/pass1/*/*AOD.root");
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(246113);

  plugin->SetSplitMaxInputFileNumber(40);
  plugin->SetExecutable("myTask.sh");
  //specify how many seconds your job may take
  plugin->SetTTL(10000);
  plugin->SetJDLName("PionHadron.jdl");
  
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetKeepLogs(kTRUE);
  plugin->SetMaxMergeStages(1);
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetGridWorkingDir("workdir");
  plugin->SetGridOutputDir("5-3-17");
  //connect the alien plugin to the manager
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

  Bool_t isTesting = kTRUE;
  if(isTesting){
    std::cout<< "Running Grid Test " << std::endl;
    plugin->SetNtestFiles(4);
    plugin->SetRunMode("test"); // 
  }
  else{
    plugin->SetRunMode("merge"); //full, merge, terminate
  }
  
  mgr->StartAnalysis("grid");
  return;
}

