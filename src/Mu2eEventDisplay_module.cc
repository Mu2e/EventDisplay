//Author: S Middleton
//Date: 2021
//Purpose: Mu2eEventDisplay Driving module

#include "TRACE/trace.h"
//Art:
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "artdaq-core/Data/ContainerFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
//#include "artdaq/DAQdata/Globals.hh"

#include "cetlib_except/exception.h"

//ROOT:
//#include "art_root_io/TFileService.h"
#include <TApplication.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>

//EVE-7
#include <ROOT/RWebWindow.hxx>
#include <ROOT/RWebWindowsManager.hxx>
#include <ROOT/REveManager.hxx>

#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>

//Offline:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic pop

//
#include "EventDisplay/inc/MainWindow.hh"
#include "EventDisplay/inc/EventDisplayManager.hh"
#include "EventDisplay/inc/CollectionFiller.hh"
#include "EventDisplay/inc/DataCollections.hh"
#include "EventDisplay/inc/GUI.hh"
#include "EventDisplay/inc/TextSelect.hh"
#include "EventDisplay/inc/DataProduct.hh"
#include "EventDisplay/inc/PrintInfo.hh"

//Ofline
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"


// Using declarations for code that might or might not be in ROOT::Experimental
#if NUMERIC_ROOT_VERSION>=6300400
using ROOT::RWebWindowsManager;
#else
using ROOT::Experimental::RWebWindowsManager;
#endif

using namespace std;

namespace mu2e
{

    class XThreadTimer : public TTimer {
        std::function<void()> foo_;
        public:
            XThreadTimer(std::function<void()> f) : foo_(f)
            {
                SetTime(0);
                R__LOCKGUARD2(gSystemMutex);
                gSystem->AddTimer(this);
            }
            Bool_t Notify() override
            {
                foo_();
                gSystem->RemoveTimer(this);
                return kTRUE;
            }
    };

    class Mu2eEventDisplay : public art::EDAnalyzer {
      public:
        struct Config{
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
          fhicl::Atom<bool> showCRV{Name("showCRV"), Comment("set false if you just want to see DS"),false};
          fhicl::Atom<bool> showPS{Name("showPS"), Comment("set false if you just want to see inside DS"),false};
          fhicl::Atom<bool> showTS{Name("showTS"), Comment("set false if you just want to see inside DS"),false};
          fhicl::Atom<bool> showDS{Name("showDS"), Comment("set false if you just want to see inside DS"),false};
          fhicl::Atom<bool> show2D{Name("show2D"), Comment(""),true};
          fhicl::Atom<bool> showST{Name("showST"), Comment(""),true};
          fhicl::Atom<bool> caloVST{Name("caloVST"), Comment(""),false};
          fhicl::Atom<bool> showSTM{Name("showSTM"), Comment(""),false};
          fhicl::Atom<bool> showCalo{Name("showCalo"), Comment(""),true};
          fhicl::Atom<bool> showTracker{Name("showTracker"), Comment(""),true};
          fhicl::Atom<bool> showCaloCrystals{Name("showCaloCrystals"), Comment(""),true};
          fhicl::Atom<bool> addErrBar{Name("addErrBar"), Comment("show combo hit err bar"),true};
          fhicl::Atom<bool> addCrystalHits{Name("addCrystalHits"), Comment("show crystal hits if presrnt"),true};
          fhicl::Atom<bool> addCRVBars{Name("addCRVBars"), Comment("show crv bars hit if presrnt"),true};
          fhicl::Atom<bool> addKalInter{Name("addKalInter"), Comment("show Kal intersections"),true};
          fhicl::Atom<bool> addTrkStrawHits{Name("addTrkStrawHits"), Comment("show Kal trk straw hits"),true};
          fhicl::Atom<bool> addTrkCaloHits{Name("addTrkCaloHits"), Comment("show Kal trk cal ohits"),true};
          fhicl::Atom<bool> specifyTag{Name("specifyTag"), Comment("to only select events of selected input tag"),false};
          fhicl::Table<CollectionFiller::Config> filler{Name("filler"),Comment("fill collections")};
          fhicl::Sequence<int>particles{Name("particles"),Comment("PDGcodes to plot")};
          fhicl::Atom<std::string>gdmlname{Name("gdmlname"),Comment("gdmlname")};
          fhicl::Atom<bool> strawdisplay{Name("strawdisplay"), Comment(""),true};
          fhicl::Atom<bool> extracted{Name("extracted"), Comment(""),false};
          fhicl::Atom<bool> showEM{Name("showEM"), Comment(""),false};
          fhicl::Atom<bool> seqMode{Name("seqMode"), Comment("turn off for go to any event functionality"),true};
        };

        typedef art::EDAnalyzer::Table<Config> Parameters;
        explicit Mu2eEventDisplay(const Parameters& conf);
        virtual ~Mu2eEventDisplay();
        virtual void beginJob() override;
        virtual void beginRun(const art::Run& run) override;
        virtual void analyze(const art::Event& e);
        virtual void endJob() override;
        void signalAppStart();
      private:

        art::ServiceHandle<art::TFileService> tfs;
        Config _conf;


        void setup_eve();
        void run_application();
        void process_single_event();
        void printOpts();
        template <class T, class S> void FillAnyCollection(const art::Event& evt, std::vector<std::shared_ptr<DataProduct>>& list, std::tuple<std::vector<std::string>, std::vector<S>>& tuple);

        // Application control
        TApplication application_{"none", nullptr, nullptr};
        std::thread appThread_{};

        // Display control
        art::EventID displayedEventID_{};
        REX::REveManager* eveMng_{nullptr};
        EventDisplayManager* eventMgr_{nullptr};

        // Control between the main thread and event-display thread
        std::condition_variable cv_{};
        std::mutex m_{};

        int  diagLevel_;
        bool showCRV_;
        bool showPS_;
        bool showTS_;
        bool showDS_;
        bool show2D_;
        bool showST_;
        bool caloVST_;
        bool showSTM_;
        bool showCalo_;
        bool showTracker_;
        bool showCaloCrystals_;
        bool addErrBar_;
        bool addCrystalHits_;
        bool addCRVBars_;
        bool addKalInter_;
        bool addTrkStrawHits_;
        bool addTrkCaloHits_;

        bool specifyTag_ = false;
        TDirectory*   directory_ = nullptr;
        CollectionFiller filler_;
        MainWindow *frame_;
        DataCollections data;
        bool firstLoop_ = true;
        bool firstLoopCalo_ = true;
        std::vector<int> particles_;
        std::string gdmlname_;
        bool strawdisplay_;
        bool extracted_;
        bool showEM_;

        // Setup Custom GUI
        GUI *fGui{nullptr};
        PrintInfo *fPrint{nullptr};
        TextSelect *fText{nullptr};
        double eventid_;
        double runid_;
        double subrunid_;
        bool seqMode_;
        int eventn;
        int runn;
        int subrunn;

        std::vector<std::shared_ptr<DataProduct>> listoflists;
        GeomOptions geomOpts;
        ConfigFileLookupPolicy configFile;
    };


  Mu2eEventDisplay::Mu2eEventDisplay(const Parameters& conf)  :
    art::EDAnalyzer(conf),
    diagLevel_(conf().diagLevel()),
    showCRV_(conf().showCRV()),
    showPS_(conf().showPS()),
    showTS_(conf().showTS()),
    showDS_(conf().showDS()),
    show2D_(conf().show2D()),
    showST_(conf().showST()),
    caloVST_(conf().caloVST()),
    showSTM_(conf().showSTM()),
    showCalo_(conf().showCalo()),
    showTracker_(conf().showTracker()),
    showCaloCrystals_(conf().showCaloCrystals()),
    addErrBar_(conf().addErrBar()),
    addCrystalHits_(conf().addCrystalHits()),
    addCRVBars_(conf().addCRVBars()),
    addKalInter_(conf().addKalInter()),
    addTrkStrawHits_(conf().addTrkStrawHits()),
    addTrkCaloHits_(conf().addTrkCaloHits()),
    specifyTag_(conf().specifyTag()),
    filler_(conf().filler()),
    particles_(conf().particles()),
    gdmlname_(configFile(conf().gdmlname())),
    strawdisplay_(conf().strawdisplay()),
    extracted_(conf().extracted()),
    showEM_(conf().showEM()),
    seqMode_(conf().seqMode())
    {
      std::cout<<"GDML file "<<gdmlname_<<std::endl;
       if(!seqMode_){
        // Take in Run, Event number
          std::cout<<" Event Number : "<<std::endl;
          cin>>eventn;
          std::cout<<" SubRun Number : "<<std::endl;
          cin>>subrunn;
          std::cout<<" Run Number : "<<std::endl;
          cin>>runn;
      }
      geomOpts.fill(showCRV_,showPS_, showTS_, showDS_, show2D_, caloVST_, showST_, extracted_, showSTM_, showCalo_, showTracker_, showCaloCrystals_, showEM_ );
    }

  Mu2eEventDisplay::~Mu2eEventDisplay() {}

  void Mu2eEventDisplay::signalAppStart()
  { 
      std::unique_lock lock{m_};
      cv_.notify_all();
  }

  void Mu2eEventDisplay::beginJob(){
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting ..."<<std::endl;
      {
            
      std::unique_lock lock{m_};

      appThread_ = std::thread{[this] { run_application(); }};

      // Wait for app init to finish ... this will process pending timer events.
      XThreadTimer sut([this]{ signalAppStart(); });
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting wait on app start"<<std::endl;
      cv_.wait(lock);
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- app start signal received, starting eve init"<<std::endl;
      //auto start1 = std::chrono::high_resolution_clock::now();

      XThreadTimer suet([this]{ setup_eve(); });
      //auto end1 = std::chrono::high_resolution_clock::now();
      //std::cout<<" time through process setup evene"<<std::chrono::duration<double, std::milli>(end1 - start1).count()<<" ms "<<std::endl;
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting wait on eve setup"<<std::endl;
      cv_.wait(lock);
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- eve setup apparently complete"<<std::endl;
      }
  }


  void Mu2eEventDisplay::beginRun(const art::Run&){

  }

  void Mu2eEventDisplay::printOpts(){
    std::cout<<"*********** REve Mu2e **************"
    <<" User Options: "
    <<" addHits : "<< filler_.addHits_
    <<" addTimeClusters : "<<filler_.addTimeClusters_
    <<" addCrvHits : "<<filler_.addCrvHits_
    <<" addCrvClusters : "<<filler_.addCrvClusters_
    <<" addClusters : "<<filler_.addClusters_
    <<" addHelices : "<<filler_.addHelixSeeds_
    <<" addTracks : "<<filler_.addKalSeeds_
    <<" addCosmicTrackSeeds : "<<filler_.addCosmicTrackSeeds_ << std::endl;
  }


  template <class T, class S> void Mu2eEventDisplay::FillAnyCollection(const art::Event& evt, std::vector<std::shared_ptr<DataProduct>>& list, std::tuple<std::vector<std::string>, std::vector<S>>& tuple){
      // get all instances of products of type T
      std::vector<art::Handle<T>> vah = evt.getMany<T>();
      std::string name;
      std::vector<std::string> alabel;
      std::vector<S> alist;
      // loop over the list of instances of products of this type
      for (auto const& ah : vah) {
        const art::Provenance* prov = ah.provenance();

        std::string fcn = prov->friendlyClassName();
        std::string modn = prov->moduleLabel();
        std::string instn = prov->processName();
        alist.push_back(ah.product());
        std::string name = fcn + "_" + prov->moduleLabel() + "_" + instn;
        alabel.push_back(name);
        if(diagLevel_ == 1){
          std::cout<<"extracting name =  "<<fcn<<" "<<modn<<" "<<instn<<std::endl;
          std::cout<<"with type =  "<<typeid(prov).name()<<std::endl;
        }
      }
      tuple = std::make_tuple(alabel,alist);

    }

  void Mu2eEventDisplay::analyze(art::Event const& event){

      //auto start = std::chrono::high_resolution_clock::now();

      //remove previous event objects;
      data.Reset();
      // Update state relevant for displaying new event.
      displayedEventID_ = event.id();
      eventid_ = event.id().event();
      runid_ = event.run();
      subrunid_ = event.subRun();

      std::vector<std::shared_ptr<DataProduct>> _chits;

      if((seqMode_) or ( runid_ == runn and subrunid_ == subrunn and eventid_ == eventn)){
        // Hand off control to display thread
        std::unique_lock lock{m_};
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- Fill collections "<<std::endl;
        //auto start1 = std::chrono::high_resolution_clock::now();
        // fill the collection lists
        if(filler_.addClusters_) {
          if(specifyTag_) filler_.FillRecoCollections(event, data, CaloClusters);
          else { FillAnyCollection<CaloClusterCollection, const CaloClusterCollection*>(event, _chits, data.calocluster_tuple);}
        }
        if(filler_.addCaloDigis_) {
          if(specifyTag_) filler_.FillRecoCollections(event, data, CaloDigis);
          else { FillAnyCollection<CaloDigiCollection, const CaloDigiCollection*>(event, _chits, data.calodigi_tuple);}
        }

        if(filler_.addHits_) {
          if(specifyTag_) { filler_.FillRecoCollections(event, data, ComboHits); }
          else { FillAnyCollection<ComboHitCollection, const ComboHitCollection*>(event, _chits, data.combohit_tuple ); }
        }

        if(filler_.addHelixSeeds_){
            if(specifyTag_) { filler_.FillRecoCollections(event, data, HelixSeeds); }
            else { FillAnyCollection<HelixSeedCollection, const HelixSeedCollection*>(event, _chits, data.helix_tuple ); }
        }

        if(filler_.addKalSeeds_) {
          if(specifyTag_) { filler_.FillRecoCollections(event, data, KalSeeds); }
          else { FillAnyCollection<KalSeedPtrCollection, const KalSeedPtrCollection*>(event, _chits, data.track_tuple ); }
        }

        if(filler_.addMCTraj_) {
          if(specifyTag_) { filler_.FillMCCollections(event, data, MCTrajectories); }
          else { FillAnyCollection<MCTrajectoryCollection, const MCTrajectoryCollection*>(event, _chits, data.mctrack_tuple ); }
        }

        if(filler_.addSurfSteps_) {
          if(specifyTag_) { filler_.FillMCCollections(event, data, SurfaceSteps); }
          else { FillAnyCollection<SurfaceStepCollection, const SurfaceStepCollection*>(event, _chits, data.surfstep_tuple ); }
        }

        if(filler_.addTimeClusters_) {
          if(specifyTag_) { filler_.FillRecoCollections(event, data, TimeClusters);}
          else { FillAnyCollection<TimeClusterCollection, const TimeClusterCollection*>(event, _chits, data.timecluster_tuple );}
        }

        if(filler_.addCrvHits_) {
          if(specifyTag_) { filler_.FillRecoCollections(event, data, CrvRecoPulses); }
          else { FillAnyCollection<CrvRecoPulseCollection, const CrvRecoPulseCollection*>(event, _chits, data.crvpulse_tuple );}
        }

        if(filler_.addCrvClusters_) {
          if(specifyTag_) { filler_.FillRecoCollections(event, data, CrvCoincidenceCluster); }
          else { FillAnyCollection<CrvCoincidenceClusterCollection, const CrvCoincidenceClusterCollection*>(event, _chits, data.crvcoin_tuple );}
        }

        if(filler_.addTrkHits_) filler_.FillRecoCollections(event, data, TrkHits);
        if(filler_.addCosmicTrackSeeds_)  filler_.FillRecoCollections(event, data, CosmicTrackSeeds);
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- Event processing started "<<std::endl;
        XThreadTimer proc_timer([this]{ process_single_event(); });

        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- transferring to TApplication thread "<<std::endl;
        cv_.wait(lock);
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- TApplication thread returning control "<<std::endl;
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] Ended Event "<<std::endl;
        seqMode_ = true;

        std::cout<<"test VALUE "<<eventMgr_->run<<std::endl;
      }
  }

    void Mu2eEventDisplay::endJob()
    {
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : EndJob] Start "<<std::endl;
      application_.Terminate(0);

      if (appThread_.joinable()) {
        appThread_.join();
      }
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : EndJob] End "<<std::endl;
    }



  // Functions invoked by the threads explicitly created by this module
  void Mu2eEventDisplay::run_application()
  {
      // Without this the startup timer might not get invoked.
      gSystem->ProcessEvents();
      application_.Run(true);
  }


  void Mu2eEventDisplay::setup_eve()
  {
      RWebWindowsManager::AssignMainThrd();
      eveMng_ = REX::REveManager::Create();
      eveMng_->AllowMultipleRemoteConnections(false, false);
      ROOT::RWebWindowsManager::SetUseSessionKey(false);
      //InitGuiInfo()
      fGui = new GUI();
      fGui->SetName("Mu2eGUI");

      fPrint = new PrintInfo();
      fText = new TextSelect();

      // call manager
      eventMgr_ = new EventDisplayManager{eveMng_, cv_, m_, fGui, fText};

      // access the world
      auto world = eveMng_->GetWorld();

      assert(world);

      frame_ = new MainWindow();
      frame_->makeGeometryScene(eveMng_, geomOpts, gdmlname_);

      //add path to the custom GUI code here, this overrides ROOT GUI
      eveMng_->AddLocation("mydir/", configFile("EventDisplay/CustomGUIv2"));
      eveMng_->SetDefaultHtmlPage("file:mydir/eventDisplay.html");

      // InitGuiInfo() cont'd
      world->AddElement(fGui);
      world->AddElement(fText);
      world->AddElement(eventMgr_);
      world->AddElement(fPrint);
      std::cout<<"[Mu2eEventDisplay::setup_eve] run in display is set to "<<eventMgr_->run<<std::endl;

      world->AddCommand("QuitRoot",  "sap-icon://log",  eventMgr_, "QuitRoot()");
      world->AddCommand("NextEvent", "sap-icon://step", eventMgr_, "NextEvent()");
      world->AddCommand("PrintMCInfo", "sap-icon://step", fPrint, "PrintMCInfo()");
      world->AddCommand("PrintRecoInfo", "sap-icon://step", fPrint, "PrintRecoInfo()");
      std::unique_lock lock{m_};
      cv_.notify_all();

  }

  // Actually interesting function responsible for drawing the current event
  void Mu2eEventDisplay::process_single_event()
    {
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] Start "<<std::endl;
      eveMng_->DisableRedraw();
      eveMng_->GetWorld()->BeginAcceptingChanges();
      eveMng_->GetScenes()->AcceptChanges(true);

      fGui->feventid = eventid_;
      fGui->fsubrunid = subrunid_;
      fGui->frunid = runid_;

      fPrint->fcalocluster_tuple = data.calocluster_tuple;
      fPrint->fmctrack_tuple = data.mctrack_tuple;
      fPrint->ftrack_tuple = data.track_tuple;

      std::cout<<"[Mu2eEventDisplay::process_single_event] display has run number set to "<<eventMgr_->run<<std::endl;
      std::cout<<"[Mu2eEventDisplay::process_single_event] value in the text class "<<fText->get()<<std::endl;

      fGui->StampObjProps();

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- extract event scene "<<std::endl;
      REX::REveElement* scene = eveMng_->GetEventScene();

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- calls to data interface "<<std::endl;

      // fill draw options
      DrawOptions drawOpts(filler_.addCosmicTrackSeeds_, filler_.addHelixSeeds_, filler_.addKalSeeds_, filler_.addCaloDigis_, filler_.addClusters_, filler_.addHits_,  filler_.addCrvHits_, filler_.addCrvClusters_, filler_.addTimeClusters_, filler_.addTrkHits_, filler_.addMCTraj_, filler_.addSurfSteps_, addErrBar_, addCrystalHits_, addCRVBars_);

      // fill kinkal options
      KinKalOptions KKOpts(addKalInter_, addTrkStrawHits_, addTrkCaloHits_);

      // call the "show events" function to add the
      frame_->showEvents(eveMng_, scene, firstLoop_, firstLoopCalo_, data, drawOpts, particles_, strawdisplay_, geomOpts, KKOpts);

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- cluster added to scene "<<std::endl;

      firstLoop_ = false;
      eveMng_->GetScenes()->AcceptChanges(false);
      eveMng_->GetWorld()->EndAcceptingChanges();
      eveMng_->EnableRedraw();

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] End "<<std::endl;
    }
  }
  using mu2e::Mu2eEventDisplay;
  DEFINE_ART_MODULE(Mu2eEventDisplay)
