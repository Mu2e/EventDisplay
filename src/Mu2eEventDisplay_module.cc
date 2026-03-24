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
#include <fstream>

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
          fhicl::Atom<bool> showCrv{Name("showCrv"), Comment("set false if you just want to see DS"),false};
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
          fhicl::Atom<bool> addCrvBars{Name("addCrvBars"), Comment("show crv bars hit if presrnt"),true};
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
          fhicl::Atom<std::string>geometry{Name("geometry"),Comment("geometry type (nominal_2020 or nominal_2025)"),"nominal_2025"};
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
        
        // Control between the main thread and event-display thread
        std::condition_variable cv_{};
        std::mutex m_{};

        int  diagLevel_;
        bool showCrv_;
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
        bool addCrvBars_;
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
        
        std::unique_ptr<GUI> fGui{nullptr};
        std::unique_ptr<TextSelect> fText{nullptr};
        std::unique_ptr<PrintInfo> fPrint{nullptr};
        REX::REveManager* eveMng_{nullptr};
        std::string geometry_;
        std::unique_ptr<EventDisplayManager> eventMgr_{nullptr};
        double eventid_;
        double runid_;
        double subrunid_;
        bool seqMode_;
        int eventn;
        int runn;
        int subrunn;
        int autoplay = 0;
        std::vector<std::shared_ptr<DataProduct>> listoflists;
        GeomOptions geomOpts;
        ConfigFileLookupPolicy configFile;
    };
}

  Mu2eEventDisplay::Mu2eEventDisplay(const Parameters& conf)  :
    art::EDAnalyzer(conf),
    diagLevel_(conf().diagLevel()),
    showCrv_(conf().showCrv()),
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
    addCrvBars_(conf().addCrvBars()),
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
    geometry_(conf().geometry()),
    seqMode_(conf().seqMode())
    {
      geomOpts.fill(showCrv_,showPS_, showTS_, showDS_, show2D_, caloVST_, showST_, extracted_, showSTM_, showCalo_, showTracker_, showCaloCrystals_, showEM_ );
    }

  Mu2eEventDisplay::~Mu2eEventDisplay() {}

  void Mu2eEventDisplay::signalAppStart()
  { 
      std::unique_lock lock{m_};
      cv_.notify_all();
  }

  void Mu2eEventDisplay::beginJob(){
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting ..."<<std::endl;
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- geometry_ = '" << geometry_ << "'"<<std::endl;
      
      // Map geometry type to GUIID
      int guiid = 4313;  // Default to 2025
      if (geometry_ == "nominal_2020") {
          guiid = 4336;
      } else if (geometry_ == "nominal_2025") {
          guiid = 4313;
      } else if (geometry_ == "extracted_2025") {
          guiid = 4285;
      } else if (geometry_ == "extracted_2020") {
          guiid = 4285;
      }
      
      // Write GUIID to config file
      std::ofstream configFile("EventDisplay/config/guiutils_current.txt");
      if (configFile.is_open()) {
          configFile << "// Automatically generated GUI configuration\n";
          configFile << "int GUIID = " << guiid << ";\n";
          configFile.close();
          std::cout << "[Mu2eEventDisplay::beginJob()] Wrote GUIID to EventDisplay/config/guiutils_current.txt: " << guiid << std::endl;
      } else {
          std::cerr << "[Mu2eEventDisplay::beginJob()] WARNING: Could not open EventDisplay/config/guiutils_current.txt for writing" << std::endl;
      }
        std::unique_lock lock{m_};
        appThread_ = std::thread{[this] { run_application(); }};
        // Wait for app init to finish ... this will process pending timer events.
        XThreadTimer sut([this]{ signalAppStart(); });
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting wait on app start"<<std::endl;
        cv_.wait(lock);
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- app start signal received, starting eve init"<<std::endl;
        XThreadTimer suet([this]{ setup_eve(); });
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- starting wait on eve setup"<<std::endl;
        cv_.wait(lock);
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : beginJob()] -- eve setup apparently complete"<<std::endl;

    }

  void Mu2eEventDisplay::beginRun(const art::Run&){}

  void Mu2eEventDisplay::printOpts(){
    // Option to print out input options
    std::cout<<"*********** REve Mu2e **************"
    <<" User Options: "
    <<" addHits : "<< filler_.addHits_
    <<" addBkgClusters : "<<filler_.addBkgClusters_
    <<" addTimeClusters : "<<filler_.addTimeClusters_
    <<" addCrvRecoPulse : "<<filler_.addCrvRecoPulse_
    <<" addCrvClusters : "<<filler_.addCrvClusters_
    <<" addClusters : "<<filler_.addClusters_
    <<" addCaloDigi : "<<filler_.addCaloDigis_
    <<" addHelices : "<<filler_.addHelixSeeds_
    <<" addTracks : "<<filler_.addKalSeeds_
    <<" addCosmicTrackSeeds : "<<filler_.addCosmicTrackSeeds_ << std::endl;
  }


template <class T, class S> 
void Mu2eEventDisplay::FillAnyCollection(const art::Event& evt, std::vector<std::shared_ptr<DataProduct>>& list, std::tuple<std::vector<std::string>, std::vector<S>>& tuple){
    
    // Use Art's getMany() to retrieve a vector of all art::Handle<T> objects 
    // that exist in the event record for the type T.
    // 
    std::vector<art::Handle<T>> vah = evt.getMany<T>();
    
    // Temporary variables to store the labels and the product pointers/handles.
    std::string name;
    std::vector<std::string> alabel; // Stores fully qualified names (label_instance_process)
    std::vector<S> alist;  // Stores the actual product pointers (const T*)
    
    // Loop over the list of art::Handle<T> retrieved from the event.
    for (auto const& ah : vah) {
        // Get the Provenance (metadata) for the product associated with this handle.
        // Provenance tells you which module created the data.
        const art::Provenance* prov = ah.provenance();

        // --- Extract Provenance Information ---
        
        // Friendly Class Name (e.g., "CaloClusterCollection")
        std::string fcn = prov->friendlyClassName();
        // Module Label (the name of the Art module that created it)
        std::string modn = prov->moduleLabel();
        // Process Name (the Art process name, used for configuration context)
        std::string instn = prov->processName();
        
        // Get the raw product pointer/handle and push it onto the list.
        // The type S is typically const T* or a similar type that references the product.
        alist.push_back(ah.product()); 
        
        // Create a unique, descriptive name string for this collection instance.
        std::string name = fcn + "_" + prov->moduleLabel() + "_" + instn;
        
        // Store the name so it can be displayed in the GUI alongside the collection.
        alabel.push_back(name);
        
        if(diagLevel_ == 1){
            std::cout<<"extracting name =  "<<fcn<<" "<<modn<<" "<<instn<<std::endl;
            // Diagnostic printout of the provenance object type ID (for debugging)
            std::cout<<"with type =  "<<typeid(prov).name()<<std::endl;
        }
    }
    
    // Assign the list of labels and the list of products to the output tuple.
    // This tuple is then used to populate the display structure.
    tuple = std::make_tuple(alabel,alist);
    }

  void Mu2eEventDisplay::analyze(art::Event const& event){
        printOpts();
      // Clear all previously stored event objects in the display data structure. 
      // This prepares the structure for the new event's data.
      data.Reset();

      // Update the internal state variables based on the current event object ('event').
      displayedEventID_ = event.id();
      eventid_ = event.id().event();
      runid_ = event.run();
      subrunid_ = event.subRun();

      // Temporary vector used during the collection filling process (often required by FillAnyCollection).
      std::vector<std::shared_ptr<DataProduct>> _chits;

      // User Input Handling (REve Command Interface) ---
      // The fText pointer (TextSelect object) holds the Run/Event numbers entered by the user
      // via the REve GUI command box.
      
      if (fText) {
          // Retrieve the user-specified Run/Event numbers from the TextSelect object.
          std::pair<int, int> user_input = fText->getRunEvent();
          int user_run = user_input.first;
          int user_event = user_input.second;
          autoplay = fText->getAutoplay();
          std::cout << "\n[Mu2eEventDisplay::analyze] -------------------------" << std::endl;
          std::cout << "[Mu2eEventDisplay::analyze] User Input Detected:" << std::endl;
          std::cout << "[Mu2eEventDisplay::analyze] Run Number:  " << user_run << std::endl;
          std::cout << "[Mu2eEventDisplay::analyze] Event Number: " << user_event << std::endl;
          std::cout << "[Mu2eEventDisplay::analyze] Autoplay set " << autoplay << std::endl;
          std::cout << "[Mu2eEventDisplay::analyze] -------------------------\n" << std::endl;
          
          // Check if valid input was provided (Run and Event are non-zero).
          if(user_run !=0 and user_event != 0){
              // Store the user-requested event number internally.
              runn = user_run;
              eventn = user_event;
              // Disable sequential processing mode. The module must stop after this event is found.
              seqMode_ = false;
          }
      }

      // Event Filtering Logic ---
      // Process the event ONLY IF:
      // 1. We are in sequential mode (seqMode_ is true, meaning iterate through all events), OR
      // 2. The current event matches the user-requested event (runn, eventn).
      if((seqMode_) or ( runid_ == runn and subrunid_ == subrunid_ and eventid_ == eventn)){
          
          // Acquire a lock on the mutex. This ensures the Art thread (analyze) waits 
          // for the REve thread (GUI) to finish processing the display before proceeding 
          // to the next event or releasing the lock.
          std::unique_lock lock{m_};
          
          if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- Fill collections "<<std::endl;
          
          // Filling Data Collections ---
          // The following blocks conditionally load various data collections (CaloClusters, ComboHits, MCTrajectories, etc.)
          // based on user configuration flags (filler_.add...).
          // FillRecoCollections and FillAnyCollection are helper methods that retrieve data products
          // from the event record and store them in the internal data structure for display.
          
          if(filler_.addClusters_) {
              // If specifyTag_ is true, use the dedicated filler method (likely based on Art module tag).
              if(specifyTag_) filler_.FillRecoCollections(event, data, CaloClusters);
              // Otherwise, use the generic collection filler (likely based on type).
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
          if(filler_.addBkgClusters_) {
              if(specifyTag_) { filler_.FillRecoCollections(event, data, BkgClusters); }
            else { FillAnyCollection<BkgClusterCollection, const BkgClusterCollection*>(event, _chits, data.bkgcluster_tuple ); }
          }
          if(filler_.addHelixSeeds_){
              if(specifyTag_) { filler_.FillRecoCollections(event, data, HelixSeeds); }
              else { FillAnyCollection<HelixSeedCollection, const HelixSeedCollection*>(event, _chits, data.helix_tuple ); }
          }

          if(filler_.addKalSeeds_) {
              if(specifyTag_) { filler_.FillRecoCollections(event, data, KalSeeds); }
              else { FillAnyCollection<KalSeedPtrCollection, const KalSeedPtrCollection*>(event, _chits, data.track_tuple ); }
          }

          // --- MC Collections (Monte Carlo) ---
          if(filler_.addMCTraj_) {
              if(specifyTag_) { filler_.FillMCCollections(event, data, MCTrajectories); }
              else { FillAnyCollection<MCTrajectoryCollection, const MCTrajectoryCollection*>(event, _chits, data.mctrack_tuple ); }
          }

          if(filler_.addSurfSteps_) {
              if(specifyTag_) { filler_.FillMCCollections(event, data, SurfaceSteps); }
              else { FillAnyCollection<SurfaceStepCollection, const SurfaceStepCollection*>(event, _chits, data.surfstep_tuple ); }
          }
          if(filler_.addSimParts_) {
              if(specifyTag_) { filler_.FillMCCollections(event, data, SimParticles); }
              else { FillAnyCollection<SimParticleCollection, const SimParticleCollection*>(event, _chits, data.sim_tuple ); }
          }
          // --- End MC Collections ---

          if(filler_.addTimeClusters_) {
              if(specifyTag_) { filler_.FillRecoCollections(event, data, TimeClusters);}
              else { FillAnyCollection<TimeClusterCollection, const TimeClusterCollection*>(event, _chits, data.timecluster_tuple );}
          }

          if(filler_.addCrvRecoPulse_) {
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

          // Thread Synchronization
          // Start the single-event processing job (updating the REve display) in the REve/ROOT thread.
          // The XThreadTimer ensures the display update happens outside the current Art/analysis thread.
          XThreadTimer proc_timer([this]{ process_single_event(); }); 

          if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- transferring to TApplication thread "<<std::endl;
          
          if (autoplay > 0) {
              // Autoplay is ON (Autoplay value is the delay in seconds)
              std::cout << "Auto play switched on.... waiting 10 s for REve display." << std::endl;
              
              auto timeout = std::chrono::seconds(5);
              cv_.wait_for(lock, timeout); 

              lock.unlock();
              
          } else {
              cv_.wait(lock);
          }
          
          // Reset these to return to sequential navigation
          seqMode_ = true;
          runn = 0;
          eventn = 0;
          

          if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] -- TApplication thread returning control "<<std::endl;
          if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : analyze()] Ended Event "<<std::endl;
          
      }
      
  }

    void Mu2eEventDisplay::endJob()
    {
        // Check if the diagnostic level is set to print informational messages.
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : EndJob] Start "<<std::endl;
        
        // --- 1. Terminate the ROOT Application ---
        // Terminate the TApplication instance. This is the main event loop 
        // running in the separate thread (appThread_). 
        // The argument '0' typically signals a normal shutdown.
        application_.Terminate(0);

        // --- 2. Join the Application Thread ---
        // Check if the thread (where the TApplication/REve display is running) 
        // is currently running and hasn't been joined yet.
        if (appThread_.joinable()) {
            // Wait for the application thread to finish its work and exit cleanly.
            // This is crucial for proper thread synchronization and preventing 
            // the main program from exiting before the thread resources are released.
            // 
            appThread_.join();
        }
        
        // Final cleanup message.
        if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : EndJob] End "<<std::endl;
    }


  void Mu2eEventDisplay::run_application()
  {
      // Without this the startup timer might not get invoked.
      // Explicitly process any waiting system events (like timer ticks or thread signals) 
      // that were queued before the TApplication event loop officially started. 
      // This ensures initialization tasks (like setting up the REve browser/GUI) are executed promptly.
      gSystem->ProcessEvents(); 

      // Start the TApplication event loop. 
      // The argument 'true' typically means that the function should return only when 
      // the application is explicitly terminated (e.g., via application_.Terminate(0) in endJob()).
      // This line blocks the appThread_ until the user closes the display or the Art job finishes.
      // 
      application_.Run(true);
  }


  void Mu2eEventDisplay::setup_eve()
  {
    // Assign the current thread as the main thread for RWebWindowsManager. 
    // This is critical for thread-safety and ensuring the display correctly manages events.
    RWebWindowsManager::AssignMainThrd(); 

    // Create the single instance of the REveManager. This is the core object 
    // that manages all display elements, scenes, and remote connections.
    eveMng_ = REX::REveManager::Create(); 

    // Configure REve to disallow multiple simultaneous browser connections.
    // The second 'false' is typically related to connection handling details.
    eveMng_->AllowMultipleRemoteConnections(false, false); 

    // Set the global WebWindowsManager to not use session keys for connections, 
    // simplifying access to the display.
    ROOT::RWebWindowsManager::SetUseSessionKey(false); 

    // --- Object Creation (Using std::unique_ptr for robust lifetime management) ---

    // Create the custom GUI object.
    fGui = std::make_unique<GUI>(); 
    fGui->SetName("Mu2eGUI");

    // Create the PrintInfo object for displaying data details.
    fPrint = std::make_unique<PrintInfo>(); 

    // Create the TextSelect object, used to capture user-input Run/Event numbers.
    fText = std::make_unique<TextSelect>(); 

    // --- Event Manager Instantiation and Singleton Link ---

    // Instantiate the core EventDisplayManager.
    // Geometry will be loaded from config file in the constructor.
    eventMgr_ = std::make_unique<EventDisplayManager>(
        eveMng_, 
        cv_, 
        m_, 
        fGui.get(),
        "2025"  // Not used; geometry loaded from config file in constructor
    );
    std::cout << "[Mu2eEventDisplay::setup_eve()] EventDisplayManager created and geometry loaded" << std::endl;

    //#eventMgr_->setGeometry("2025"); 
    // --- Scene Setup ---

    // Get the top-level scene, the "World," which is the container for all elements.
    auto world = eveMng_->GetWorld(); 
    assert(world); // Ensure the world was successfully retrieved.

    // Instantiate the custom MainWindow (likely handles the GDML geometry).
    frame_ = new MainWindow(); 
    // Load the geometry into the main window scene.
    frame_->makeGeometryScene(eveMng_, geomOpts, gdmlname_); 

    // --- Custom GUI and HTML Page Setup ---

    // Add a location path for custom files (like JavaScript or CSS).
    eveMng_->AddLocation("mydir/", configFile("EventDisplay/CustomGUIv2")); 
    // Set the default HTML page that the web browser will load.
    eveMng_->SetDefaultHtmlPage("file:mydir/eventDisplay.html"); 

    // --- Add Elements to the World (Element ID Assignment occurs here) ---

    // Elements must be added to a scene/world before their Element ID is assigned and valid.
    world->AddElement(fText.get());     
    world->AddElement(eventMgr_.get()); 
    world->AddElement(fPrint.get());    
    //world->AddElement(fGui.get());

    // --- Final Linking and Command Registration ---

    //#eventMgr_->setTextSelectId(fText->GetElementId()); 

    // Register commands that the GUI buttons will execute. The command is routed to the 
    // method on the specified element (eventMgr_.get() or fPrint.get()).
    world->AddCommand("QuitRoot", "sap-icon://log", eventMgr_.get(), "QuitRoot()");
    world->AddCommand("NextEvent", "sap-icon://step", eventMgr_.get(), "NextEvent()");
    world->AddCommand("PrintMCInfo", "sap-icon://step", fPrint.get(), "PrintMCInfo()");
    world->AddCommand("PrintRecoInfo", "sap-icon://step", fPrint.get(), "PrintRecoInfo()");
    //#eventMgr_->setid(fText->GetElementId() );
    // --- Signal Art Thread to Proceed ---

    // Acquire a lock on the mutex.
    std::unique_lock lock{m_}; 
    // Signal the waiting Art thread (in analyze() or beginRun()) to continue processing 
    // now that the REve display is fully initialized.
    cv_.notify_all();
    std::cout << "[DEBUG] TextSelect object ID assigned: " << fText->GetElementId() << std::endl; 
}
  
  

  // Actually interesting function responsible for drawing the current event
  void Mu2eEventDisplay::process_single_event()
  {
      
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] Start "<<std::endl;
      
      // --- 1. Disable Redrawing and Start Change Tracking ---
      
      // Temporarily disable the browser redraw to prevent the display from updating 
      // multiple times during the process, leading to flickering and slow performance.
      eveMng_->DisableRedraw(); 
      
      // Tell the World Scene to start tracking changes. All element additions/modifications 
      // are batched until EndAcceptingChanges() is called.
      eveMng_->GetWorld()->BeginAcceptingChanges(); 
      
      // Tell all scenes managed by the REveManager to start accepting batched changes.
      eveMng_->GetScenes()->AcceptChanges(true); 

      // --- 2. Update GUI and PrintInfo Data Structures ---

      // Transfer the core event identifiers to the GUI object for display.
      fGui->feventid = eventid_;
      fGui->fsubrunid = subrunid_;
      fGui->frunid = runid_;

      // Transfer the relevant collection data to the PrintInfo object. 
      // This allows the PrintInfo command to access and display details when triggered.
      fPrint->fcalocluster_tuple = data.calocluster_tuple;
      fPrint->fmctrack_tuple = data.mctrack_tuple;
      fPrint->ftrack_tuple = data.track_tuple;

      // Update the custom properties displayed in the REve GUI (e.g., in a sidebar table).
      fGui->StampObjProps(); 

      // --- 3. Prepare the Scene and Data Options ---
      
      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- extract event scene "<<std::endl;
      // Get the dedicated scene for event data (where tracks, hits, etc., will be drawn).
      REX::REveElement* scene = eveMng_->GetEventScene(); 

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- calls to data interface "<<std::endl;

      // Create a structure defining which data products should be drawn (based on module configuration).
      DrawOptions drawOpts(filler_.addCosmicTrackSeeds_, filler_.addHelixSeeds_, filler_.addKalSeeds_, filler_.addCaloDigis_, filler_.addClusters_, filler_.addHits_, filler_.addBkgClusters_, filler_.addCrvRecoPulse_, filler_.addCrvClusters_, filler_.addTimeClusters_, filler_.addTrkHits_, filler_.addMCTraj_, filler_.addSurfSteps_, filler_.addSimParts_, addErrBar_, addCrystalHits_, addCrvBars_);

      // Create a structure defining visualization options specific to Kinematic/Kalman fitting results.
      KinKalOptions KKOpts(addKalInter_, addTrkStrawHits_, addTrkCaloHits_);

      // --- 4. Draw the Event ---

      // Call the core function responsible for converting data collections into REve elements (lines, points, clusters) 
      // and adding them to the event scene.
      frame_->showEvents(eveMng_, scene, firstLoop_, firstLoopCalo_, data, drawOpts, particles_, strawdisplay_, geomOpts, KKOpts);

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] -- cluster added to scene "<<std::endl;

      // --- 5. Finalize Changes and Enable Redraw ---

      // Set the flag to false, indicating subsequent events are not the "first loop" (allowing optimization).
      firstLoop_ = false; 
      
      // Stop tracking changes in all scenes.
      eveMng_->GetScenes()->AcceptChanges(false); 
      
      // Finalize the batch of changes for the World scene. This triggers the update signal.
      eveMng_->GetWorld()->EndAcceptingChanges(); 
      
      // Re-enable the browser redraw. The browser now performs a single, optimized refresh 
      // using all the batched changes.
      eveMng_->EnableRedraw(); 

      if(diagLevel_ == 1) std::cout<<"[Mu2eEventDisplay : process_single_event] End "<<std::endl;
    }

  using mu2e::Mu2eEventDisplay;
  DEFINE_ART_MODULE(Mu2eEventDisplay)
