#ifndef CollectionFiller_hh
#define CollectionFiller_hh

//Art:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include <TObject.h>
#include <TROOT.h>
#include <TGComboBox.h>
#include <TGListBox.h>

#include "EventDisplay/inc/DataCollections.hh"
namespace mu2e{

    // --- Enumerations for Data Product Selection ---
    
    /**
     * @brief Enum listing all reconstructed data products supported by the display.
     * Used as a code in FillRecoCollections to retrieve a specific product type.
     */
    enum RecoDataProductName {
        ComboHits, BkgClusters, CrvRecoPulses, TimeClusters, CaloClusters, HelixSeeds, 
        KalSeeds, CosmicTrackSeeds, TrkHits, CrvCoincidenceCluster, CaloDigis
    };
    
    /**
     * @brief Enum listing all Monte Carlo data products supported by the display.
     * Used as a code in FillMCCollections.
     */
    enum MCDataProductName {
        MCTrajectories, SurfaceSteps, SimParticles
    };
    
    /**
     * @brief Class responsible for configuring data product retrieval from the Art event.
     */
    class CollectionFiller
    {
    public:
        // --- FHiCL Configuration Struct ---
        /**
         * @brief Defines the structure of the configuration parameters for CollectionFiller
         * as read from the FHiCL file.
         */
        struct Config{
            using Name=fhicl::Name;
            using Comment=fhicl::Comment;
            
            // Debugging level
            fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
            
            // FHiCL Sequences for product tags (Art InputTags) - multiple instances allowed
            fhicl::Sequence<art::InputTag>chTag{Name("ComboHitCollection"),Comment("chTag")};
            fhicl::Sequence<art::InputTag>bcTag{Name("BkgClusterCollection"),Comment("bcTag")};
            fhicl::Sequence<art::InputTag>tcTag{Name("TimeClusterCollection"),Comment("ttcTag")};
            fhicl::Sequence<art::InputTag>crvrecoTag{Name("CrvRecoPulseCollection"),Comment("crvTag")};
            fhicl::Sequence<art::InputTag>crvcoinTag{Name("CrvCoincidenceClusterCollection"),Comment("crvcoinTag")};
            fhicl::Sequence<art::InputTag>calodigTag{Name("CaloDigiCollection"),Comment("calodigTag")};
            fhicl::Sequence<art::InputTag>cluTag{Name("CaloClusterCollection"),Comment("cluTag")};
            fhicl::Sequence<art::InputTag>helixSeedTag{Name("HelixSeedCollection"),Comment("helixseedTag")};
            fhicl::Sequence<art::InputTag>kalSeedTag{Name("KalSeedPtrCollection"),Comment("kalseedTag")};
            
            // FHiCL Atom for a single product tag (one instance expected)
            fhicl::Atom<art::InputTag>cosmicTrackSeedTag{Name("CosmicTrackSeedCollection"),Comment("cosmicTrackSeedTag")};
            
            // MC product tags
            fhicl::Sequence<art::InputTag>MCTrajTag{Name("MCTrajectoryCollection"),Comment("MCTrajTag")};
            fhicl::Sequence<art::InputTag>SurfStepsTag{Name("SurfaceStepCollection"),Comment("SurfaceSteps Collection Tag")};
            fhicl::Sequence<art::InputTag>SimTag{Name("SimParticleCollection"),Comment("SimTag")};
            
            // Boolean flags to enable/disable loading of specific collection types
            fhicl::Atom<bool> addHits{Name("addHits"), Comment("set to add the hits"),false}; // Corresponds to ComboHits
            fhicl::Atom<bool> addBkgClusters{Name("addBkgClusters"), Comment("set to add the bkg clusters"),false}; // Corresponds to bkg clusters
            fhicl::Atom<bool> addCrvRecoPulse{Name("addCrvRecoPulse"), Comment("set to add crv hits"),false}; // Corresponds to CrvRecoPulses
            fhicl::Atom<bool> addCrvClusters{Name("addCrvClusters"), Comment("set to add crv clusters"),false}; // Corresponds to CrvCoincidenceCluster
            fhicl::Atom<bool> addCrvTrack{Name("addCrvTrack"), Comment("set to add crv track intersections"),false};
            fhicl::Atom<bool> addTimeClusters{Name("addTimeClusters"), Comment("set to add the Crv hits"),false};
            fhicl::Atom<bool> addTrkHits{Name("addTrkHits"), Comment("set to add the Trk hits"),false}; // Alias for ComboHits/TrkHits
            fhicl::Atom<bool> addCaloDigis{Name("addCaloDigis"), Comment("set to add calodigis"),false};
            fhicl::Atom<bool> addClusters{Name("addClusters"), Comment("set to add caloclusters"),false};
            fhicl::Atom<bool> addHelixSeeds{Name("addHelixSeeds"), Comment("set to add helixseeds"),false};
            fhicl::Atom<bool> addKalSeeds{Name("addKalSeeds"), Comment("set to add kalseeds"),false};
            fhicl::Atom<bool> addCosmicTrackSeeds{Name("addCosmicTrackSeeds"), Comment("set to add cosmic track seeds"),false};
            fhicl::Atom<bool> addMCTraj{Name("addMCTraj"), Comment("set to add MCTrajectories"),false};
            fhicl::Atom<bool> addSurfSteps{Name("addSurfSteps"), Comment("set to add SurfaceStep MC"),false};
            fhicl::Atom<bool> addSimParts{Name("addSimParts"), Comment("set to add SimParticles MC"),false};
            fhicl::Atom<bool> addTrackerHist{Name("addTrackerHist"), Comment("set to add tracker histogram"),false};
            fhicl::Atom<bool> addCaloHist{Name("addCaloHist"), Comment("set to add calorimeter histogram"),false};
            
            // Global flag to attempt retrieval of all collections found in the event
            fhicl::Atom<bool> FillAll{Name("FillAll"), Comment("to see all available products"), false};
        };

        // --- Constructors and Operators ---
        
        /**
         * @brief Primary constructor; accepts and initializes member variables from the FHiCL Config struct.
         * @param conf The FHiCL configuration object.
         */
        explicit CollectionFiller(const Config& conf);
        
        // Disable copy constructor and assignment operator to prevent issues with Art pointers/handles
        CollectionFiller(const CollectionFiller &);
        CollectionFiller& operator=(const CollectionFiller &);

        // --- Input Tag Members (Filled from FHiCL config) ---
        // These store the actual InputTags used by art::Event::getValidHandle().
        std::vector<art::InputTag> chTag_;
        std::vector<art::InputTag> bcTag_;
        std::vector<art::InputTag> tcTag_;
        std::vector<art::InputTag> crvrecoTag_;
        std::vector<art::InputTag> crvcoinTag_;
        std::vector<art::InputTag> calodigTag_;
        std::vector<art::InputTag> cluTag_;
        std::vector<art::InputTag> helixSeedTag_;
        std::vector<art::InputTag> kalSeedTag_;
        art::InputTag cosmicTrackSeedTag_;
        std::vector<art::InputTag> MCTrajTag_;
        std::vector<art::InputTag> SurfStepsTag_;
        std::vector<art::InputTag> SimTag_;
        
        // --- Event/Run Pointers (Optional: may be set externally for easy access) ---
        art::Event *_event;
        art::Run *_run;
        
        // --- Boolean Control Flags (Copied from FHiCL Config) ---
         bool addHits_, addBkgClusters_, addCrvRecoPulse_, addCrvClusters_, addCrvTrack_, addTimeClusters_, addTrkHits_, addCaloDigis_, 
              addClusters_, addHelixSeeds_, addKalSeeds_, addCosmicTrackSeeds_, addMCTraj_,
               addSurfSteps_, addSimParts_, addTrackerHist_, addCaloHist_, FillAll_;

             
        // --- Collection Retrieval Methods ---
        
        /**
         * @brief Fills the DataCollections structure with requested reconstructed data products.
         * @param evt The current Art event.
         * @param data The structure to populate with product pointers and labels.
         * @param code The specific RecoDataProductName to retrieve (if FillAll_ is false).
         */
        void FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName code);
        
        /**
         * @brief Fills the DataCollections structure with requested Monte Carlo data products.
         * @param evt The current Art event.
         * @param data The structure to populate with product pointers and labels.
         * @param code The specific MCDataProductName to retrieve (if FillAll_ is false).
         */
        void FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName code);
        
        // Default virtual destructor
        virtual ~CollectionFiller(){};

    private:
        // Copy of the FHiCL configuration object (if needed internally after construction)
        Config _conf;

    };
}

#endif
