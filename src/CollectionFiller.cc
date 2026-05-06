#include <TObject.h>
#include <TSystem.h>
#include <TFile.h>
#include "EventDisplay/inc/CollectionFiller.hh"
#include "art/Framework/Principal/SubRun.h"

using namespace mu2e;
namespace mu2e{

  CollectionFiller::CollectionFiller(const Config& conf) :
    chTag_(conf.chTag()),
    bcTag_(conf.bcTag()),
    tcTag_(conf.tcTag()),
    crvrecoTag_(conf.crvrecoTag()),
    crvcoinTag_(conf.crvcoinTag()),
    calodigTag_(conf.calodigTag()),
    cluTag_(conf.cluTag()),
    helixSeedTag_(conf.helixSeedTag()),
    kalSeedTag_(conf.kalSeedTag()),
    cosmicTrackSeedTag_(conf.cosmicTrackSeedTag()),
    MCTrajTag_(conf.MCTrajTag()),
    SurfStepsTag_(conf.SurfStepsTag()),
    SimTag_(conf.SimTag()),
    addHits_(conf.addHits()),
    addBkgClusters_(conf.addBkgClusters()),
    addCrvRecoPulse_(conf.addCrvRecoPulse()),
    addCrvClusters_(conf.addCrvClusters()),
    addCrvTrack_(conf.addCrvTrack()),
    addTimeClusters_(conf.addTimeClusters()),
    addTrkHits_(conf.addTrkHits()),
    addCaloDigis_(conf.addCaloDigis()),
    addClusters_(conf.addClusters()),
    addHelixSeeds_(conf.addHelixSeeds()),
    addKalSeeds_(conf.addKalSeeds()),
    addCosmicTrackSeeds_(conf.addCosmicTrackSeeds()),
    addMCTraj_(conf.addMCTraj()),
    addSurfSteps_(conf.addSurfSteps()),
    addSimParts_(conf.addSimParts()),
    addTrackerHist_(conf.addTrackerHist()),
    FillAll_(conf.FillAll())
  {}

  /**
 * @brief Converts any object or variable that can be streamed (via operator<<) 
 * into a standard C++ string.
 *
 * This function uses a std::ostringstream to perform the conversion in a type-safe 
 * and generic manner, relying on the type T's definition of operator<<.
 *
 * @tparam T The type of the value to be converted (e.g., int, double, a custom class).
 * @param value The value of type T to convert.
 * @return std::string The string representation of the input value.
 */
  template <typename T> 
  std::string TurnNameToString( const T& value )
  {
      // Create a std::ostringstream object. 
      std::ostringstream ss;

      // Insert the input 'value' into the string stream.
      ss << value;

      // Extract the final string from the string stream buffer and return it.
      return ss.str();
  }

  /**
 * @brief Retrieves specified reconstructed data products from the Art event.
 * * It iterates over a list of Art tags (labels) for each product type to handle multiple
 * instances of the same collection produced by different modules in the job.
 * * @param evt The current art::Event object.
 * @param data The DataCollections structure where pointers and labels are stored.
 * @param CollectionName An enum value used to selectively fill one type of collection.
 */
void CollectionFiller::FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName CollectionName) {
    
    // --- 1. ComboHits (Reconstructed Straw Tracker Hits) ---
    if(FillAll_ or (CollectionName == ComboHits)){
        for(const auto &tag : chTag_){
            auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(tag);
            
            data.chcol = chH.product();
            data.combohit_list.push_back(data.chcol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting ComboHit Instance: "<<name<<std::endl;
            data.combohit_labels.push_back(name);
        }
        data.combohit_tuple = std::make_tuple(data.combohit_labels,data.combohit_list);
    }

    // --- 2. BkgClusters (Background hit clusters) ---
    if(FillAll_  or (CollectionName == BkgClusters)){
      for(const auto &tag : bcTag_){
        auto bcH = evt.getValidHandle<mu2e::BkgClusterCollection>(tag);
        data.bccol = bcH.product();
        data.bkgcluster_list.push_back(data.bccol);
        std::string name = TurnNameToString(tag);
        std::cout<<"Plotting BkgCluster Instance: "<<name<<std::endl;
        data.bkgcluster_labels.push_back(name);
      }
      data.bkgcluster_tuple = std::make_tuple(data.bkgcluster_labels,data.bkgcluster_list);
    }

    // --- 3. CrvRecoPulses (Crv Raw Reconstructed Hits) ---
    if(FillAll_ or (addCrvRecoPulse_ and CollectionName==CrvRecoPulses)){
        for(const auto &tag : crvrecoTag_){
            auto chH = evt.getValidHandle<mu2e::CrvRecoPulseCollection>(tag);
            data.crvrecocol = chH.product();
            data.crvpulse_list.push_back(data.crvrecocol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting Crv Instance: "<<name<<"  "<<data.crvpulse_list.size()<<std::endl;
            data.crvpulse_labels.push_back(name);
        }
        data.crvpulse_tuple = std::make_tuple(data.crvpulse_labels,data.crvpulse_list);
    }

    
    // --- 4. CrvCoincidenceClusters (Crv Coincidence Clusters) ---
    if(FillAll_ or (addCrvClusters_ and CollectionName==CrvCoincidenceCluster)){
        std::cout << "Fill CrvClusters " << std::endl;
        for(const auto &tag : crvcoinTag_){
            std::cout << "Searching for tag " << tag << std::endl;
            auto chH = evt.getValidHandle<mu2e::CrvCoincidenceClusterCollection>(tag);
            data.crvcoincol = chH.product();
            data.crvcoin_list.push_back(data.crvcoincol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting Crv Instance: "<<name<<"  "<<data.crvcoin_list.size()<<std::endl;
            data.crvcoin_labels.push_back(name);
        }
        data.crvcoin_tuple = std::make_tuple(data.crvcoin_labels,data.crvcoin_list);
    }
    
    // --- 5. TimeClusters (Seed Time Clusters) ---
    if(FillAll_ or (CollectionName == TimeClusters)){
        for(const auto &tag : tcTag_){
            auto chH = evt.getValidHandle<mu2e::TimeClusterCollection>(tag);
            data.tccol = chH.product();
            data.timecluster_list.push_back(data.tccol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting TimeCluster Instance: "<<name<<std::endl;
            data.timecluster_labels.push_back(name);
        }
        data.timecluster_tuple = std::make_tuple(data.timecluster_labels,data.timecluster_list);
    }
    
    // --- 6. TrkHits (ComboHits treated as Tracking Hits) ---
    //FIXME - is this the same as the ComboHits version?
    if(FillAll_ or (addTrkHits_ and CollectionName == TrkHits)){
        for(const auto &tag : chTag_){
            auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(tag);
            data.chcol = chH.product();
            data.combohit_list.push_back(data.chcol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting TrkHit Instance: "<<name<<std::endl;
            data.combohit_labels.push_back(name);
        }
        data.combohit_tuple = std::make_tuple(data.combohit_labels,data.combohit_list);
    }
    
    // --- 7. CaloDigis (Raw Calorimeter Signal Samples) ---
    if(FillAll_ or (CollectionName == CaloDigis)){
        for(const auto &tag : calodigTag_){
            auto chH = evt.getValidHandle<mu2e::CaloDigiCollection>(tag);
            data.calodigicol = chH.product();
            data.calodigi_list.push_back(data.calodigicol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting CaloDigi Instance: "<<name<<std::endl;
            data.calodigi_labels.push_back(name);
        }
        data.calocluster_tuple = std::make_tuple(data.calocluster_labels,data.calocluster_list);
    }
    
    // --- 8. CaloClusters (Reconstructed Calorimeter Energy Clusters) ---
    if(FillAll_ or (CollectionName == CaloClusters)){
        for(const auto &tag : cluTag_){
            auto chH = evt.getValidHandle<mu2e::CaloClusterCollection>(tag);
            data.clustercol = chH.product();
            data.calocluster_list.push_back(data.clustercol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting CaloCluster Instance: "<<name<<std::endl;
            data.calocluster_labels.push_back(name);
        }
        data.calocluster_tuple = std::make_tuple(data.calocluster_labels,data.calocluster_list);
    }
    
    // --- 9. HelixSeeds (Pre-tracking Helix Fits) ---
    if(FillAll_ or (CollectionName==HelixSeeds)){
        for(const auto &tag : helixSeedTag_){
            auto chH = evt.getValidHandle<mu2e::HelixSeedCollection>(tag);
            data.helixSeedcol = chH.product();
            data.helix_list.push_back(data.helixSeedcol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting HelixSeed Instance: "<<name<<std::endl;
            data.helix_labels.push_back(name);
        }
        data.helix_tuple = std::make_tuple(data.helix_labels,data.helix_list);
    }
    
    // --- 10. KalSeeds (Final Kalman Filter Tracks) ---
    if(FillAll_ or (CollectionName==KalSeeds)){
        for(const auto &tag : kalSeedTag_){
            auto chH = evt.getValidHandle<mu2e::KalSeedPtrCollection>(tag); 
            data.kalSeedcol = chH.product();
            data.track_list.push_back(data.kalSeedcol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting KalSeed Instance: "<<name<<std::endl;
            data.track_labels.push_back(name);
        }
        data.track_tuple = std::make_tuple(data.track_labels,data.track_list);
    }
    
    // --- 11. CosmicTrackSeeds (Cosmic Ray Track Candidates) --- // FIXME - is this needed?
    if(FillAll_ or (CollectionName == CosmicTrackSeeds)){
        auto chH = evt.getValidHandle<mu2e::CosmicTrackSeedCollection>(cosmicTrackSeedTag_);
        data.CosmicTrackSeedcol = chH.product();
    }
}

/**
 * @brief Retrieves specified Monte Carlo data products from the Art event.
 * * It supports retrieving multiple instances of the same collection produced by
 * different simulation stages or modules by iterating over Art tags.
 * @param evt The current art::Event object.
 * @param data The DataCollections structure where pointers and labels are stored.
 * @param CollectionName An enum value used to selectively fill one type of MC collection.
 */
void CollectionFiller::FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName CollectionName){

    // --- 1. MCTrajectories (Full Monte Carlo Particle Paths) ---
    if(FillAll_ or (CollectionName==MCTrajectories)){

        for(const auto &tag : MCTrajTag_){
            auto chH = evt.getValidHandle<mu2e::MCTrajectoryCollection>(tag);
            
            data.mctrajcol = chH.product();
            data.mctrack_list.push_back(data.mctrajcol);

            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting MCTrajectory Instance: "<<name<<std::endl;
            data.mctrack_labels.push_back(name);
        }
        data.mctrack_tuple = std::make_tuple(data.mctrack_labels,data.mctrack_list);
    }

    // --- 2. SurfaceSteps (Steps taken by particles near critical surfaces/boundaries) ---
    if(FillAll_ or (CollectionName==SurfaceSteps)){

        for(const auto &tag : SurfStepsTag_){
            auto chH = evt.getValidHandle<mu2e::SurfaceStepCollection>(tag);
            
            data.surfstepcol = chH.product();
            data.surfstep_list.push_back(data.surfstepcol);

            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting SurfaceStep Instance: "<<name<<std::endl;
            data.surfstep_labels.push_back(name);
        }
        data.surfstep_tuple = std::make_tuple(data.surfstep_labels,data.surfstep_list);
    }
    
    // --- 3. SimParticles (Primary information about each simulated particle) ---
    if(FillAll_ or (CollectionName==SimParticles)){
        for(const auto &tag : SimTag_){
            auto chH = evt.getValidHandle<mu2e::SimParticleCollection>(tag);
            
            data.simcol = chH.product();
            data.sim_list.push_back(data.simcol);

            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting SimParticle Instance: "<<name<<std::endl;
            data.sim_labels.push_back(name);
        }
        data.sim_tuple = std::make_tuple(data.sim_labels,data.sim_list);
    }
}

}
