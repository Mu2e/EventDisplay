#include <TObject.h>
#include <TSystem.h>
#include <TFile.h>
#include "REve/inc/CollectionFiller.hh"
#include "art/Framework/Principal/SubRun.h"

using namespace mu2e;
namespace mu2e{

    CollectionFiller::CollectionFiller(const Config& conf) :
      chTag_(conf.chTag()),
      tcTag_(conf.tcTag()),
      crvrecoTag_(conf.crvrecoTag()),
      crvcoinTag_(conf.crvcoinTag()),
      calodigTag_(conf.cluTag()),
      cluTag_(conf.cluTag()),
      helixSeedTag_(conf.helixSeedTag()),
      kalSeedTag_(conf.kalSeedTag()),
      cosmicTrackSeedTag_(conf.cosmicTrackSeedTag()),
      MCTrajTag_(conf.MCTrajTag()),
      addHits_(conf.addHits()),
      addCrvHits_(conf.addCrvHits()),
      addCrvClusters_(conf.addCrvClusters()),
      addTimeClusters_(conf.addTimeClusters()),
      addTrkHits_(conf.addTrkHits()),
      addCaloDigis_(conf.addCaloDigis()),
      addClusters_(conf.addClusters()),
      addHelixSeeds_(conf.addHelixSeeds()),
      addKalSeeds_(conf.addKalSeeds()),
      addCosmicTrackSeeds_(conf.addCosmicTrackSeeds()),
      addMCTraj_(conf.addMCTraj()),
      FillAll_(conf.FillAll())
      {}

    template <typename T> std::string TurnNameToString( const T& value )
    {
        std::ostringstream ss;
        ss << value;
        return ss.str();
    }
   
    void CollectionFiller::FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName CollectionName)  {
        if(FillAll_  or (CollectionName == ComboHits)){ 
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
        if(FillAll_ or (addCrvHits_ and CollectionName==CRVRecoPulses)){
          for(const auto &tag : crvrecoTag_){ 
            auto chH = evt.getValidHandle<mu2e::CrvRecoPulseCollection>(tag);
            data.crvrecocol = chH.product();
            data.crvpulse_list.push_back(data.crvrecocol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting CRV Instance: "<<name<<"  "<<data.crvpulse_list.size()<<std::endl;
            data.crvpulse_labels.push_back(name);
          }
            data.crvpulse_tuple = std::make_tuple(data.crvpulse_labels,data.crvpulse_list);
        }
        if(FillAll_ or (addCrvClusters_ and CollectionName==CRVCoincidenceCluster)){
          for(const auto &tag : crvcoinTag_){ 
            auto chH = evt.getValidHandle<mu2e::CrvCoincidenceClusterCollection>(tag);
            data.crvcoincol = chH.product();
            data.crvcoin_list.push_back(data.crvcoincol);
            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting CRV Instance: "<<name<<"  "<<data.crvcoin_list.size()<<std::endl;
            data.crvpulse_labels.push_back(name);
          }
            data.crvpulse_tuple = std::make_tuple(data.crvpulse_labels,data.crvpulse_list);
        }
        if(FillAll_  or (CollectionName == TimeClusters)){ 
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
        if(FillAll_  or (CollectionName == CaloDigis)){
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
        if(FillAll_  or (CollectionName == CaloClusters)){
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
        if(FillAll_ or (CollectionName==KalSeeds)){
          
            for(const auto &tag : kalSeedTag_){
                auto chH = evt.getValidHandle<mu2e::KalSeedCollection>(tag);
                data.kalSeedcol = chH.product();
                data.track_list.push_back(data.kalSeedcol);

                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting KalSeed Instance: "<<name<<std::endl;
                data.track_labels.push_back(name);

            }
            data.track_tuple = std::make_tuple(data.track_labels,data.track_list);
            
        }
        if(FillAll_  or (CollectionName == CosmicTrackSeeds)){
            auto chH = evt.getValidHandle<mu2e::CosmicTrackSeedCollection>(cosmicTrackSeedTag_);
            data.CosmicTrackSeedcol = chH.product();
        }
    }
    
    void CollectionFiller::FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName CollectionName){

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
        
    }

}
