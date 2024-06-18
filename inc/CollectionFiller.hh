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

#include "REve/inc/DataCollections.hh"

namespace mu2e{

 enum RecoDataProductName {ComboHits,  CRVRecoPulses, TimeClusters, CaloClusters, HelixSeeds, KalSeeds, CosmicTrackSeeds, TrkHits, CRVCoincidenceCluster, CaloDigis};
 enum MCDataProductName {MCTrajectories};
	class CollectionFiller
	{
      public:
        struct Config{
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
          fhicl::Sequence<art::InputTag>chTag{Name("ComboHitCollection"),Comment("chTag")};
          fhicl::Sequence<art::InputTag>tcTag{Name("TimeClusterCollection"),Comment("ttcTag")};
          fhicl::Sequence<art::InputTag>crvrecoTag{Name("CrvRecoPulseCollection"),Comment("crvTag")};
          fhicl::Sequence<art::InputTag>crvcoinTag{Name("CrvCoincidenceClusterCollection"),Comment("crvcoinTag")};
          fhicl::Sequence<art::InputTag>calodigTag{Name("CaloDigiCollection"),Comment("caldigTag")};
          fhicl::Sequence<art::InputTag>cluTag{Name("CaloClusterCollection"),Comment("cluTag")};
          fhicl::Sequence<art::InputTag>helixSeedTag{Name("HelixSeedCollection"),Comment("helixseedTag")};
          fhicl::Sequence<art::InputTag>kalSeedTag{Name("KalSeedCollection"),Comment("kalseedTag")};
          fhicl::Atom<art::InputTag>cosmicTrackSeedTag{Name("CosmicTrackSeedCollection"),Comment("cosmicTrackSeedTag")};
          fhicl::Sequence<art::InputTag>MCTrajTag{Name("MCTrajectoryCollection"),Comment("MCTrajTag")};
          fhicl::Atom<bool> addHits{Name("addHits"), Comment("set to add the hits"),false};
          fhicl::Atom<bool> addCrvHits{Name("addCrvHits"), Comment("set to add crv hits"),false};
          fhicl::Atom<bool> addCrvClusters{Name("addCrvClusters"), Comment("set to add crv clusters"),false};
          fhicl::Atom<bool> addTimeClusters{Name("addTimeClusters"), Comment("set to add the CRV hits"),false};
          fhicl::Atom<bool> addTrkHits{Name("addTrkHits"), Comment("set to add the Trk hits"),false};
          fhicl::Atom<bool> addCaloDigis{Name("addCaloDigis"), Comment("set to add calodigis"),false};
          fhicl::Atom<bool> addClusters{Name("addClusters"), Comment("set to add caloclusters"),false};
          fhicl::Atom<bool> addHelixSeeds{Name("addHelixSeeds"), Comment("set to add helixseeds"),false};
          fhicl::Atom<bool> addKalSeeds{Name("addKalSeeds"), Comment("set to add kalseeds"),false};
          fhicl::Atom<bool> addCosmicTrackSeeds{Name("addCosmicTrackSeeds"), Comment("set to add cosmic track seeds"),false};
          fhicl::Atom<bool> addMCTraj{Name("addMCTraj"), Comment("set to add add MC information"),false};
          fhicl::Atom<bool> FillAll{Name("FillAll"), Comment("to see all available products"), false};
        };

        explicit CollectionFiller(const Config& conf);
        CollectionFiller(const CollectionFiller &);
        CollectionFiller& operator=(const CollectionFiller &);

        std::vector<art::InputTag> chTag_;
        std::vector<art::InputTag> tcTag_; 
        std::vector<art::InputTag> crvrecoTag_;
        std::vector<art::InputTag> crvcoinTag_;
        std::vector<art::InputTag> calodigTag_;
        std::vector<art::InputTag> cluTag_;
        std::vector<art::InputTag> helixSeedTag_;
        std::vector<art::InputTag> kalSeedTag_;
        art::InputTag cosmicTrackSeedTag_;
        std::vector<art::InputTag> MCTrajTag_;
        art::Event *_event;
        art::Run *_run;
	bool addHits_,  addCrvHits_, addCrvClusters_, addTimeClusters_, addTrkHits_, addCaloDigis_, addClusters_, addHelixSeeds_, addKalSeeds_, addCosmicTrackSeeds_, addMCTraj_, FillAll_; 
	      void FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName code);
        void FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName code);
        //template <class L> void test(const L& item);
        //template <class T, class V> int FillAllCollections(const art::Event& evt, std::vector<std::shared_ptr<V>>& list);
        virtual ~CollectionFiller(){};

      private:
        Config _conf;
      
};

}

#endif 
