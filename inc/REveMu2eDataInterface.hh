#ifndef REveMu2eDataInterface_hh
#define REveMu2eDataInterface_hh
#include <ROOT/REveElement.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveLine.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveBox.hxx>
#include <ROOT/REveScene.hxx>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <ROOT/REveTrackPropagator.hxx>
#include <ROOT/REveTrans.hxx>
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace mu2e;
namespace REX = ROOT::Experimental;
namespace mu2e{
    class REveMu2eDataInterface {
        public:
          static int const mstyle = 1;
          static int const msize = 5;
          explicit REveMu2eDataInterface(){};
          explicit REveMu2eDataInterface(const REveMu2eDataInterface &);
          REveMu2eDataInterface& operator=(const REveMu2eDataInterface &);
          virtual ~REveMu2eDataInterface() = default;
          #ifndef __CINT__

          inline constexpr double pointmmTocm(double mm){ return mm/10; };
          void AddComboHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple, REX::REveElement* &scene, bool strawdisplay, bool AddErrBar);
          void AddTimeClusters(REX::REveManager *&eveMng, bool firstloop, std::tuple<std::vector<std::string>, std::vector<const TimeClusterCollection*>> timecluster_tuple, REX::REveElement* &scene);
          void AddTrkHits(REX::REveManager *&eveMng, bool firstLoop_,  std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple, std::tuple<std::vector<std::string>, std::vector<const KalSeedPtrCollection*>> track_tuple, REX::REveElement* &scene);
          void AddCaloDigis(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloDigiCollection*>> calodigi_tuple, REX::REveElement* &scene, bool addCrystalDraw);
          void AddCaloClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloClusterCollection*>> calocluster_tuple, REX::REveElement* &scene, bool addCrystalDraw, double t1, double t2);
          void AddCRVInfo(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvRecoPulseCollection*>> crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars);
          void AddCRVClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvCoincidenceClusterCollection*>>  crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars);
          void AddHelixSeedCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const HelixSeedCollection*>> helix_tuple, REX::REveElement* &scene);
          void AddKalIntersection(KalSeed const& kalseed, REX::REveElement* &scene);
          template<class KTRAJc> void AddTrkStrawHit(KalSeed const& kalseed, REX::REveElement* &scene, std::unique_ptr<KTRAJc> &trajectory);
          void AddTrkCaloHit(KalSeed const& kalseed, REX::REveElement* &scene);
          void AddKalSeedPtrCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const KalSeedPtrCollection*>> track_tuple, REX::REveElement* &scene );
          template<class KTRAJ> void AddKinKalTrajectory( std::unique_ptr<KTRAJ> &trajectory, REX::REveElement* &scene, unsigned int j, std::string kaltitle, double& t1, double& t2);
          void FillKinKalTrajectory(REX::REveManager *&eveMng, bool firstloop, REX::REveElement* &scene, std::tuple<std::vector<std::string>, std::vector<const KalSeedPtrCollection*>> track_tuple, bool kalinter, bool hits, double& t1, double& t2);
          void AddCosmicTrackFit(REX::REveManager *&eveMng, bool firstLoop_, const mu2e::CosmicTrackSeedCollection *cosmiccol, REX::REveElement* &scene);

          #endif
          ClassDef(REveMu2eDataInterface, 0);
      };
}


#endif
