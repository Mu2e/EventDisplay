#ifndef _REveMu2ePrintInfo_hh
#define _REveMu2ePrintInfo_hh

#include <ROOT/REveElement.hxx>
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include <iostream>
#include <vector>
#include <tuple>
#include <string>

namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class REveMu2ePrintInfo : public ROOT::Experimental::REveElement
{
  public:
     REveMu2ePrintInfo() = default;
     const mu2e::CaloClusterCollection* clustercol = 0;
     std::vector<const mu2e::CaloClusterCollection*> calocluster_list;
     std::vector<std::string> calocluster_labels;
     std::tuple<std::vector<std::string>, std::vector<const mu2e::CaloClusterCollection*>> fcalocluster_tuple;

     const mu2e::MCTrajectoryCollection *mctrajcol = 0;
     std::vector<const mu2e::MCTrajectoryCollection*> mctrack_list;
     std::vector<std::string> mctrack_labels;
     std::tuple<std::vector<std::string>, std::vector<const mu2e::MCTrajectoryCollection*>> fmctrack_tuple;

     const mu2e::KalSeedCollection* kalSeedcol = 0;
     std::vector<const mu2e::KalSeedCollection*> track_list;
     std::vector<std::string> track_labels;
     std::tuple<std::vector<std::string>, std::vector<const mu2e::KalSeedCollection*>> ftrack_tuple;
     
     const mu2e::CrvCoincidenceClusterCollection* crvcoincol = 0;
     std::vector<const mu2e::CrvCoincidenceClusterCollection*> crvcoin_list;
     std::vector<std::string> crvcoin_labels;
     std::tuple<std::vector<std::string>, std::vector<const mu2e::CrvCoincidenceClusterCollection*>> fcrvcoin_tuple;

     void PrintRecoInfo(); //prints everything which is there (tracks, clusters, CRV coin, not hits)
     void PrintMCInfo(); // print MCs and SimParticles
     void PrintSimInfo();
     void PrintKalInfo();
     void PrintCaloInfo();
     void PrintCRVInfo();
};
}

#endif
