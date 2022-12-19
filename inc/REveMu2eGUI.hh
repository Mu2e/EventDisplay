#ifndef _REveMu2eGUI_hh
#define _REveMu2eGUI_hh

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

#include <iostream>
#include <vector>
#include <tuple>
#include <string>

#include <ROOT/REveElement.hxx>
#include "nlohmann/json.hpp"
namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class REveMu2eGUI : public ROOT::Experimental::REveElement
{
  public:
 
     //int fCount{0};
     //int fTotal{10};
     int feventid{0};
     int fsubrunid{0};
     int frunid{0};
     
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

     int WriteCoreJson(nlohmann::json &j, int rnr_offset) override;
     void PrintEveInfo();
     void PrintSimInfo();
     void PrintKalInfo();
     void PrintCaloInfo();

};
}

#endif
