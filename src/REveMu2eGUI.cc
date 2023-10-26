#include "REve/inc/REveMu2eGUI.hh"

using namespace mu2e;


int REveMu2eGUI::WriteCoreJson(nlohmann::json &j, int rnr_offset)
{
  j["path"] = "Event/SubRun/Run"; 
  //j["count"] = fCount;
  //j["total"] = fTotal;
  j["eventid"] = feventid;
  j["subrunid"] = fsubrunid;
  j["runid"] = frunid;
  j["UT_PostStream"] = "UT_refresh_event_info";
  return ROOT::Experimental::REveElement::WriteCoreJson(j, 0);
}

void REveMu2eGUI::PrintEveInfo(){
   std::cout<<"Run : "<<frunid<<" SubRun : "<<fsubrunid<<" Event : "<<feventid<<std::endl;
   PrintSimInfo();
   PrintKalInfo();
   PrintCaloInfo();
  }

void REveMu2eGUI::PrintSimInfo(){
 // SimParticle information
  std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(fmctrack_tuple);
  for(unsigned int j=0; j< track_list.size(); j++){
    const MCTrajectoryCollection* trajcol = track_list[j];
    if(trajcol !=0){
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      std::cout<<"SIM PARTICLE INFORMATION"<<std::endl;
      std::cout<<"Number of SimParticles = "<<trajcol->size()<<std::endl;
      std::cout<<"  "<<std::endl;
      std::cout<<" ID  "<<"   PDGID   "<<" Energy    "<<"    p0x     "<<"     p0y      "<<"       p0z    "
      <<"      posx    "<<"     posy      "<<"    posz    "<<"  t0      "<<"     p1x   "<<"    p1y   "<<"   p1z    "<<"    t1"<<std::endl;
      for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++){
         std::string pdgID = std::to_string(trajectoryIter->first->pdgId());
	 std::string t0 = std::to_string(trajectoryIter->first->startGlobalTime());
	 std::string p0x = std::to_string(trajectoryIter->first->startMomentum().x());
         std::string p0y = std::to_string(trajectoryIter->first->startMomentum().y());
         std::string p0z = std::to_string(trajectoryIter->first->startMomentum().z());
	 std::string energy = std::to_string(trajectoryIter->first->startMomentum().e());
         std::string posx = std::to_string(trajectoryIter->first->startPosition().x());
         std::string posy = std::to_string(trajectoryIter->first->startPosition().y());
         std::string posz = std::to_string(trajectoryIter->first->startPosition().z());
         std::string t1 = std::to_string(trajectoryIter->first->endGlobalTime());
         std::string p1x = std::to_string(trajectoryIter->first->endMomentum().x());
         std::string p1y = std::to_string(trajectoryIter->first->endMomentum().y());
         std::string p1z = std::to_string(trajectoryIter->first->endMomentum().z());
         std::string pos1x = std::to_string(trajectoryIter->first->endPosition().x());
         std::string pos1y = std::to_string(trajectoryIter->first->endPosition().y());
         std::string pos1z = std::to_string(trajectoryIter->first->endPosition().z());
	       std::string id = std::to_string(trajectoryIter->first->id().asInt());

        std::cout<<" "<<id<<"      "<<pdgID<<"     "<<energy<<"    "<<p0x<<"     "<<p0y<<"     "<<p0z
        <<"  "<<posx<<"  "<<posy<<"   "<<posz<<"  "<<t0<<"   "<<p1x<<"  "<<p1y<<"  "<<p1z<<"   "<<t1<<std::endl;	
      }    	
    }
  }
}

void  REveMu2eGUI::PrintKalInfo(){

 // KalSeed info
  std::vector<const KalSeedCollection*> ktrack_list = std::get<1>(ftrack_tuple);
  for(unsigned int j=0; j< ktrack_list.size(); j++){
    const KalSeedCollection* seedcol = ktrack_list[j];
    std::cout<<" "<<std::endl;
    std::cout<<"KALSEED INFORMATION"<<std::endl;  
    if(seedcol->size() !=0){
       for(unsigned int i = 0; i < seedcol->size(); i++){
         mu2e::KalSeed const  &kseed= (*seedcol)[i];
         const std::vector<mu2e::TrkStrawHitSeed>* hots = &kseed.hits();
         int n_krep_hits = hots->size();
	 std::string kt0 = std::to_string(kseed.t0().t0());
	 const std::vector<mu2e::KalSegment> &segments = kseed.segments();
         unsigned int nSegments=segments.size();
         std::cout<<" t0 = "<<kt0<<" hits = "<<n_krep_hits<<" segments = "<<nSegments<<std::endl; 
         // Print more details about the track hits
         /*for(int ih=0; ih<n_krep_hits; ++ih) {
           const mu2e::TrkStrawHitSeed* hit =  &hots->at(ih);
           if((hit != nullptr) and (hit->flag().hasAnyProperty(mu2e::StrawHitFlagDetail::active))) {
	     mu2e::StrawId sid = hit->strawId(); 
	     // const mu2e::Straw* straw = &tracker->straw(sid);
             hit->driftRadius();
	     // const CLHEP::Hep3Vector* v1 = &straw->getMidPoint();
	     std::cout<<sid<<"    "<<hit->driftRadius()<<"    "<<hit->_edep<<"    "<<hit->_tottdrift<<std::endl; 
	   }
	 }*/                
       }
    }
  }
}

void REveMu2eGUI::PrintCaloInfo(){
 // CaloClusterCollection info
  std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(fcalocluster_tuple);
  if(calocluster_list.size()!=0){
    for(unsigned int j = 0; j< calocluster_list.size(); j++){
      const CaloClusterCollection* clustercol = calocluster_list[j];
      std::cout<<" "<<std::endl;
      std::cout<<"CALO CLUSTER INFORMATION"<<std::endl;
      std::cout<<"  Energy  "<<"     Time  "<<"        X       "<<"      Y    "<<std::endl;
      if(clustercol->size() != 0){
       for(unsigned int i = 0; i < clustercol->size(); i++){
         mu2e::CaloCluster const  &cluster= (*clustercol)[i];
         // Info to print:
         std::string cluster_energy = std::to_string(cluster.energyDep());
         std::string cluster_time = std::to_string(cluster.time());
         std::string cluster_x = std::to_string(cluster.cog3Vector().x());
         std::string cluster_y = std::to_string(cluster.cog3Vector().y());
         std::cout<<" "<<cluster_energy<<"   "<<cluster_time<<"   "<<cluster_x<<"    "<<cluster_y<<std::endl;
       }
      }
    }
  }
}		
