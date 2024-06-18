#include "REve/inc/REveMu2ePrintInfo.hh"

using namespace mu2e;

void REveMu2ePrintInfo::PrintMCInfo(){
  PrintSimInfo();
}

void REveMu2ePrintInfo::PrintRecoInfo(){
  PrintKalInfo();
  PrintCaloInfo();
  PrintCRVInfo();
}


void REveMu2ePrintInfo::PrintSimInfo(){

  std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(fmctrack_tuple);
  if(track_list.size() > 0){
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
}

void  REveMu2ePrintInfo::PrintKalInfo(){
    std::vector<const KalSeedCollection*> ktrack_list = std::get<1>(ftrack_tuple);
    std::vector<std::string> names = std::get<0>(ftrack_tuple);
    if(ktrack_list.size() > 0){
      for(unsigned int j=0; j< ktrack_list.size(); j++){
        const KalSeedCollection* seedcol = ktrack_list[j];
        std::cout<<" "<<std::endl;
        if(seedcol->size() !=0){
          std::cout<<"KALSEED INFORMATION"<<std::endl;
          for(unsigned int k = 0; k < seedcol->size(); k++){
            mu2e::KalSeed kseed = (*seedcol)[k];
            const std::vector<mu2e::KalSegment> &segments = kseed.segments();
            if(kseed.loopHelixFit())
            {
              auto trajectory=kseed.loopHelixFitTrajectory();
              auto lh = segments[0].loopHelix();
              std::string kaltitle = "Loop Helix KalSeed tag : " + names[j] +  '\n'
                  + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
                  + " t0 " + std::to_string(lh.t0()) + "ns " +  '\n'
                  + " lam "  + std::to_string(lh.lam() ) +  '\n'
                  + " rad "  + std::to_string(lh.rad() ) +  '\n'
                  + " cx "  + std::to_string(lh.cx() ) +  '\n'
                  + " cy "  + std::to_string(lh.cy() ) +  '\n'
                  + " phi0 "  + std::to_string(lh.phi0() );
                  std::cout<<kaltitle<<std::endl;
            }
            else if(kseed.centralHelixFit())
            {
              auto ch = segments[0].centralHelix();
              std::string kaltitle = "Central Helix KalSeed tag : " + names[j] +  '\n'
                + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
                + " t0 " + std::to_string(ch.t0()) + "ns " +  '\n'
                + " tandip " + std::to_string(ch.tanDip() ) +  '\n'
                + " d0 " + std::to_string(ch.d0() ) +  '\n'
                + " z0 " + std::to_string(ch.z0() ) +  '\n'
                + " phi0 " + std::to_string(ch.phi0() ) +  '\n'
                + " omega " + std::to_string(ch.omega() );;
                std::cout<<kaltitle<<std::endl;
            }
            else if(kseed.kinematicLineFit())
            {
               auto kl = segments[0].kinematicLine();
               std::string kaltitle = "KinematicLine KalSeed tag : " + names[j]
                  + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
                  + " t0 " + std::to_string(kl.t0()) + "ns " +  '\n'
                  + " d0 " + std::to_string(kl.d0() ) +  '\n'
                  + " z0 " + std::to_string(kl.z0() ) +  '\n'
                  + " phi0 " + std::to_string(kl.phi0() ) +  '\n'
                  + " theta " + std::to_string(kl.theta() ) ;
                  std::cout<<kaltitle<<std::endl;
            }
          }
      }
    }
  }
}

void REveMu2ePrintInfo::PrintCaloInfo(){
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

void REveMu2ePrintInfo::PrintCRVInfo(){
  std::cout<<"CRV COINCIDENCE INFORMATION"<<std::endl;
  std::vector<const CrvCoincidenceClusterCollection*> crvcoin_list = std::get<1>(fcrvcoin_tuple);
  if(crvcoin_list.size()!=0){
  for(unsigned int j = 0; j< crvcoin_list.size(); j++){
      const CrvCoincidenceClusterCollection* crvcoincol = crvcoin_list[j];
      std::cout<<"  avTime  "<<"        avX       "<<"      avY    "<<"      avZ       "<<std::endl;
      if(crvcoincol->size() != 0){
       for(unsigned int i = 0; i < crvcoincol->size(); i++){
         mu2e::CrvCoincidenceCluster const  &cluster= (*crvcoincol)[i];
         std::string coin_avtime = std::to_string(cluster.GetAvgHitTime());
         std::string coin_avX = std::to_string(cluster.GetAvgHitPos().x());
         std::string coin_avY = std::to_string(cluster.GetAvgHitPos().y());
         std::string coin_avZ = std::to_string(cluster.GetAvgHitPos().z());
         std::cout<<"    "<<coin_avtime<<"              "<<coin_avX<<"          "<<coin_avY<<"             "<<coin_avZ<<std::endl;
         }
       }
     }
  }
}
