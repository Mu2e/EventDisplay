#include "REve/inc/REveMu2eDataInterface.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
using namespace mu2e;
namespace REX = ROOT::Experimental;

std::string drawfilename("REve/config/drawutils.txt");
SimpleConfig drawconfig(drawfilename);


void REveMu2eDataInterface::AddCaloDigis(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloDigiCollection*>> calodigi_tuple, REX::REveElement* &scene, bool addCrystalDraw){

  std::cout<<"[REveMu2eDataInterface] AddCaloDigi "<<std::endl;
  std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(calodigi_tuple);
  std::vector<std::string> names = std::get<0>(calodigi_tuple);
  for(unsigned int j = 0; j< calodigi_list.size(); j++){

    const CaloDigiCollection* calodigicol = calodigi_list[j];
    if(calodigicol->size() != 0){
      if(!firstLoop_){
        scene->DestroyElements();;
      }
      mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
      for(unsigned int i = 0; i < calodigicol->size(); i++){
        mu2e::CaloDigi const  &digi= (*calodigicol)[i];
        int sipmID = digi.SiPMID();
        //bool isEven = sipmID % 2;
        int cryID = sipmID/2;
        //std::cout<<"crystal ID "<<cryID<<" "<<isEven<<" "<<digi.t0()<<std::endl;
        Crystal const &crystal = cal.crystal(cryID);
        double crystalXLen = pointmmTocm(crystal.size().x());
        double crystalYLen = pointmmTocm(crystal.size().y());
        double crystalZLen = pointmmTocm(crystal.size().z());


        GeomHandle<DetectorSystem> det;
        double zpos = 0;
        double diskID = 0;
        if(cryID < 674) {
          zpos = 186.53;
          diskID = 0;
        }
        if(cryID >= 674) {
          zpos = 256.53;
          diskID = 1;
        }
        CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(diskID,crystal.position()) ;

        std::string crytitle =  "Crystal ID = " + std::to_string(cryID) +  '\n'
          + " Time = " + std::to_string(digi.t0())+" ns ";
        char const *crytitle_c = crytitle.c_str();

        // plot the crystals which are present in this event:
        auto b = new REX::REveBox("crystal",crytitle_c);
        b->SetMainColor(602);
        double width = crystalXLen/2;
        double height = crystalYLen/2;
        double thickness = crystalZLen/2;
        b->SetVertex(0, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())- thickness + 2*thickness + zpos);//---
        b->SetVertex(1, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height, pointmmTocm(crystalPos.z())- thickness + 2*thickness  +zpos);//-+-
        b->SetVertex(2, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height ,pointmmTocm(crystalPos.z())- thickness + 2*thickness  + zpos);//++-
        b->SetVertex(3, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())-thickness  + 2*thickness + zpos);//+--
        b->SetVertex(4, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())+ thickness + 2*thickness + zpos);//--+
        b->SetVertex(5, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z())+ thickness  + 2*thickness + zpos);//-++
        b->SetVertex(6, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z()) + thickness + 2*thickness+zpos); //+++
        b->SetVertex(7,pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())+ thickness + 2*thickness + zpos);//+-+
        scene->AddElement(b);
      }
    }
  }
}
void REveMu2eDataInterface::AddCaloClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloClusterCollection*>> calocluster_tuple, REX::REveElement* &scene, bool addCrystalDraw, double t1, double t2){

  std::cout<<"[REveMu2eDataInterface] AddCaloClusters "<<std::endl;
  std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(calocluster_tuple);
  std::vector<std::string> names = std::get<0>(calocluster_tuple);
  for(unsigned int j = 0; j< calocluster_list.size(); j++){

    // Extract cluster list
    const CaloClusterCollection* clustercol = calocluster_list[j];
    if(clustercol->size() != 0){
      if(!firstLoop_){
        scene->DestroyElements();;
      }
      Color_t color = kWhite;
      mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
      GeomHandle<DetectorSystem> det;
      // Add crystals
      double maxE = 1e-6;
      double minE = 1000;
      std::list<CaloCluster> cluList;
      //std::array<Color_t, 10> colorList {kViolet, kBlue, kBlue, kGreen, kGreen, kYellow, kYellow, kOrange, kOrange, kRed};
      //std::vector<Color_t> colors;

      for(unsigned int i = 0; i < clustercol->size(); i++){
        mu2e::CaloCluster const  &cluster= (*clustercol)[i];
        cluList.push_back(cluster);
        if(cluster.energyDep() > maxE) maxE = cluster.energyDep();
        if(cluster.energyDep() < minE) minE = cluster.energyDep();

      // sort in energy
      cluList.sort([] (CaloCluster lhs, CaloCluster rhs) {return lhs.energyDep() > rhs.energyDep();} );
      }

      /*if(cluList.size() <= 2){
          Color_t color;
          for(auto const& cluster : cluList){
            if(cluster.energyDep() == minE){
              color = kViolet;
              colors.push_back(color);
            }
            if(cluster.energyDep() == maxE){
              color = kRed;
              colors.push_back(color);
            }
        }
      }
      unsigned int c = -1;
      if(cluList.size() > 2){
        for(auto const& cluster : cluList){
          c+=1;
          Color_t color;
          if(cluster.energyDep() == minE){
            color = kViolet;
            colors.push_back(color);
          }
          if(cluster.energyDep() == maxE){
            color = kRed;
            colors.push_back(color);
          }
          if (c > 1 and c < colorList.size()-1){
            color = colorList.at(c);
            colors.push_back(color);
          }
          if (c > 1 and c >= colorList.size()-1){
            color = kRed;
            colors.push_back(color);
          }
        }
      } */

      //unsigned int i = -1;
      for(auto const& cluster : cluList){
        //mu2e::CaloCluster const  &cluster= cluList.at(i);
        //i += 1;

        // Info for label:
        std::string cluster_energy = std::to_string(cluster.energyDep());
        std::string cluster_time = std::to_string(cluster.time());
        std::string cluster_x = std::to_string(cluster.cog3Vector().x());
        std::string cluster_y = std::to_string(cluster.cog3Vector().y());
        // Extract center of gravity, convert to coord sys
        CLHEP::Hep3Vector COG(cluster.cog3Vector().x(),cluster.cog3Vector().y(), cluster.cog3Vector().z());
        CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDisk(cluster.diskID(),COG);
        CLHEP::Hep3Vector pointInMu2e = det->toMu2e(crystalPos);
        // Info for label
        std::string cluster_z = std::to_string(abs(pointInMu2e.z()));
        // Make label and REve objects
        std::string label =   " Instance = " + names[0] +  '\n'
          + " Energy Dep. = "+cluster_energy+" MeV "+   '\n'
          + " Time = "+cluster_time+" ns " +  '\n'
          + " Pos =  ("+cluster_x+","+cluster_y+","+cluster_z+") mm";
        std::string name = "disk" + std::to_string(cluster.diskID());
        auto ps1 = new REX::REvePointSet(name, "CaloClusters Disk 1: "+label,0);
        auto ps2 = new REX::REvePointSet(name, "CaloClusters Disk 2: "+label,0);


        // Set positions of clusters
        if(cluster.diskID() == 0) ps1->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));
        if(cluster.diskID() == 1) ps2->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));

        // Set draw options
        /*TColor color1;
        color1.SetPalette(1,0);
        Color_t color = color1.GetNumber();// colors[i];*/

        if(abs(cluster.time() - t2) < 20) color = kRed;
        if(abs(cluster.time() - t2) >=  20 and abs(cluster.time() - t2) < 100) color = kOrange;
        if(abs(cluster.time() - t2) >= 100 and abs(cluster.time() - t2) < 200) color = kYellow;
        if(abs(cluster.time() - t2) >= 200 and abs(cluster.time() - t2) < 400) color = kGreen;
        if(abs(cluster.time() - t2) >= 400 and abs(cluster.time() - t2) < 600) color = kBlue;
        if(abs(cluster.time() - t2) >= 600 ) color = kViolet;

        ps1->SetMarkerColor(color);
        ps1->SetMarkerStyle(REveMu2eDataInterface::mstyle);
        ps1->SetMarkerSize(REveMu2eDataInterface::mstyle);

        ps2->SetMarkerColor(color);
        ps2->SetMarkerStyle(REveMu2eDataInterface::mstyle);
        ps2->SetMarkerSize(REveMu2eDataInterface::mstyle);

        // Add to REve world
        scene->AddElement(ps1);
        scene->AddElement(ps2);

        // Add crystals
        if(addCrystalDraw){
          for(unsigned h =0 ; h < cluster.caloHitsPtrVector().size();h++)     {

            art::Ptr<CaloHit>  crystalhit = cluster.caloHitsPtrVector()[h];
            int cryID = crystalhit->crystalID();

            //int diskID = cal.crystal(crystalhit->crystalID()).diskID();
            Crystal const &crystal = cal.crystal(cryID);
            double crystalXLen = pointmmTocm(crystal.size().x());
            double crystalYLen = pointmmTocm(crystal.size().y());
            double crystalZLen = pointmmTocm(crystal.size().z());

            GeomHandle<DetectorSystem> det;
            CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(cluster.diskID(),crystal.position()) ;

            //  make title
            std::string crytitle =   "disk"+std::to_string(cal.crystal(crystalhit->crystalID()).diskID()) + " Crystal ID = " + std::to_string(cryID) +  '\n'
              + " Energy Dep. = "+std::to_string(crystalhit->energyDep())+" MeV "+   '\n'
              + " Time = "+std::to_string(crystalhit->time())+" ns ";
            char const *crytitle_c = crytitle.c_str();
            auto b = new REX::REveBox(crytitle_c,crytitle_c);

            // plot the crystals which are present in this event in lego:
            b->SetMainColor(color);//s[i]);


            double width = crystalXLen/2;
            double height = crystalYLen/2;
            double thickness = crystalhit->energyDep()/maxE * crystalZLen/2; //length proportional to energy

            b->SetVertex(0, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())- thickness  + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//---
            b->SetVertex(1, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height, pointmmTocm(crystalPos.z())- thickness   +abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//-+-
            b->SetVertex(2, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height ,pointmmTocm(crystalPos.z())- thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//++-
            b->SetVertex(3, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())-thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//+--
            b->SetVertex(4, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())+ thickness  + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//--+
            b->SetVertex(5, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z())+ thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//-++
            b->SetVertex(6, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z()) + thickness +abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2); //+++
            b->SetVertex(7,pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())+ thickness + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//+-+
            scene->AddElement(b);
          }
        }

      }
    }
  }
}

void REveMu2eDataInterface::AddComboHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple, REX::REveElement* &scene, bool strawdisplay, bool AddErrorBar_){

  std::vector<const ComboHitCollection*> combohit_list = std::get<1>(combohit_tuple);
  std::vector<std::string> names = std::get<0>(combohit_tuple);

  // Loop over hit lists
  for(unsigned int j = 0; j < combohit_list.size(); j++){
    const ComboHitCollection* chcol = combohit_list[j];
    int colour = (j+3);
    if(chcol->size() !=0 ){
      // Loop over hits
      for(unsigned int i=0; i< chcol->size(); i++){
        mu2e::ComboHit const  &hit= (*chcol)[i];
        // Display hit straws if selected too in FCL file
        if(strawdisplay){
          mu2e::GeomHandle<mu2e::Tracker> tracker;
          const auto& allStraws = tracker->getStraws();
          int sid = hit._sid.asUint16();
          CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
          const mu2e::Straw& s = allStraws[sid];
          const CLHEP::Hep3Vector& p = s.getMidPoint();
          const CLHEP::Hep3Vector& d = s.getDirection();
          double x = p.x();
          double y = p.y();
          double z = p.z();
          double l = s.halfLength();
          double st=sin(d.theta());
          double ct=cos(d.theta());
          double sp=sin(d.phi()+TMath::Pi()/2.0);
          double cp=cos(d.phi()+TMath::Pi()/2.0);
          if(sid < drawconfig.getInt("maxSID")){
            double x1=x+l*st*sp;
            double y1=y-l*st*cp;
            double z1=z+l*ct;
            double x2=x-l*st*sp;
            double y2=y+l*st*cp;
            double z2=z-l*ct;
            std::string strawtitle;
            int idStraw =  s.id().getStraw();
            int idPanel =  s.id().getPanel();
            int idPlane =  s.id().getPlane();
            int colorid = idPlane + idPanel;
            strawtitle =Form("Straw %i Panel %i  Plane %i",idStraw,idPanel,idPlane);
            sposi.set(x1, y1, z1);
            sposf.set(x2, y2, z2);
            if(sposi.x()!=0){
              auto strawline = new REX::REveLine("StrawHit ",strawtitle,2);
              strawline->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y()),pointmmTocm(sposi.z()));
              strawline->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y()),pointmmTocm(sposf.z()));
              strawline->SetLineWidth(1);
              strawline->SetLineColor(colorid);
              if(strawline->GetSize() !=0 ) scene->AddElement(strawline);
            }
          }
        }

        if(AddErrorBar_){
          //XY
          auto const& p = hit.pos();
          auto w = hit.uDir();
          auto const& s = hit.wireRes();
          double x1 = (p.x()+s*w.x());
          double x2 = (p.x()-s*w.x());
          double z1 = (p.z()+s*w.z());
          double z2 = (p.z()-s*w.z());
          double y1 = (p.y()+s*w.y());
          double y2 = (p.y()-s*w.y());

          std::string errorbar = Form("ErrorBar Length: %f, %f, %f, %f, %f, %f",x1,y1,z1,x2,y2,z2);
          auto error = new REX::REveLine("errors",errorbar.c_str(),2);
          error->SetPoint(0, pointmmTocm(x1),pointmmTocm(y1),pointmmTocm(z1));
          error->SetNextPoint(pointmmTocm(x2),pointmmTocm(y2),pointmmTocm(z2));
          //std::cout<<"points "<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
          //std::cout<<"errors "<<x1<<" "<<y1<<" "<<z1<<std::endl;
          error->SetLineColor(kRed);
          error->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
          scene->AddElement(error);
        }
        CLHEP::Hep3Vector HitPos(pointmmTocm(hit.pos().x()), pointmmTocm(hit.pos().y()), pointmmTocm(hit.pos().z()));
        std::string chtitle = "ComboHits tag = "
          + (names[j])  +  '\n'
          + " position : x "  + std::to_string(hit.pos().x())  +  '\n'
          + " y " + std::to_string(hit.pos().y())  +  '\n'
          + " z " + std::to_string(hit.pos().z())  +  '\n'
          + " time :" + std::to_string(hit.time()) +  '\n'
          + " energy dep : "
          + std::to_string(hit.energyDep())  +
          + "MeV";
        auto ps1 = new REX::REvePointSet("ComboHit", chtitle,0);
        ps1->SetNextPoint(HitPos.x(), HitPos.y() , HitPos.z());
        ps1->SetMarkerColor(colour);
        ps1->SetMarkerStyle(REveMu2eDataInterface::mstyle);
        ps1->SetMarkerSize(REveMu2eDataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }
}

/*------------Function to add CRV information to the display:-------------*/
void REveMu2eDataInterface::AddCRVInfo(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvRecoPulseCollection*>>  crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars){
  std::cout<<"[ REveDataInterface::AddCRVInfo() ]"<<std::endl;
  std::vector<const CrvRecoPulseCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
  std::vector<std::string> names = std::get<0>(crvpulse_tuple);
  GeomHandle<CosmicRayShield> CRS;
  GeomHandle<DetectorSystem> det;
  if(crvpulse_list.size() !=0){
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
      const CrvRecoPulseCollection* crvRecoPulse = crvpulse_list[i];
      if(crvRecoPulse->size() !=0){
        std::string crvtitle = " tag : " + names[i];
        auto ps1 = new REX::REvePointSet("CRVRecoPulse", crvtitle,0);
        for(unsigned int j=0; j< crvRecoPulse->size(); j++){
          mu2e::CrvRecoPulse const &crvpulse = (*crvRecoPulse)[j];
          const CRSScintillatorBarIndex &crvBarIndex = crvpulse.GetScintillatorBarIndex();
          const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
          const CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
          CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
          CLHEP::Hep3Vector HitPos(crvCounterPos.x(), crvCounterPos.y(), crvCounterPos.z());

          CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvCounterPos);
          CLHEP::Hep3Vector sibardetails(barDetail.getHalfLengths()[0],barDetail.getHalfLengths()[1],barDetail.getHalfLengths()[2]);

          if(addCRVBars){
            if(!extracted){
              auto b = new REX::REveBox("box","label");
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b->	SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){ //CRV D, CRV U -- WORKS
                b->SetVertex(0,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height); //---
                b->SetVertex(1,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height);//-+-
                b->SetVertex(2,  -1*length, pointmmTocm(pointInMu2e.y()) + width , pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//++-
                b->SetVertex(3,  -1*length, pointmmTocm(pointInMu2e.y()) + width  ,pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - width,   pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height  );//--+
                b->SetVertex(5,  length, pointmmTocm(pointInMu2e.y()) - width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//+++
                b->SetVertex(7, length,pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height );//+-+

              }

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //CRV TS -- WORKS
                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length ); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//+-+
              }

              if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){ //CRV T -- orientation DONE
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width); //---
                b->SetVertex(1, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//--+
                b->SetVertex(5, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+-+
              }

               if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 2){ //CRV R, CRV L -- orientation DONE
                b->SetVertex(0, pointmmTocm(pointInMu2e.y()) - height, -1*length, pointmmTocm(pointInMu2e.x()) - width );//---
                b->SetVertex(1, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.y()) - height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.y()) - height, length, pointmmTocm(pointInMu2e.x()) - width);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.y()) + height, length, pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.y()) + height, length,pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.y()) - height, length,pointmmTocm(pointInMu2e.x()) + width );//+-+
              }
              scene->AddElement(b);

            }
            if(extracted){ //TODO same for nominal geom

              // CRV hit scintillation bars highlighted
              // std::string const& base;
              // std::string bartitle = crvCounter.name( base );
              // char const *bartitle_c = base.c_str(); //TODO title

              // Draw "bars hit" in red:
              auto b = new REX::REveBox("box");
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b->	SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //T1

                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, -1*length); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y())  + height, -1*length);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, -1*length);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, -1*length);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height,  length );//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, length );//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height ,length );//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, length );//+-+
              } else { //EX, T2
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//---
                b->SetVertex(1,-1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//-+-
                b->SetVertex(2, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width  + pointmmTocm(pointInMu2e.z()));//++-
                b->SetVertex(3, length, pointmmTocm(pointInMu2e.y()) - height , pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()) );//+--
                b->SetVertex(4, -1*length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z())); //--+
                b->SetVertex(5,-1*length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()) ); // -++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+++
                b->SetVertex(7,  length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+-+

              }
              scene->AddElement(b);
            }
          }
          // Add Reco Pulse to collection
          ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()) , pointmmTocm(pointInMu2e.z()));

        }
        // Draw reco pulse collection
        ps1->SetMarkerColor(i+3);
        ps1->SetMarkerStyle(REveMu2eDataInterface::mstyle);
        ps1->SetMarkerSize(REveMu2eDataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }

}


/*------------Function to add CRV information to the display:-------------*/
void REveMu2eDataInterface::AddCRVClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvCoincidenceClusterCollection*>>  crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars){

  std::vector<const CrvCoincidenceClusterCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
  std::vector<std::string> names = std::get<0>(crvpulse_tuple);
  GeomHandle<CosmicRayShield> CRS;
  GeomHandle<DetectorSystem> det;
  if(crvpulse_list.size() !=0){
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
      const CrvCoincidenceClusterCollection* crvClusters = crvpulse_list[i];
      if(crvClusters->size() !=0){
        std::string crvtitle = " tag : " + names[i];
        auto ps1 = new REX::REvePointSet("CRVCoincidenceClusters", crvtitle,0);
        for(unsigned int j=0; j< crvClusters->size(); j++){
          mu2e::CrvCoincidenceCluster const &crvclu = (*crvClusters)[j];
          CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvclu.GetAvgHitPos());
          ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()) , pointmmTocm(pointInMu2e.z()));
          for(unsigned h =0 ; h < crvclu.GetCrvRecoPulses().size();h++)     {

            art::Ptr<CrvRecoPulse>  crvpulse = crvclu.GetCrvRecoPulses()[h];
            const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvpulse->GetScintillatorBarIndex();
            const mu2e::CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
            const mu2e::CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
            CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
            CLHEP::Hep3Vector HitPos(crvCounterPos.x(), crvCounterPos.y(), crvCounterPos.z());

            CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvCounterPos);
            CLHEP::Hep3Vector sibardetails(barDetail.getHalfLengths()[0],barDetail.getHalfLengths()[1],barDetail.getHalfLengths()[2]);

            if(addCRVBars){
            if(!extracted){
              auto b = new REX::REveBox("box","label");
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b->	SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){ //CRV D, CRV U -- WORKS
                b->SetVertex(0,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height); //---
                b->SetVertex(1,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height);//-+-
                b->SetVertex(2,  -1*length, pointmmTocm(pointInMu2e.y()) + width , pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//++-
                b->SetVertex(3,  -1*length, pointmmTocm(pointInMu2e.y()) + width  ,pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - width,   pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height  );//--+
                b->SetVertex(5,  length, pointmmTocm(pointInMu2e.y()) - width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//+++
                b->SetVertex(7, length,pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height );//+-+

              }

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //CRV TS -- WORKS
                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length ); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//+-+

              }

              if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){ //CRV T -- WORKS
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width); //---
                b->SetVertex(1, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//--+
                b->SetVertex(5, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+-+

              }

               if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 2){ //CRV R, CRV L -- orientation DONE
                b->SetVertex(0, pointmmTocm(pointInMu2e.y()) - height, -1*length, pointmmTocm(pointInMu2e.x()) - width ); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.y()) - height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.y()) - height, length, pointmmTocm(pointInMu2e.x()) - width);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.y()) + height, length, pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.y()) + height, length,pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.y()) - height, length,pointmmTocm(pointInMu2e.x()) + width );//+-+

              }
              scene->AddElement(b);

            }
            if(extracted){ //TODO same for nominal geom

              // CRV hit scintillation bars highlighted
              // std::string const& base;
              // std::string bartitle = crvCounter.name( base );
              // char const *bartitle_c = base.c_str(); //TODO title

              // Draw "bars hit" in red:
              auto b = new REX::REveBox("box");
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b->	SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //T1

                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, -1*length); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y())  + height, -1*length);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, -1*length);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, -1*length);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height,  length );//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, length );//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height ,length );//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, length );//+-+
              } else { //EX, T2
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//---
                b->SetVertex(1,-1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//-+-
                b->SetVertex(2, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width  + pointmmTocm(pointInMu2e.z()));//++-
                b->SetVertex(3, length, pointmmTocm(pointInMu2e.y()) - height , pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()) );//+--
                b->SetVertex(4, -1*length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z())); //--+
                b->SetVertex(5,-1*length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()) ); // -++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+++
                b->SetVertex(7,  length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+-+

              }
              scene->AddElement(b);
            }
          }
         }
        }

        ps1->SetMarkerColor(drawconfig.getInt("CRVHitColor"));
        ps1->SetMarkerStyle(REveMu2eDataInterface::mstyle);
        ps1->SetMarkerSize(REveMu2eDataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }

}

/*------------Function to add TimeCluster Collection in 3D and 2D displays:-------------*/
void REveMu2eDataInterface::AddTimeClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const TimeClusterCollection*>>  timecluster_tuple, REX::REveElement* &scene){

  std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(timecluster_tuple);
  std::vector<std::string> names = std::get<0>(timecluster_tuple);

  if(timecluster_list.size() !=0){
    for(unsigned int i = 0; i < timecluster_list.size(); i++){
      const TimeClusterCollection* tccol = timecluster_list[i];
      if(tccol->size() != 0){
        if(!firstLoop_){
          scene->DestroyElements();;
        }


        for(size_t j=0; j<tccol->size();j++){
          mu2e::TimeCluster const  &tclust= (*tccol)[j];
          std::string tctitle = "Time Cluster tag: " + names[i] + '\n'
          + "t0 " + std::to_string(tclust.t0().t0()) + " +/- " + std::to_string(tclust.t0().t0Err()) + " ns " + '\n' ;
          auto ps1 = new REX::REvePointSet("TimeClusters", tctitle, 0);
          CLHEP::Hep3Vector HitPos(tclust._pos.x(), tclust._pos.y(), tclust._pos.z());
          ps1->SetNextPoint(pointmmTocm(HitPos.x()), pointmmTocm(HitPos.y()) , pointmmTocm(HitPos.z()));
          ps1->SetMarkerColor(i+6);
          ps1->SetMarkerStyle(kOpenCircle);
          ps1->SetMarkerSize(REveMu2eDataInterface::msize);
          if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
        }
      }
    }
  }
}


void REveMu2eDataInterface::AddHelixSeedCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const HelixSeedCollection*>> helix_tuple, REX::REveElement* &scene){
  std::cout<<"[REveMu2eDataInterface] AddHelices "<<std::endl;
  std::vector<const HelixSeedCollection*> helix_list = std::get<1>(helix_tuple);
  std::vector<std::string> names = std::get<0>(helix_tuple);
  for(unsigned int j=0; j< helix_list.size(); j++){
    const HelixSeedCollection* seedcol = helix_list[j];
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k++){
        mu2e::HelixSeed hseed = (*seedcol)[k];
        const ComboHitCollection& hhits = hseed.hits();
        unsigned int nhhits = hhits.size();
        auto line = new REX::REveLine(names[j], names[j],nhhits);
        // std::cout<<"Helix hits = "<<nhhits<<" radius = "<<hseed.helix()._radius<<" x0 = "<<hseed.helix()._rcent*cos(hseed.helix()._fcent)<<std::endl;
        float helrad = hseed.helix()._radius;
        float circphi = 0.0;
        float xc = hseed.helix()._rcent*cos(hseed.helix()._fcent);
        float yc = hseed.helix()._rcent*sin(hseed.helix()._fcent);
        if(hhits.size() !=0 ){
          for(int i=-1500; i < 1500; i +=100){ //TODO Remove hard-coded numbers
            float z = i;
            if(hseed.helix()._lambda !=0.0) circphi = hseed.helix()._fz0 + z/hseed.helix()._lambda;
            float x= xc + helrad*cos(circphi);
            float y = yc + helrad*sin(circphi);
            //std::cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<std::endl;
            CLHEP::Hep3Vector HelPos(x, y, z);
            if(circphi !=0.0)line->SetNextPoint(pointmmTocm(HelPos.x()),pointmmTocm(HelPos.y()) ,pointmmTocm(HelPos.z()));
          }
        }
        line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
        line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
        scene->AddElement(line);
      }
    }
  }
}

void REveMu2eDataInterface::AddKalIntersection(KalSeed kalseed, REX::REveElement* &scene){
  //Plot intersecitons:
  std::vector<mu2e::KalIntersection> inters = kalseed.intersections();

  for(unsigned int i = 0; i < inters.size(); i++){
  KinKal::VEC3 posKI = inters[i].position3();
  std::string title = "KalIntersection position : x "  + std::to_string(posKI.x())  +  '\n'
                  + " y " + std::to_string(posKI.y())  +  '\n'
                  + " z " + std::to_string(posKI.z())  +  '\n'
                  + " time :" + std::to_string(inters[i].time()) +  '\n'
                  + " mom : "
                  + std::to_string(inters[i].mom())  +
                  + "MeV/c " + '\n'
                  + "Surface ID: " + std::to_string(inters[i].surfaceId().index()) ;
    auto interpoint = new REX::REvePointSet("KalIntersections", title,1);

    interpoint->SetMarkerStyle(REveMu2eDataInterface::mstyle);
    interpoint->SetMarkerSize(REveMu2eDataInterface::msize);
    interpoint->SetMarkerColor(kViolet);
    interpoint->SetNextPoint(pointmmTocm(posKI.x()),pointmmTocm(posKI.y()) ,pointmmTocm(posKI.z()));
    if(interpoint->GetSize() !=0 ) scene->AddElement(interpoint);
  }

}

template<class KTRAJc> void REveMu2eDataInterface::AddTrkStrawHit(KalSeed kalseed, REX::REveElement* &scene,  std::unique_ptr<KTRAJc> &lhptr){
  std::cout<<"[REveMu2eDataInterface::AddTrkStrawHit]"<<std::endl;
  //Plot trk straw hits
  mu2e::GeomHandle<mu2e::Tracker> tracker;
  const std::vector<mu2e::TrkStrawHitSeed> &hits = kalseed.hits();

  for(unsigned int i = 0; i < hits.size(); i++){
    const mu2e::TrkStrawHitSeed &tshs = hits.at(i);
    auto const& straw = tracker->straw(tshs.strawId());
    mu2e::WireHitState whs(mu2e::WireHitState::State(tshs._ambig),
        mu2e::StrawHitUpdaters::algorithm(tshs._algo),
        tshs._kkshflag);
    bool active = whs.active();
    bool usedrift = whs.driftConstraint();
    if(active){ // maybe draw inactive hits but with a different color? TODO
      // then find the position at the reference POCA: start with the position on the wire
      auto tshspos = XYZVectorF(straw.wirePosition(tshs._rupos));
      double herr(0.0);
      //find the direction along DOCA
      auto tdir = lhptr->direction(tshs._ptoca);
      auto wdir = XYZVectorF(straw.wireDirection(tshs._rupos));
      auto ddir = wdir.Cross(tdir).Unit()*whs.lrSign();
      if(usedrift){
        // move the position out on the signed drift direction.
        tshspos += tshs._rdrift*ddir;
        herr = tshs._sderr;
      } else {
        herr = tshs._uderr;
      }
      // set the line length to be N sigma.  1 may be too short to see, TODO
      double nsigma(2.0);
      auto end1 = tshspos + nsigma*herr*ddir;
      auto end2 = tshspos - nsigma*herr*ddir;
      std::string err_title = "+/-"+std::to_string(nsigma) +"sigma";
      auto line = new REX::REveLine("TrkStrawHit Error",err_title, 1);
      line->SetNextPoint(pointmmTocm(end1.x()),pointmmTocm(end1.y()) ,pointmmTocm(end1.z()));
      line->SetNextPoint(pointmmTocm(end2.x()),pointmmTocm(end2.y()) ,pointmmTocm(end2.z()));
      //goes along that same line (ddir)
      std::string title = "TrkStrawHitSeed : x "  + std::to_string(pointmmTocm(tshspos.x()))  +  '\n'
                  + " y " + std::to_string(pointmmTocm(tshspos.y()))  +  '\n'
                  + " z " + std::to_string(pointmmTocm(tshspos.z()))  +  '\n'
                  + " time :" + std::to_string(tshs.hitTime())+  '\n'
                  + " energyDep :" + std::to_string(tshs.energyDep())+ "MeV";
      auto trkstrawpoint = new REX::REvePointSet("TrkStrawHitSeed", title,1);
      trkstrawpoint->SetMarkerStyle(REveMu2eDataInterface::mstyle);
      trkstrawpoint->SetMarkerSize(REveMu2eDataInterface::msize);
      trkstrawpoint->SetMarkerColor(drawconfig.getInt("TrkHitColor"));
      if(!usedrift)trkstrawpoint->SetMarkerColor(drawconfig.getInt("TrkNoHitColor"));
      trkstrawpoint->SetNextPoint(pointmmTocm(tshspos.x()),pointmmTocm(tshspos.y()) ,pointmmTocm(tshspos.z()));
      scene->AddElement(trkstrawpoint);
      scene->AddElement(line);
    }
  }
}

void REveMu2eDataInterface::AddTrkCaloHit(KalSeed kalseed, REX::REveElement* &scene){
  /*mu2e::TrkCaloHitSeed caloseed = kalseed.caloHit();
    art::Ptr<mu2e::CaloCluster> cluster = caloseed.caloCluster(); // TODO - use this to access x,y,z,t
    */
}

/*------------Function to color code the Tracker hits -------------*/
void REveMu2eDataInterface::AddTrkHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple,std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, REX::REveElement* &scene){
  std::cout<<"[REveMu2eDataInterface] AddTrkHits  "<<std::endl;
  std::vector<const ComboHitCollection*> combohit_list = std::get<1>(combohit_tuple);
  std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);

  GeomHandle<DetectorSystem> det;
  std::vector<StrawId> trksid(drawconfig.getInt("maxStrawID"));
  unsigned int trkhitsize=0;
  //Save the hit straw IDs of the KalSeed hits
  for(unsigned int j = 0; j< track_list.size(); j++){
    const KalSeedCollection* seedcol = track_list[j];
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k++){
        KalSeed kseed = (*seedcol)[k];
        const std::vector<mu2e::TrkStrawHitSeed> &hits = kseed.hits();
        trkhitsize = hits.size();
        for(unsigned int i = 0; i <trkhitsize; i++){
          const mu2e::TrkStrawHitSeed &hit = hits.at(i);
          trksid[i] = hit._sid;
        }
        std::vector<StrawId> usedtrksid(trkhitsize);
        std::vector<unsigned int> usedid(trkhitsize);
        //Compare the straw IDs of the Kal seed hits with the hits in the ComboHit Collection
        for(unsigned int j=0; j< combohit_list.size(); j++){
          const ComboHitCollection* chcol = combohit_list[j]; //TODO shouldnt rely on ComboHits
          if(chcol!=0){
            for(unsigned int i=0; i<chcol->size();i++){
              ComboHit hit = (*chcol)[i];
              for(unsigned int q=0; q<trkhitsize; q++){
                if(hit._sid == trksid[q]){
                  usedtrksid[q]=hit._sid;//Save the Straw ID if the KalSeed and Combo hit ID matches
                  usedid[q]=q;
                  CLHEP::Hep3Vector HitPos(hit.pos().x(), hit.pos().y(), hit.pos().z());
                  std::string chtitle = "TrkSeedHit position : x "  + std::to_string(hit.pos().x())  +  '\n'
                  + " y " + std::to_string(hit.pos().y())  +  '\n'
                  + " z " + std::to_string(hit.pos().z())  +  '\n'
                  + " time :" + std::to_string(hit.time()) +  '\n'
                  + " energy dep : "
                  + std::to_string(hit.energyDep())  +
                  + "MeV";

                  auto trkhit = new REX::REvePointSet("TrkHits", chtitle,0);
                  trkhit ->SetMarkerStyle(REveMu2eDataInterface::mstyle);
                  trkhit ->SetMarkerSize(REveMu2eDataInterface::msize);
                  // trkhit ->SetMarkerColor(drawconfig.getInt("RecoTrackColor")-4);
                  trkhit ->SetMarkerColor(drawconfig.getInt("TrkHitColor"));
                  trkhit ->SetNextPoint(pointmmTocm(HitPos.x()),pointmmTocm(HitPos.y()) ,pointmmTocm(HitPos.z()));
                  if(trkhit->GetSize() !=0 ) scene->AddElement(trkhit);
                  // std::cout<<"TrkHit = "<<HitPos.x()<<"  "<<HitPos.y()<<" "<<HitPos.z()<<std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
}

using LHPT = KalSeed::LHPT;
using CHPT = KalSeed::CHPT;
using KLPT = KalSeed::KLPT;
template<class KTRAJ> void REveMu2eDataInterface::AddKinKalTrajectory( std::unique_ptr<KTRAJ> &trajectory, REX::REveElement* &scene, unsigned int j, std::string kaltitle, double& t1, double& t2){
  t1=trajectory->range().begin();
  t2=trajectory->range().end();

  double x1=trajectory->position3(t1).x();
  double y1=trajectory->position3(t1).y();
  double z1=trajectory->position3(t1).z();

  auto line = new REX::REveLine(kaltitle,kaltitle, 100);
  line->SetPoint(0,pointmmTocm(x1), pointmmTocm(y1) , pointmmTocm(z1));
  for(double t=t1; t<=t2; t+=0.1)
  {
    const auto &p = trajectory->position3(t);
    double xt=p.x();
    double yt=p.y();
    double zt=p.z();
    line->SetNextPoint(pointmmTocm(xt), pointmmTocm(yt) , pointmmTocm(zt));
  }
  line->SetLineColor(j+6);
  line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
  scene->AddElement(line);
}

void REveMu2eDataInterface::FillKinKalTrajectory(REX::REveManager *&eveMng, bool firstloop, REX::REveElement* &scene, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, bool plotKalIntersection, bool addTrkHits, double& t1, double& t2){

  std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);
  std::vector<std::string> names = std::get<0>(track_tuple);

  for(unsigned int j=0; j< track_list.size(); j++){
    const KalSeedCollection* seedcol = track_list[j];
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k++){
        mu2e::KalSeed kseed = (*seedcol)[k];
        const std::vector<mu2e::KalSegment> &segments = kseed.segments();
        if(plotKalIntersection) AddKalIntersection(kseed, scene);
        if(kseed.loopHelixFit())
        {
          auto trajectory=kseed.loopHelixFitTrajectory();
          //auto lh = trajectory.loopHelix();
          auto lh = segments[0].loopHelix();
          std::string kaltitle = "KalSeed tag : " + names[j] +  '\n'
              + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
              + " t0 " + std::to_string(lh.t0()) + "ns " +  '\n'
              + " lam "  + std::to_string(lh.lam() ) +  '\n'
              + " rad "  + std::to_string(lh.rad() ) +  '\n'
              + " cx "  + std::to_string(lh.cx() ) +  '\n'
              + " cy "  + std::to_string(lh.cy() ) +  '\n'
              + " phi0 "  + std::to_string(lh.phi0() )+  '\n'
              + " track arrival time " + std::to_string(t1);
          AddKinKalTrajectory<LHPT>(trajectory,scene,j, kaltitle, t1, t2);
          if(addTrkHits) AddTrkStrawHit<LHPT>(kseed, scene, trajectory);
        }
        else if(kseed.centralHelixFit())
        {
          auto ch = segments[0].centralHelix();
          std::string kaltitle = "KalSeed tag : " + names[j] +  '\n'
            + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
                  + " t0 " + std::to_string(ch.t0()) + "ns " +  '\n'
            + " tandip " + std::to_string(ch.tanDip() ) +  '\n'
            + " d0 " + std::to_string(ch.d0() ) +  '\n'
            + " z0 " + std::to_string(ch.z0() ) +  '\n'
            + " phi0 " + std::to_string(ch.phi0() ) +  '\n'
            + " omega " + std::to_string(ch.omega() )+  '\n'
              + " track arrival time " + std::to_string(t1);
          auto trajectory=kseed.centralHelixFitTrajectory();
          AddKinKalTrajectory<CHPT>(trajectory,scene,j, kaltitle, t1, t2);
          if(addTrkHits) AddTrkStrawHit<CHPT>(kseed, scene, trajectory);
        }
        else if(kseed.kinematicLineFit())
        {
           auto kl = segments[0].kinematicLine();
           std::string kaltitle = "KalSeed tag : " + names[j]
              + " mom " + std::to_string(segments[0].mom()) + "MeV/c"+ '\n'
                  + " t0 " + std::to_string(kl.t0()) + "ns " +  '\n'
              + " d0 " + std::to_string(kl.d0() ) +  '\n'
              + " z0 " + std::to_string(kl.z0() ) +  '\n'
              + " phi0 " + std::to_string(kl.phi0() ) +  '\n'
              + " theta " + std::to_string(kl.theta() ) +  '\n'
              + " track arrival time " + std::to_string(t1);
          auto trajectory=kseed.kinematicLineFitTrajectory();
          AddKinKalTrajectory<KLPT>(trajectory,scene,j, kaltitle, t1, t2);
          if(addTrkHits) AddTrkStrawHit<KLPT>(kseed, scene, trajectory);
        }

      }
    }
  }
}


void REveMu2eDataInterface::AddKalSeedCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, REX::REveElement* &scene){
  //bool useKinKal = true; //TODO this will eventually be the default

  std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);
  std::vector<std::string> names = std::get<0>(track_tuple);

  for(unsigned int j=0; j< track_list.size(); j++){
    const KalSeedCollection* seedcol = track_list[j];
    //bool useKinKal = names[j].find("KK") != std::string::npos;
    //if(useKinKal) {
    //  FillKinKalTrajectory(eveMng, firstloop, scene,  track_tuple);
    //}
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k = k + 20){
        mu2e::KalSeed kseed = (*seedcol)[k];
        //if(plotKalIntersection) AddKalIntersection(kseed, scene);
        const std::vector<mu2e::KalSegment> &segments = kseed.segments();
        unsigned int nSegments=segments.size();

        if(nSegments==0) continue;
        const mu2e::KalSegment &segmentFirst = kseed.segments().front();
        const mu2e::KalSegment &segmentLast = kseed.segments().back();
        double fltLMin=segmentFirst.fmin();
        double fltLMax=segmentLast.fmax();

        for(unsigned int m = 0; m<nSegments; m++){

          const mu2e::KalSegment &segment = segments.at(m);
          fltLMin=segment.fmin();
          fltLMax=segment.fmax();
          std::string kaltitle = " tag : " + names[j] ;

          if(kseed.loopHelixFit()){
            auto lh = segment.loopHelix();
            kaltitle = "KalSeed tag : " + names[j] +  '\n'
              + " mom " + std::to_string(segment.mom()) +  '\n'
              + " MeV/c t0 " + std::to_string(lh.t0()) +  '\n'
              + " lam "  + std::to_string(lh.lam() ) +  '\n'
              + " rad "  + std::to_string(lh.rad() ) +  '\n'
              + " cx "  + std::to_string(lh.cx() ) +  '\n'
              + " cy "  + std::to_string(lh.cy() ) +  '\n'
              + " phi0 "  + std::to_string(lh.phi0() );
          }
          if(kseed.centralHelixFit() ){
            auto ch = segment.centralHelix();
            kaltitle = "KalSeed tag : " + names[j] +  '\n'
              + " mom " + std::to_string(segment.mom()) +  '\n'
              + " MeV/c t0 " + std::to_string(ch.t0()) +  '\n'
              + " tandip " + std::to_string(ch.tanDip() ) +  '\n'
              + " d0 " + std::to_string(ch.d0() ) +  '\n'
              + " z0 " + std::to_string(ch.z0() ) +  '\n'
              + " phi0 " + std::to_string(ch.phi0() ) +  '\n'
              + " omega " + std::to_string(ch.omega() );
          }
          if(kseed.kinematicLineFit() ){
            auto kl = segment.kinematicLine();
            kaltitle = "KalSeed tag : " + names[j]
              + " mom " + std::to_string(segment.mom())
              + " MeV/c t0 " + std::to_string(kl.t0()) +  '\n'
              + " d0 " + std::to_string(kl.d0() ) +  '\n'
              + " z0 " + std::to_string(kl.z0() ) +  '\n'
              + " phi0 " + std::to_string(kl.phi0() ) +  '\n'
              + " theta " + std::to_string(kl.theta() );
          }
          auto line = new REX::REveLine(kaltitle,kaltitle,1);
          if(m>0){
            double fltLMaxPrev=segments.at(m-1).fmax();
            fltLMin=(fltLMin+fltLMaxPrev)/2.0;
          }
          if(m+1<nSegments){
            double fltLMinNext=segments.at(m+1).fmin();
            fltLMax=(fltLMax+fltLMinNext)/2.0;
          }
          for(double fltL=fltLMin; fltL<=fltLMax; fltL+=1.0){
            XYZVectorF pos;
            segment.helix().position(fltL,pos);
            CLHEP::Hep3Vector p =  GenVector::Hep3Vec(pos);
            line->SetNextPoint(pointmmTocm(p.x()), pointmmTocm(p.y()) , pointmmTocm(p.z()));
          }
          int col = j+3;
          if(col > 9 ){  col += col+10; }
          line->SetLineColor(col);
          line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
          scene->AddElement(line);
        }

      }
    }
  }
  //}
}

void REveMu2eDataInterface::AddCosmicTrackFit(REX::REveManager *&eveMng, bool firstLoop_, const mu2e::CosmicTrackSeedCollection *cosmiccol, REX::REveElement* &scene){

  std::cout<<"[REveMu2eDataInterface] AddCosmicTrackSeed "<<std::endl;
  if(cosmiccol!=0){
    auto line = new REX::REveLine("Cosmic","Cosmic",2);
    for(unsigned int i=0; i< cosmiccol->size(); i++){
      mu2e::CosmicTrackSeed const  &sts= (*cosmiccol)[i];
      mu2e::CosmicTrack st = sts._track;
      double ty1 = sts._straw_chits[0].pos().y();
      double ty2 = sts._straw_chits[sts._straw_chits.size()-1].pos().y();
      double tx1 = st.MinuitParams.A0  - st.MinuitParams.A1*ty1;
      double tx2 = st.MinuitParams.A0  - st.MinuitParams.A1*ty2;
      double tz1 = st.MinuitParams.B0  - st.MinuitParams.B1*ty1;
      double tz2 = st.MinuitParams.B0  - st.MinuitParams.B1*ty2;
      line->SetNextPoint(pointmmTocm(tx1), pointmmTocm(ty1) , pointmmTocm(tz1));
      line->SetNextPoint(pointmmTocm(tx2), pointmmTocm(ty2) , pointmmTocm(tz2));
    }

    line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
    line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
    scene->AddElement(line);
  }

}

