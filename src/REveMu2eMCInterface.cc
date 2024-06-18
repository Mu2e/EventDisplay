#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "REve/inc/REveMu2eMCInterface.hh"
using namespace mu2e;
namespace REX = ROOT::Experimental;

  int REveMu2eMCInterface::Contains(std::vector<int> v, int x)
  {
    return std::count(v.begin(), v.end(), abs(x));
  }
  
  const char* REveMu2eMCInterface::GetParticleName(int PDGCode){ 
    const char* pid = "pid";
    switch(PDGCode) {
        case PDGCode::e_minus:
            pid = "electron -";
            break;
        case PDGCode::e_plus:
            pid = "positron +";
            break;
        case PDGCode::mu_minus:
            pid = "muon - ";
            break;
        case PDGCode::mu_plus:
            pid = "muon + ";
            break;
        case PDGCode::pi_minus:
            pid = "pion -";
            break;
        case PDGCode::pi_plus:
            pid = "pion +";
            break;
        case PDGCode::pi0:
            pid = "pion 0";
            break;
        case PDGCode::proton:
            pid = "proton";
            break;
        case PDGCode::n0:
            pid = "neutron";
            break;
        case PDGCode::gamma:
            pid = "gamma";
            break;
        default:
            pid = "other";
            break;
        }
    return pid;
  }
  
  void REveMu2eMCInterface::SetLineColorPID(int PDGCode,REX::REveLine *line){
    Color_t color;
    switch(PDGCode) {
      case PDGCode::e_minus:
        color = kRed;
        break;
      case PDGCode::e_plus:
        color = kGreen;
        break;
      case PDGCode::mu_minus:
        color = kBlack;
        break;
      case PDGCode::mu_plus:
        color = kViolet;
        break;
      case PDGCode::pi_minus:
        color = kMagenta;
        break;
      case PDGCode::pi_plus:
        color = kRed-7;
        break;
      case PDGCode::pi0:
        color = kGreen-7;
        break;
      case PDGCode::proton:
        color = kBlue;
        break;
      case PDGCode::n0:
        color = kViolet-2;
        break;
      case PDGCode::gamma:
        color = kOrange;
        break;
      default:
        color = kCyan-4;
        break;
    }
    line->SetLineColor(color);
  }
  
  
  void REveMu2eMCInterface::AddMCTrajectoryCollection(REX::REveManager *&eveMng, bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const MCTrajectoryCollection *>> mctrack_tuple, REX::REveElement* &scene, std::vector<int> particleIds, bool extracted){
  
    std::string drawfilename("REve/config/drawutils.txt");
    SimpleConfig drawconfig(drawfilename);
    
    // eEtract the track and input tag:
    std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(mctrack_tuple);
    std::vector<std::string> names = std::get<0>(mctrack_tuple);
    
    // Loop over tracks:
    for(unsigned int j=0; j< track_list.size(); j++){
      const MCTrajectoryCollection* trajcol = track_list[j];
      
      if(trajcol!=0){
          std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
          for(unsigned int k = 0; k < trajcol->size(); k++){ 

            for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
            { 
              // Check user defined list of particles to plot
              int x = Contains(particleIds,trajectoryIter->first->pdgId()); 
              // get particle name
              const char* particlename = GetParticleName(trajectoryIter->first->pdgId());
              if(x == 1){
                const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
                // Make label
                std::string energy = std::to_string(points[0].kineticEnergy());
                std::string mctitle = " MCTrajectory tag : particle = " + std::string(particlename)  +  '\n'
                + " energy = " + energy + "MeV" +  '\n'
                + " creation code = " + std::to_string(trajectoryIter->first->creationCode()) +  '\n'
                + " stopping code = " + std::to_string(trajectoryIter->first->stoppingCode())  +  '\n'
                + " end global time = " + std::to_string(trajectoryIter->first->endGlobalTime())
                + " ns";
                // create line with the above label
                auto line = new REX::REveLine(mctitle,mctitle,1);
                // add points
                for(unsigned int i=0; i < points.size();i++){
                  CLHEP::Hep3Vector Pos(points[i].x(), points[i].y(), points[i].z());
                    GeomHandle<DetectorSystem> det;
                    CLHEP::Hep3Vector HitPos = det->toDetector(Pos);
                    if(pointmmTocm(HitPos.y()) < 800){ // a reasonable height above the CRV
                      line->SetNextPoint(pointmmTocm((HitPos.x())),pointmmTocm((HitPos.y())),pointmmTocm(HitPos.z()));
                    }
                }
                // set line colour
                SetLineColorPID(trajectoryIter->first->pdgId(),line);
                line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
                scene->AddElement(line); 
              } else std::cout<<"Warning: No Particles of User-Specified Type In File "<<std::endl;
            }
          }
        }
      }
   }



