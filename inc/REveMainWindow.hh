#ifndef _REveMainWindow_hh
#define _REveMainWindow_hh

#include <vector>
#include <string>
#include <iostream>
#include "TGeoManager.h"
#include "TClass.h"
#include "TRandom.h"
#include "TGeoTube.h"
#include "TGeoNode.h"
#include "TGeoSphere.h"
#include "TParticle.h"
#include "TApplication.h"
#include "TMatrixDSym.h"
#include "TVector.h"
#include "TMatrixDEigen.h"
#include "TTimer.h"
#include "TGeoMatrix.h"
#include "TError.h"
#include <TApplication.h>
#include <TSystem.h>
#include <TGTextEntry.h>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveUtil.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveJetCone.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/RWebDisplayArgs.hxx>
#include <ROOT/REveTrack.hxx>
#include <ROOT/REveTrackPropagator.hxx>
#include <ROOT/REveEllipsoid.hxx>

#include "art/Framework/Principal/Event.h"
#include "REve/inc/DataCollections.hh"
#include "REve/inc/REveMu2eDataInterface.hh"
#include "REve/inc/REveMu2eMCInterface.hh"
#include "REve/inc/REveGeomUtil.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
namespace REX = ROOT::Experimental;
 
namespace mu2e {
    
    
  struct GeomOptions{
      // geom options
      bool showCRV = false;
      bool showPS = false;
      bool showTS = false; 
      bool showDS = false;
      bool show2D = true;
      bool caloVST = false;
      bool showST = true;
      bool extracted = false;
      bool showSTM = false;
      bool showCalo = true;
      bool showTracker = true;
      GeomOptions(){};
      GeomOptions(bool crv, bool ps, bool ts, bool ds, bool twodim, bool cVST, bool st, bool ext, bool stm, bool calo, bool trk) 
      : showCRV(crv), showPS(ps), showTS(ts), showDS(ds), show2D(twodim), caloVST(cVST), showST(st), extracted(ext), showSTM(stm), showCalo(calo), showTracker(trk) {};
      void fill(bool crv, bool ps, bool ts, bool ds, bool twodim, bool cVST, bool st, bool ext, bool stm, bool cal, bool trk) {
        showCRV = (crv);
        showPS = (ps);
        showTS = (ts);
        showDS = (ds);
        show2D = (twodim);
        caloVST = (cVST);
        showST = (st);
        extracted = (ext);
        showSTM = (stm);
        showCalo = (cal);
        showTracker = (trk);
      }
      void print(){
        std::cout<<"***** Geom Options ****** "<<'\n'
        <<" show CRV : "<<showCRV <<'\n'
        <<" show PS : "<<showPS <<'\n'
        <<" show TS : "<<showTS <<'\n'
        <<" show DS : "<<showDS <<'\n'
        <<" show 2D : "<<show2D <<'\n'
        <<" show ST : "<<showST <<'\n'
        <<" show STM : "<<showSTM <<'\n'
        <<" show Calo : "<<showCalo <<'\n'
        <<" show Trk : "<<showTracker <<'\n'
        <<" show Extracted : "<<extracted <<'\n'
        <<"************************ "<<std::endl;
      }
     };
     
  struct DrawOptions{
      // data options
      bool addCosmicTracks = false;
      bool addHelices = false;
      bool addTracks = false;
      bool addClusters = false; 
      bool addComboHits = false;
      bool addCRVInfo = false;
      bool addCRVClusters = false;
      bool addTimeClusters = false;
      bool addTrkHits = false;
      bool addMCTrajectories = false;
      bool addTrkErrBar = true;
      bool addCrystalDraw = false;
      DrawOptions(){};
      DrawOptions(bool cosmictracks, bool helices, bool tracks, bool clusters, bool combohits, bool crv, bool crvclu, bool timeclusters, bool trkhits, bool mctraj, bool errbar, bool crys) 
      : addCosmicTracks(cosmictracks), addHelices(helices), addTracks(tracks), addClusters(clusters), addComboHits(combohits), addCRVInfo(crv), addCRVClusters(crvclu), addTimeClusters(timeclusters), addTrkHits(trkhits), addMCTrajectories(mctraj), addTrkErrBar(errbar), addCrystalDraw(crys) {};
     };
     
    class REveMainWindow  : public REX::REveElement {

        public :
            explicit REveMainWindow() { SetErrorHandler(DefaultErrorHandler); }
            virtual ~REveMainWindow() {}
            #ifndef __CINT__
            REveMu2eDataInterface *pass_data;
            REveMu2eMCInterface *pass_mc;

            void makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int j, bool crys1, bool crys2, std::string name, int color);
            void changeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int j, bool crys1, bool crys2, std::string name);
            void showNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal, std::vector<double> shift, bool print, bool single, int color);
            void changeNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal, std::vector<double> shift, bool print, bool single);
            void GeomDrawer(TGeoNode* node, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, GeomOptions geomOpts); 
            void makeGeometryScene(REX::REveManager *eveMng,  GeomOptions geomOpts, std::string filename);
            void showEvents(REX::REveManager *eveMng,  REX::REveElement* &eventScene, bool firstLoop, DataCollections &data, DrawOptions drawOpts, std::vector<int> particleIds, bool strawdisplay, GeomOptions geomOpts);
            void changeEveGeoShape(TGeoNode* node, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level);
            void createProjectionStuff(REX::REveManager *eveMng, GeomOptions geomOpt);
            void AddTrackerProjection(REX::REveManager *eveMng);
            void projectScenes(REX::REveManager *eveMng, bool geomp, bool eventp);
            void projectEvents(REX::REveManager *eveMng, GeomOptions geomOpt);
            REX::REveProjectionManager *mngTrackerXY = nullptr;
            REX::REveProjectionManager *mngXYCaloDisk1 = nullptr;
            REX::REveProjectionManager *mngXYCaloDisk2 = nullptr;
            REX::REveProjectionManager *mngRhoZ   = nullptr;
            REX::REveProjectionManager *mngCRVXY = nullptr;
            REX::REveProjectionManager *mngCRVRhoZ   = nullptr;
            REX::REveScene  *TrackerXYGeomScene = nullptr, *TrackerXYEventScene = nullptr;
            REX::REveScene  *XYCaloDisk1GeomScene = nullptr, *XYCaloDisk1EventScene = nullptr;
            REX::REveScene  *XYCaloDisk2GeomScene = nullptr, *XYCaloDisk2EventScene = nullptr;
            REX::REveScene  *CRVXYGeomScene = nullptr, *CRVXYEventScene = nullptr;
            REX::REveScene  *CRVrhoZGeomScene = nullptr, *CRVrhoZEventScene = nullptr;
            REX::REveScene  *rhoZGeomScene = nullptr, *rhoZEventScene = nullptr;
            REX::REveViewer *TrackerXYView = nullptr;
            REX::REveViewer *XYCaloDisk1View = nullptr;
            REX::REveViewer *XYCaloDisk2View = nullptr;
            REX::REveViewer *rhoZView = nullptr;
            REX::REveViewer *CRVXYView = nullptr;
            REX::REveViewer *CRVrhoZView = nullptr;
            std::vector<int> crystals_hit;
            
            #else
                ClassDef(REveMainWindow, 0);
            #endif

    };

}
#endif
