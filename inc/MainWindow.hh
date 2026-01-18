#ifndef _MainWindow_hh
#define _MainWindow_hh

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
#include <ROOT/REveTableInfo.hxx>
#include <ROOT/REveViewContext.hxx>
#include "art/Framework/Principal/Event.h"
#include "EventDisplay/inc/DataCollections.hh"
#include "EventDisplay/inc/DataInterface.hh"
#include "EventDisplay/inc/MCInterface.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"

#include <utility>
namespace REX = ROOT::Experimental;

namespace mu2e {


  struct GeomOptions{
      // geom options
      bool showCrv = false;
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
      bool showCaloCrystals = true;
      bool showEM = false;
      GeomOptions(){};
      GeomOptions(bool crv, bool ps, bool ts, bool ds, bool twodim, bool cVST, bool st, bool ext, bool stm, bool calo, bool trk, bool crys, bool em)
      : showCrv(crv), showPS(ps), showTS(ts), showDS(ds), show2D(twodim), caloVST(cVST), showST(st), extracted(ext), showSTM(stm), showCalo(calo), showTracker(trk), showCaloCrystals(crys), showEM(em) {};
      void fill(bool crv, bool ps, bool ts, bool ds, bool twodim, bool cVST, bool st, bool ext, bool stm, bool cal, bool trk, bool crys, bool em) {
        showCrv = (crv);
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
        showCaloCrystals = (crys);
        showEM = (em);
      }
      void print(){
        std::cout<<"***** Geom Options ****** "<<'\n'
        <<" show Crv : "<<showCrv <<'\n'
        <<" show PS : "<<showPS <<'\n'
        <<" show TS : "<<showTS <<'\n'
        <<" show DS : "<<showDS <<'\n'
        <<" show 2D : "<<show2D <<'\n'
        <<" show ST : "<<showST <<'\n'
        <<" show STM : "<<showSTM <<'\n'
        <<" show Calo : "<<showCalo <<'\n'
        <<" show Trk : "<<showTracker <<'\n'
        <<" show crystals : "<<showCaloCrystals <<'\n'
        <<" show Extracted : "<<extracted <<'\n'
         <<" show Ext Mon : "<<showEM <<'\n'
        <<"************************ "<<std::endl;
      }
     };

  struct KinKalOptions{
    bool addKalInter = false;
    bool addTrkStrawHits = false;
    bool addTrkCaloHits = false;
    KinKalOptions(){};
    KinKalOptions(bool kalinter, bool trkstrawhits, bool trkcalohits)
      : addKalInter(kalinter), addTrkStrawHits(trkstrawhits), addTrkCaloHits(trkcalohits){};

  };


  struct DrawOptions{
      // data options
      bool addCosmicTracks = false;
      bool addHelices = false;
      bool addTracks = false;
      bool addCaloDigis = false;
      bool addClusters = false;
      bool addComboHits = false;

      bool addBkgClusters = false;
      bool addCrvRecoPulse = false;
      bool addCrvClusters = false;

      bool addTimeClusters = false;
      bool addTrkHits = false; // legacy
      bool addMCTrajectories = false;
      bool addSurfaceSteps = false;
      bool addSimParts = false;
      bool addTrkErrBar = true;
      bool addCrystalDraw = false;
      bool addCrvBars = true;
      DrawOptions(){};

      DrawOptions(bool cosmictracks, bool helices, bool tracks, bool calodigis, bool clusters, bool combohits, bool bkgclusters, bool crv, bool crvclu, bool timeclusters, bool trkhits, bool mctraj, bool surfsteps, bool simparts, bool errbar, bool crys, bool crvbars)
      : addCosmicTracks(cosmictracks), addHelices(helices), addTracks(tracks), addCaloDigis(calodigis), addClusters(clusters), addComboHits(combohits), addBkgClusters(bkgclusters), addCrvRecoPulse(crv), addCrvClusters(crvclu), addTimeClusters(timeclusters), addTrkHits(trkhits), addMCTrajectories(mctraj), addSurfaceSteps(surfsteps), addSimParts(simparts), addTrkErrBar(errbar), addCrystalDraw(crys), addCrvBars(crvbars) {};

     };

    class MainWindow  : public REX::REveElement {

        public :
            explicit MainWindow() { SetErrorHandler(DefaultErrorHandler); }
            virtual ~MainWindow() {}
            #ifndef __CINT__
            DataInterface *pass_data;
            MCInterface *pass_mc;

            void makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int j, bool crys1, bool crys2, std::string name, int color);
            void changeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int j, bool crys1, bool crys2, std::string name);
            void showNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal, std::vector<double> shift, bool print, bool single, int color);
            void getOffsets(TGeoNode* n, const std::string& str, REX::REveTrans& trans, int maxlevel, int level, std::vector<std::pair<std::string, std::vector<float>>> & offsets);
            void changeNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal, std::vector<double> shift, bool print, bool single);
            void GeomDrawer(TGeoNode* node, REX::REveTrans& trans,  REX::REveElement* bholder, REX::REveElement* trholder, REX::REveElement* cholder,REX::REveElement* crholder, REX::REveElement* vholder, REX::REveElement* tholder,int maxlevel, int level, GeomOptions geomOpts);

            void GeomDrawerNominal(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* bholder, REX::REveElement* trholder, REX::REveElement* cholder,REX::REveElement* crholder, REX::REveElement* vholder, REX::REveElement* tholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets);
            void GeomDrawerSol(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets);
            void GeomDrawerExtracted(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* bholder, REX::REveElement* trholder, REX::REveElement* cholder,REX::REveElement* crholder, REX::REveElement* vholder, REX::REveElement* tholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets);
            void makeGeometryScene(REX::REveManager *eveMng,  GeomOptions geomOpts, std::string filename);
            void showEvents(REX::REveManager *eveMng,  REX::REveElement* &eventScene, bool firstLoop,  bool firstLoopCalo, DataCollections &data, DrawOptions drawOpts, std::vector<int> particleIds, bool strawdisplay, GeomOptions geomOpts, KinKalOptions KKOpts);
            void changeEveGeoShape(TGeoNode* node, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level);
            void createProjectionStuff(REX::REveManager *eveMng);
            void AddTrackerProjection(REX::REveManager *eveMng);
            void projectScenes(REX::REveManager *eveMng, bool geomp, bool eventp);
            void projectEvents(REX::REveManager *eveMng);
            void maketable(REX::REveManager *eveMng);

            REX::REveProjectionManager *mngTrackerXY = nullptr;
            REX::REveProjectionManager *mngXYCaloDisk0 = nullptr;
            REX::REveProjectionManager *mngXYCaloDisk1 = nullptr;
            REX::REveProjectionManager *mngRhoZ   = nullptr;
            REX::REveScene  *TrackerXYGeomScene = nullptr, *TrackerXYEventScene = nullptr;
            REX::REveScene  *XYCaloDisk0GeomScene = nullptr, *XYCaloDisk0EventScene = nullptr;
            REX::REveScene  *XYCaloDisk1GeomScene = nullptr, *XYCaloDisk1EventScene = nullptr;
            REX::REveScene  *rhoZGeomScene = nullptr, *rhoZEventScene = nullptr;
            REX::REveViewer *TrackerXYView = nullptr;
            REX::REveViewer *XYCaloDisk0View = nullptr;
            REX::REveViewer *XYCaloDisk1View = nullptr;
            REX::REveViewer *rhoZView = nullptr;

            #else
                ClassDef(MainWindow, 0);
            #endif

    };

}
#endif
