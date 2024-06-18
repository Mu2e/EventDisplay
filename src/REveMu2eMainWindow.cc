#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "REve/inc/REveMu2eMainWindow.hh"

namespace REX = ROOT::Experimental;
using namespace std;
using namespace mu2e;

// Get geom information from configs:
std::string calfilename("REve/config/geomutils.txt");
std::string drawoptfilename("REve/config/drawutils.txt");
SimpleConfig geomconfig(calfilename);
SimpleConfig drawconfigf(drawoptfilename);

double disk1_center = 0;
double disk2_center = 0;
double nCrystals = 674;

SimpleConfig TConfig("Offline/Mu2eG4/geom/tracker_v6.txt");
double trackerz0 = TConfig.getDouble("tracker.z0");//10171
SimpleConfig ExConfig("Offline/Mu2eG4/geom/geom_common_extracted.txt");
double trackerz0_extracted = ExConfig.getDouble("tracker.z0");//24175
double motherhalflength = TConfig.getDouble("tracker.mother.halfLength",motherhalflength);

SimpleConfig CConfig("Offline/Mu2eG4/geom/calorimeter_CsI.txt");
double caloz0 = CConfig.getDouble("calorimeter.caloMotherZ0");//11842;
double caloz1 = CConfig.getDouble("calorimeter.caloMotherZ1");//13220;
double dz = (caloz1 - caloz0)/2 + caloz0 - trackerz0;

SimpleConfig SConfig("Offline/Mu2eG4/geom/stoppingTarget_CD3C_34foils.txt");
double STz0 = SConfig.getDouble("stoppingTarget.z0InMu2e");
double STz = STz0 - motherhalflength;//4236

SimpleConfig PTConfig("Offline/Mu2eG4/geom/ProductionTargetInPS.txt");
SimpleConfig PSConfig("Offline/Mu2eG4/geom/ProductionSolenoid_v02.txt");
SimpleConfig PTHConfig("Offline/Mu2eG4/geom/ProductionTarget_Hayman_v2_0.txt");
double PScoilRefZ  = PSConfig.getDouble("PS.coilRefZ"); // from PS
double PTz0 = PTConfig.getDouble("productionTarget.zNominal");
double PTHL = 120;//PTHConfig.getDouble("targetPS_halfHaymanLength ");//120;
double PTz = PTz0;// - PTHL;
double FrontTracker_gdmltag;

void REveMu2eMainWindow::makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int val, bool crystal1, bool crystal2, std::string name, int color)
{
  bool isMother = name.find("ProductionTargetMother") != string::npos;
  auto gss = n->GetVolume()->GetShape();
  auto b1s = new REX::REveGeoShape(n->GetName());

  b1s->InitMainTrans();
  b1s->RefMainTrans().SetFrom(trans.Array());
  b1s->SetShape(gss);

  int transpar = drawconfigf.getInt("BLtrans");
  if(name.find("CRS") != string::npos) transpar = drawconfigf.getInt("CRVtrans");
  if(name.find("Tracker") != string::npos) transpar = drawconfigf.getInt("TRKtrans") ;
  if(name.find("caloDisk") != string::npos) transpar = drawconfigf.getInt("CALtrans") ;
  if(name.find("caloCrystal") != string::npos) transpar = drawconfigf.getInt("CRYtrans") ;

  b1s->SetMainTransparency(transpar);
  b1s->SetMainColor(color);
  b1s->SetEditMainColor(true);
  b1s-> SetEditMainTransparency(true);
  if(!isMother) holder->AddElement(b1s);

  // make 2D projections
  if( crystal1 ){ // add crystals to correct disk
    mngXYCaloDisk0->ImportElements(b1s, XYCaloDisk0GeomScene);
  }if( crystal2 ){
    mngXYCaloDisk1->ImportElements(b1s, XYCaloDisk1GeomScene);
  }

  //show only first tracker plane
  if( val == FrontTracker_gdmltag ){
    mngTrackerXY->ImportElements(b1s, TrackerXYGeomScene); //shows only one plane for 2D view for simplicity
  }

  bool isCRV = name.find("CRS") != string::npos;
  // remove CRV from the YZ view as it blocks tracker/calo view (will have its own view)
  if(!isCRV) { mngRhoZ->ImportElements(b1s, rhoZGeomScene); }
}

int j = 0;
int fp =0;

void REveMu2eMainWindow::showNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal,  std::vector<double> shift, bool print, bool single, int color) {
  ++level;
  if (level > maxlevel){
    return;
  }
  std::string name(n->GetName());
  j++;
  if(print) std::cout<<j<<" "<<name<<std::endl;
  bool cry1 = false; // these help us know which disk and to draw crystals
  bool cry2 = false;
  int ndau = n->GetNdaughters();
  if (name.find(str)!= std::string::npos and single){
    REX::REveTrans ctrans;
    ctrans.SetFrom(trans.Array());
    {
      TGeoMatrix     *gm = n->GetMatrix();
      const Double_t *rm = gm->GetRotationMatrix();
      const Double_t *tv = gm->GetTranslation();
      REX::REveTrans t;
      t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
      t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
      t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
      t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1]  + shift[1]; t(3,4) = tv[2] + shift[2];
      ctrans *= t;
    }
    n->ls();
    makeEveGeoShape(n, ctrans, holder, j, cry1, cry2, name, color);
  }
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * pn = n->GetDaughter(i);
    if (name.find(str)!= std::string::npos and i==0 ){ //To make the geometry translucant we only add i==0
      REX::REveTrans ctrans;
      ctrans.SetFrom(trans.Array());
      {
        TGeoMatrix     *gm = n->GetMatrix();
        const Double_t *rm = gm->GetRotationMatrix();
        const Double_t *tv = gm->GetTranslation();
        REX::REveTrans t;
        t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
        t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
        t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
        t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1]  + shift[1]; t(3,4) = tv[2] + shift[2];
        //std::cout<<name<<"  "<<tv[0] + shift[0]<<" "<<tv[1]  + shift[1] << " "<< tv[2] + shift[2]<<std::endl;

        if(name.find("TrackerPlaneEnvelope_00") != string::npos) {
          FrontTracker_gdmltag = j;
        }
        if(name.find("caloDisk_00") != string::npos) {
          disk1_center = tv[2] ;
        }
        if(name.find("caloDisk_10") != string::npos) {
          disk2_center = tv[2] ;
        }
        if(caloshift){
          fp++;
          double d = 0;
          if(fp < nCrystals) { d = disk1_center; }
          else { d = disk2_center; }
          if(crystal) { t(1,4) = tv[0]; t(2,4) = tv[1] ; t(3,4) = tv[2] + dz/10 + d; }
          else { t(1,4) = tv[0]; t(2,4) = tv[1]; t(3,4) = tv[2] + dz/10; }
          if (fp < nCrystals and crystal) cry1 = true;
          if (fp >= nCrystals and crystal) cry2 = true;
        }
        ctrans *= t;
      }
      n->ls();
      makeEveGeoShape(n, ctrans, holder, j, cry1, cry2, name, color);
    }
    showNodesByName( pn, str, onOff, _diagLevel, trans,  holder, maxlevel, level, caloshift, crystal, shift, print, single, color);
  }
}
/* TODO this is a work in progress
   void REveMu2eMainWindow::maketable(REX::REveManager *eveMng){
   REX::REveViewContext          *m_viewContext {nullptr};
   REX::REveTableViewInfo *tableInfo = new REX::REveTableViewInfo();

   tableInfo->table("RecHit").
   column("pt",     1, "i.fPt");

   m_viewContext->SetTableViewInfo(tableInfo);
   for (auto &c : eveMng->GetScenes()->RefChildren()) {
   if (!strncmp(c->GetCName(),"Table", 5))
   c->AddElement(m_viewContext->GetTableViewInfo());
//m_collections = eveMng->SpawnNewScene("Collections", "Collections");
}
}
*/
/* function to hide all elements which are not PS,TS, DS */
void REveMu2eMainWindow::GeomDrawer(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, REX::REveElement* trackerholder, REX::REveElement* caloholder, REX::REveElement* crystalsholder, REX::REveElement* crvholder, REX::REveElement* targetholder,int maxlevel, int level, GeomOptions geomOpt){

  std::vector<double> shift;

  shift.push_back(0);
  shift.push_back(0);
  shift.push_back(0);
  //geomOpt.print();

  if(geomOpt.showEM){
    static std::vector <std::string> substrings_em  {"ExtMon"};
    shift.at(0) = geomconfig.getDouble("psts_x")/10;
    shift.at(1) = geomconfig.getDouble("psts_y")/10;
    shift.at(2) = geomconfig.getDouble("psts_z")/10;
    for(auto& i: substrings_em){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }
    static std::vector <std::string> substrings_ps  {"ProductionTarget"};
    shift.at(0) = 780.85; //GDML
    shift.at(1) = -0.06; //GDML
    shift.at(2) = 1*PTz/10 -1*trackerz0/10 + PTHL/10;
    for(auto& i: substrings_ps){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showTS){
    static std::vector <std::string> substrings_ts  {"TS"};
    shift.at(0) = geomconfig.getDouble("psts_x")/10; //TODO these numbers need to be better defined
    shift.at(1) = geomconfig.getDouble("psts_y")/10;
    shift.at(2) = geomconfig.getDouble("psts_z")/10;
    for(auto& i: substrings_ts){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level,  false, false, shift, false, false, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showPS){
    static std::vector <std::string> substrings_ps  {"PSVacuum"};
    shift.at(0) = geomconfig.getDouble("psts_x")/10;
    shift.at(1) = geomconfig.getDouble("psts_y")/10;
    shift.at(2) = geomconfig.getDouble("psts_z")/10;
    for(auto& i: substrings_ps){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level,  false, false, shift, false, false, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showPS){
    static std::vector <std::string> substrings_ps  {"ProductionTarget"};
    shift.at(0) = 780.85; //GDML
    shift.at(1) = -0.06; //GDML
    shift.at(2) = 1*PTz/10 -1*trackerz0/10 + PTHL/10;
    for(auto& i: substrings_ps){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showDS){
    static std::vector <std::string> substrings_ds {"DS1Vacuum","DS2Vacuum","DS3Vacuum","StoppingTarget_Al"};
    for(auto& i: substrings_ds){
      shift.at(0) = geomconfig.getDouble("psts_x")/10;
      shift.at(1) = geomconfig.getDouble("psts_y")/10;
      shift.at(2) = geomconfig.getDouble("psts_z")/10;
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, false, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showCRV and !geomOpt.extracted){
    static std::vector <std::string> substrings_crv {"CRSmother_CRV"};
    shift.at(0) = 481.75;//from GDML, look at CRS layer 18, center is where 0 should be.
    //6.45 to shift by comparing to Offline, need to understand.
    shift.at(1) = 585.31  + 6.45;//from GDML, look at layer 16, 17 these should be centered at 0
    shift.at(2) = -431   + 6.45 ;//1280.07-849 first is from GDML latter is from Offline
    for(auto& i: substrings_crv){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, false, drawconfigf.getInt("CRVColor"));
    }
  }
  if(geomOpt.showCRV and geomOpt.extracted){
    SimpleConfig CRVConfig("Offline/Mu2eG4/geom/crv_counters_extracted_v01.txt");
    std::vector<double> firstCounterEX;
    std::vector<double> firstCounterT1;
    std::vector<double> firstCounterT2;
    CRVConfig.getVectorDouble("crs.firstCounterEX", firstCounterEX);
    CRVConfig.getVectorDouble("crs.firstCounterT1", firstCounterT1);
    CRVConfig.getVectorDouble("crs.firstCounterT2", firstCounterT2);

    static std::vector <std::string> substrings_ex {"CRSscintLayer_0","CRSmotherLayer_CRV_EX"};//
    shift.at(0) = 0;
    shift.at(1) = firstCounterEX[1]/10;
    shift.at(2) =  82.8/2 + 124.2 + firstCounterEX[2]/10 - trackerz0_extracted/10; //tracker - first counter pos plus GDML position of first + wisth of module
    //std::cout<<"Position of CRV "<<firstCounterEX[2] - trackerz0_extracted<<std::endl;
    for(auto& i: substrings_ex){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }

    static std::vector <std::string> substrings_t1  {"CRSscintLayer_1","CRSmotherLayer_CRV_T1"};//"CRSScintillatorBar_1_1",
    shift.at(0) = 0;
    shift.at(1) = firstCounterT1[1]/10;
    shift.at(2) = firstCounterT1[2]/10 - trackerz0_extracted/10; //tracker - first counter pos
    for(auto& i: substrings_t1){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }

    static std::vector <std::string> substrings_t2  {"CRSscintLayer_2","CRSmotherLayer_CRV_T2"};//"CRSScintillatorBar_1_2",note for the bar, maybe first number changes?
    shift.at(0) = 0;
    shift.at(1) = firstCounterT2[1]/10;
    shift.at(2) = 82.8/2 +41.4 + firstCounterT2[2]/10 - trackerz0_extracted/10; //tracker - first counter pos
    for(auto& i: substrings_t2){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));

    }
  }
  if(geomOpt.showCalo){
    static std::vector <std::string> substrings_disk  {"caloDisk"};
    for(auto& i: substrings_disk){
      shift.at(0) = 0;
      shift.at(1) = 0;
      shift.at(2) = 0;
      showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, false, false, drawconfigf.getInt("CALColor") );
    }

    if(geomOpt.showCaloCrystals){
      static std::vector <std::string> substrings_crystals  {"caloCrystal"};
      for(auto& i: substrings_crystals){
        shift.at(0) = 0;
        shift.at(1) = 0;
        shift.at(2) = 0;
        showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, false, false, drawconfigf.getInt("CALColor"));
      }
    }
  }
  if(geomOpt.showTracker){
    static std::vector <std::string> substring_tracker  {"TrackerPlaneEnvelope"};
    for(auto& i: substring_tracker){
      shift.at(0) = 0;
      shift.at(1) = 0;
      shift.at(2) = 0 ;
      showNodesByName(node,i,kFALSE, 0, trans, trackerholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("TRKColor"));
    }
  }if(geomOpt.showST and !geomOpt.extracted){
    static std::vector <std::string> substrings_stoppingtarget  {"TargetFoil"};
    shift.at(0) = 0;
    shift.at(1) = 0;
    shift.at(2) = -1*STz/10;
    for(auto& i: substrings_stoppingtarget){
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
    }
  }
  if(geomOpt.showST and !geomOpt.extracted){
    static std::vector <std::string> substrings_stoppingtarget  {"protonabs"};
    shift.at(0) = 0;
    shift.at(1) = 0;
    shift.at(2) = -1*STz/10;
    for(auto& i: substrings_stoppingtarget){
      showNodesByName(node,i,kFALSE, 0, trans, targetholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
    }
  }

  if(geomOpt.showSTM){
    static std::vector <std::string> substrings_stm  {"stmDet1Can","stmDet1","stmDet2Can","stmDet2"};
    for(auto& i: substrings_stm){
      shift.at(0) = 0;
      shift.at(1) = 0;
      shift.at(2) = 0;
      showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, false, drawconfigf.getInt("BLColor"));
    }
  }
}

void REveMu2eMainWindow::projectEvents(REX::REveManager *eveMng)
{
  for (auto &ie : eveMng->GetEventScene()->RefChildren())
  {
    TrackerXYView->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    XYCaloDisk0View->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    XYCaloDisk1View->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    rhoZView->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    mngTrackerXY->ImportElements(ie, TrackerXYEventScene);

    mngRhoZ  ->ImportElements(ie, rhoZEventScene);

    if(ie->GetName().find("disk0") != string::npos ) mngXYCaloDisk0->ImportElements(ie, XYCaloDisk0EventScene);
    if(ie->GetName().find("disk1") != string::npos ) mngXYCaloDisk1->ImportElements(ie, XYCaloDisk1EventScene);
  }
}

void REveMu2eMainWindow::createProjectionStuff(REX::REveManager *eveMng)
{
  // -------------------Tracker XY View ----------------------------------
  TrackerXYGeomScene  = eveMng->SpawnNewScene("TrackerXY Geometry","TrackerXY");
  TrackerXYEventScene = eveMng->SpawnNewScene("TrackerXY Event Data","TrackerXY");

  mngTrackerXY = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  TrackerXYView = eveMng->SpawnNewViewer("TrackerXY View", "");
  TrackerXYView->AddScene(TrackerXYGeomScene);
  TrackerXYView->AddScene(TrackerXYEventScene);

  // --------------------Tracker + Calo YZ View ------------------------------

  rhoZGeomScene  = eveMng->SpawnNewScene("ZY Detector Geometry", "ZY");
  rhoZEventScene = eveMng->SpawnNewScene("ZY Event Data","YZ");

  mngRhoZ = new REX::REveProjectionManager(REX::REveProjection::kPT_ZY );

  rhoZView = eveMng->SpawnNewViewer("ZY Detector View", "");
  rhoZView->AddScene(rhoZGeomScene);
  rhoZView->AddScene(rhoZEventScene);

  // ---------------------Calo Disk 1 XY View ----------------------------

  XYCaloDisk0GeomScene  = eveMng->SpawnNewScene("XYCaloDisk0 Geometry", "XYCaloDisk0");
  XYCaloDisk0EventScene = eveMng->SpawnNewScene("XYCaloDisk0 Event Data","XYCaloDisk0");

  mngXYCaloDisk0 = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  XYCaloDisk0View = eveMng->SpawnNewViewer("XYCaloDisk0 View", "");
  XYCaloDisk0View->AddScene(XYCaloDisk0GeomScene);
  XYCaloDisk0View->AddScene(XYCaloDisk0EventScene);

  // -------------------- Calo Disk 2 XY View ------------------------

  XYCaloDisk1GeomScene  = eveMng->SpawnNewScene("XYCaloDisk1 Geometry", "XYCaloDisk1");
  XYCaloDisk1EventScene = eveMng->SpawnNewScene("XYCaloDisk1 Event Data","XYCaloDisk1");

  mngXYCaloDisk1 = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  XYCaloDisk1View = eveMng->SpawnNewViewer("XYCaloDisk1 View", "");
  XYCaloDisk1View->AddScene(XYCaloDisk1GeomScene);
  XYCaloDisk1View->AddScene(XYCaloDisk1EventScene);
}


void REveMu2eMainWindow::showEvents(REX::REveManager *eveMng, REX::REveElement* &eventScene, bool firstLoop, bool firstLoopCalo, DataCollections &data, DrawOptions drawOpts, std::vector<int> particleIds, bool strawdisplay, GeomOptions geomOpts, KinKalOptions KKOpts){
  if(!firstLoop){
    eventScene->DestroyElements();
  }
  //...addReco:
  // for start and end of track
  double t1 = 0.;
  double t2 = 1696.;
  std::vector<const KalSeedPtrCollection*> track_list = std::get<1>(data.track_tuple);
  if(drawOpts.addTracks and track_list.size() !=0) {
    if(drawOpts.useBTrk) { pass_data->AddKalSeedPtrCollection(eveMng, firstLoop, data.track_tuple, eventScene); }
    if(!drawOpts.useBTrk){ pass_data->FillKinKalTrajectory(eveMng, firstLoop, eventScene, data.track_tuple, KKOpts.addKalInter,  KKOpts.addTrkStrawHits, t1, t2); }
  }
  if(drawOpts.addComboHits) {
    std::vector<const ComboHitCollection*> combohit_list = std::get<1>(data.combohit_tuple);
    if(combohit_list.size() !=0 ) pass_data->AddComboHits(eveMng, firstLoop, data.combohit_tuple, eventScene, strawdisplay, drawOpts.addTrkErrBar);
  }

  if(drawOpts.addCRVInfo){
    std::vector<const CrvRecoPulseCollection*> crvpulse_list = std::get<1>(data.crvpulse_tuple);
    if(crvpulse_list.size() !=0) pass_data->AddCRVInfo(eveMng, firstLoop, data.crvpulse_tuple, eventScene, geomOpts.extracted, drawOpts.addCRVBars);
  }

  if(drawOpts.addCRVClusters){
    std::vector<const CrvCoincidenceClusterCollection*> crvcoin_list = std::get<1>(data.crvcoin_tuple);
    if(crvcoin_list.size() !=0) pass_data->AddCRVClusters(eveMng, firstLoop, data.crvcoin_tuple, eventScene, geomOpts.extracted, drawOpts.addCRVBars);
  }

  if(drawOpts.addCaloDigis){
    std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(data.calodigi_tuple);
    if(calodigi_list.size() !=0 ) pass_data->AddCaloDigis(eveMng, firstLoop, data.calodigi_tuple, eventScene, drawOpts.addCrystalDraw);
  }

  if(drawOpts.addClusters){
    //auto start1 = std::chrono::high_resolution_clock::now();
    std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(data.calocluster_tuple);
    if(calocluster_list.size() !=0 ) pass_data->AddCaloClusters(eveMng, firstLoopCalo, data.calocluster_tuple, eventScene, drawOpts.addCrystalDraw, t1, t2);
    //auto end1 = std::chrono::high_resolution_clock::now();
    // std::cout<<" time through fill caloclusters"<<std::chrono::duration<double, std::milli>(end1 - start1).count()<<" ms "<<std::endl;
  }

  std::vector<const HelixSeedCollection*> helix_list = std::get<1>(data.helix_tuple);
  if(drawOpts.addHelices and helix_list.size() !=0) {
    pass_data->AddHelixSeedCollection(eveMng, firstLoop, data.helix_tuple, eventScene);
  }


  if(drawOpts.addCosmicTracks){
    pass_data->AddCosmicTrackFit(eveMng, firstLoop, data.CosmicTrackSeedcol, eventScene);
  }
  if(drawOpts.addTimeClusters){
    std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(data.timecluster_tuple);
    if(timecluster_list.size() !=0) pass_data->AddTimeClusters(eveMng, firstLoop, data.timecluster_tuple, eventScene);
  }

  //... add MC:
  std::vector<const MCTrajectoryCollection*> mctrack_list = std::get<1>(data.mctrack_tuple);
  if(drawOpts.addMCTrajectories and mctrack_list.size() !=0){
    //auto start1 = std::chrono::high_resolution_clock::now();
    pass_mc->AddMCTrajectoryCollection(eveMng, firstLoop,  data.mctrack_tuple, eventScene, particleIds, geomOpts.extracted);
    //auto end1 = std::chrono::high_resolution_clock::now();
    //  std::cout<<" time through fill mctraj"<<std::chrono::duration<double, std::milli>(end1 - start1).count()<<" ms "<<std::endl;
  }

  // ... project these events onto 2D geometry:
  projectEvents(eveMng);
}

void REveMu2eMainWindow::makeGeometryScene(REX::REveManager *eveMng, GeomOptions geomOpt, std::string gdmlname)
{

  TGeoManager *geom = TGeoManager::Import(gdmlname.c_str());
  TGeoVolume* topvol = geom->GetTopVolume();
  gGeoManager->SetTopVolume(topvol);
  gGeoManager->SetTopVisible(kFALSE);
  TGeoNode* topnode = gGeoManager->GetTopNode();
  createProjectionStuff(eveMng);
  std::string name(topnode->GetName());
  {
    // make holders for different parts of geometry
    auto trackerholder = new REX::REveElement("Tracker");
    auto caloholder = new REX::REveElement("Calorimeter");
    auto crystalsholder = new REX::REveElement("CalorimeterCrystals");
    auto crvholder = new REX::REveElement("CRV");
    auto targetholder = new REX::REveElement("TargetElements");
    auto beamlineholder = new REX::REveElement("BeamlineElements");
    eveMng->GetGlobalScene()->AddElement(trackerholder);
    eveMng->GetGlobalScene()->AddElement(caloholder);
    eveMng->GetGlobalScene()->AddElement(crvholder);
    eveMng->GetGlobalScene()->AddElement(targetholder);
    eveMng->GetGlobalScene()->AddElement(beamlineholder);
    eveMng->GetGlobalScene()->AddElement(crystalsholder);

    REX::REveTrans trans;
    GeomDrawer(topnode, trans, beamlineholder, trackerholder, caloholder, crystalsholder, crvholder, targetholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt);
  }

  try {
    eveMng->Show();
  } catch (...) {
    std::cout<< "Blocking window" << std::endl;
  }

}

#ifdef __CINT__
ClassImp(REveMu2eMainWindow)
#endif
