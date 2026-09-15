#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "canvas/Persistency/Provenance/EventID.h"
#include "EventDisplay/inc/MainWindow.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include <ROOT/RWebWindow.hxx>

namespace REX = ROOT::Experimental;
using namespace std;
using namespace mu2e;

//positions of components for DS component offsets
double x_d0 = 0; double x_d1 = 0; double x_cal = 0;double x_trk = 0;double x_ds3 = 0;double x_st = 0; double x_ds2 = 0; double x_ipa =0;
double y_d0 = 0; double y_d1 = 0; double y_cal  = 0;double y_trk = 0;double y_ds3 = 0;double y_st = 0; double y_ds2 = 0; double y_ipa =0;
double z_d0 = 0; double z_d1 = 0; double z_cal = 0;double z_trk = 0;double z_ds3 = 0;double z_st = 0; double z_ds2 = 0; double z_ipa =0;


double nCrystals = 674;
double FrontTracker_gdmltag;

// Get geom information from Event Display configs:
std::string drawoptfilename("EventDisplay/config/drawutils.txt");
SimpleConfig drawconfigf(drawoptfilename);



void MainWindow::makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int val, bool crystal1, bool crystal2, std::string name, int color)
{
  bool isMother = name.find("ProductionTargetMother") != string::npos;
  auto gss = n->GetVolume()->GetShape();
  auto b1s = new REX::REveGeoShape(n->GetName());

  b1s->InitMainTrans();
  b1s->RefMainTrans().SetFrom(trans.Array());
  b1s->SetShape(gss);

  int transpar = drawconfigf.getInt("BLtrans");
  if(name.find("CRS") != string::npos) transpar = drawconfigf.getInt("Crvtrans");
  if(name.find("Tracker") != string::npos) transpar = drawconfigf.getInt("TRKtrans") ;
  if(name.find("CaloDisk") != string::npos) transpar = drawconfigf.getInt("CALtrans") ;
  if(name.find("CaloWrapper") != string::npos ) transpar = drawconfigf.getInt("CRYtrans") ;

  b1s->SetMainTransparency(transpar);
  b1s->SetMainColor(color);
  b1s->SetEditMainColor(true);
  b1s-> SetEditMainTransparency(true);
  if(!isMother) holder->AddElement(b1s);
}

int j = 0;
int fp =0;

void MainWindow::showNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal,  std::vector<double> shift, bool print, bool single, int color) {
  ++level;
  if (level > maxlevel){
    return;
  }
  std::string name(n->GetName());
  j++;
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
      REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
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
        REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
        t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
        t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
        t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
        t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1]  + shift[1]; t(3,4) = tv[2] + shift[2];
       
        if(name.find("TrackerPlaneEnvelope_00") != string::npos) {
          FrontTracker_gdmltag = j;
        }
        if(caloshift){
          fp++;
          double d = 0;
          if(fp < nCrystals) { d = z_d0- z_cal ; }//
          else { d = z_d1- z_cal ; }
          if(crystal) { t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1] + shift[1] ; t(3,4) = tv[2] + d + shift[2]; }
          else { t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1] + shift[1]; t(3,4) = tv[2]+ shift[2]; }
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

void MainWindow::getOffsets(TGeoNode* n,const std::string& str, REX::REveTrans& trans, int maxlevel, int level, std::vector<std::pair<std::string, std::vector<float>>> & offsets) {
  ++level;
  if (level > maxlevel){
    return;
  }
  std::string name(n->GetName());
  j++;
  REX::REveTrans ctrans;
  ctrans.SetFrom(trans.Array());
  TGeoMatrix     *gm = n->GetMatrix();
  const Double_t *rm = gm->GetRotationMatrix();
  const Double_t *tv = gm->GetTranslation();
  REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
  t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
  t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
  t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
  t(1,4) = tv[0]  ; t(2,4) = tv[1]  ; t(3,4) = tv[2] ;

  ctrans *= t;
  std::vector<float> pos;
  std::pair<std::string, std::vector<float>> offset;
  pos.push_back(tv[0]);
  pos.push_back(tv[1]);
  pos.push_back(tv[2]);
  offset = make_pair(name, pos);
  offsets.push_back(offset);
  n->ls();
  int ndau = n->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * pn = n->GetDaughter(i);
      REX::REveTrans ctrans;
      ctrans.SetFrom(trans.Array());
      {
        TGeoMatrix     *gm = n->GetMatrix();
        const Double_t *rm = gm->GetRotationMatrix();
        const Double_t *tv = gm->GetTranslation();
        REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
        t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
        t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
        t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
        t(1,4) = tv[0]  ; t(2,4) = tv[1]  ; t(3,4) = tv[2] ;

        ctrans *= t;
        std::vector<float> pos;
        std::pair<std::string, std::vector<float>> offset;
        pos.push_back(tv[0]);
        pos.push_back(tv[1]);
        pos.push_back(tv[2]);
        offset = make_pair(name, pos);
        offsets.push_back(offset);
        n->ls();

      }
    getOffsets( pn, str, trans,  maxlevel, level, offsets);
  }
}

void MainWindow::GeomDrawerSol(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets) {
  // Note: Solenoid component positions (PS, TS1, TS5, DS3, DS2) still obtained from offsets
  // as they require hierarchical transformations from the geometry tree.
  // CRV coordinates could use GeometryService if needed for display alignment.
  double x_ps = 0; double y_ps=0; double z_ps=0;
  double x_ts1 = 0; double y_ts1=0; double z_ts1=0;
  double x_ts5 = 0; double y_ts5=0; double z_ts5=0;
  double x_ptm = 0; double y_ptm=0; double z_ptm=0;
  double x_ds3 = 0; double y_ds3 = 0; double z_ds3 = 0;
  for(unsigned int i = 0; i < offsets.size(); i++){

     if(offsets[i].first.find("DS3Vacuum") != string::npos)      {
        x_ds3 = offsets[i].second[0];
        y_ds3 = offsets[i].second[1];
        z_ds3 = offsets[i].second[2];
      }
     if(offsets[i].first.find("DS2Vacuum") != string::npos)
     {
        x_ds2 = offsets[i].second[0];
        y_ds2 = offsets[i].second[1];
        z_ds2 = offsets[i].second[2];
      }
      if(offsets[i].first.find("TS1Vacuum") != string::npos)
      {
        x_ts1= offsets[i].second[0];
        y_ts1 = offsets[i].second[1];
        z_ts1 = offsets[i].second[2];
      }
      if(offsets[i].first.find("TS5Vacuum") != string::npos)
      {
        x_ts5= offsets[i].second[0];
        y_ts5 = offsets[i].second[1];
        z_ts5 = offsets[i].second[2];
      }
      if(offsets[i].first.find("PSVacuum") != string::npos)
      {
        x_ps = offsets[i].second[0];
        y_ps = offsets[i].second[1];
        z_ps = offsets[i].second[2];
      }
      if(offsets[i].first.find("ProductionTargetMother") != string::npos)
      {
        x_ptm = x_ps + offsets[i].second[0];
        y_ptm = y_ps + offsets[i].second[1];
        z_ptm = z_ps + offsets[i].second[2];
      }

    }
    std::vector<double> shift;
    shift.push_back(-x_ds3);
    shift.push_back(-y_ds3);
    shift.push_back(-1*z_ds3-z_trk);
    // set tracker to be central in the frame (shift to 0,0,0)
    if(geomOpt.showPS){

        // outer vessel
        static std::vector <std::string> substring_psvac_vessel {"PSVacum",
        "PSEnclosureFlange",
        "PSEnclosureWindow_1",
        "PSEnclosureWindow_2",
        "PSEnclosureWindowFrameInside1"
        ,"PSEnclosureWindowFrameOutside0",
        "PSEnclosureWindowFrameInside0",
        "PSEnclosureWindowFrameOutside1"
        ,"PSEnclosureWindow_3",
        "PSEnclosurePipeHole",
        "PSEnclosureWindowFrameInside2",
        "PSEnclosureWindowFrameOutside2"
        ,"PSVacVesselInner",
        "PSVacVesselOuter",
        "PSRing1",
        "PSRing2",
        "PSVacVesselEndPlateD",
        "PSVacVesselEndPlateU"
        ,"psVacuumVesselVacuum",
        "PSEnclosureEndPlate",
        "TS1UpstreamEndwall",
        "PSEnclosureShell"};
         for(auto& i: substring_psvac_vessel){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
         
        static std::vector <std::string> substring_ptm  {"ProductionTargetMother"};
        for(auto& i: substring_ptm){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
        /*
        "PbarAbsTS1InPeg2","PbarAbsTS1InPeg3","PbarAbsTS1InPeg1",
        "PSShieldShell1",
"PSShieldShell2",
"PSShieldShell3",
"PSShieldShell4",
"PSShieldEndRing1",
"PSShieldEndRing2",
"PSShieldEndRing3",
"PbarAbsTS1In",
"PbarAbsTS1InSup",
"PbarAbsTS1InFrameUp",
"PbarAbsTS1InFrameDown",
"PSVacuumPbarAbsTS1InSupTab1",
"PbarAbsTS1InSupTab2",
"PbarAbsTS1InSupTab3"
        */
        
        static std::vector <std::string> substring_pt  {
            "ProductionTargetStartingCoreSection",
            "ProductionTargetCoreSection",
            "ProductionTargetFinTopSection",
            "ProductionTargetFinTopStartingSection",
            "ProductionTargetFinStartingSection",
            "ProductionTargetNegativeEndRing",
            "ProductionTarget"};
        for(auto& i: substring_pt){
          std::vector<double> shift_pt;
          shift_pt.push_back(x_ptm-x_ds3);
          shift_pt.push_back(y_ptm-y_ds3 );
          shift_pt.push_back(z_ptm-z_ds3-z_trk);
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift_pt, false, true, drawconfigf.getInt("BLColor"));
        }
      }
      if(geomOpt.showTS){
        static std::vector <std::string> substring_ts     {"TS1Vacuum",
        "TS1CryoInsVac",
        "TS1InnerCryoShell",
        "TS1OuterCryoShell",
        "TS1DownstreamEndwall",
        "TS2Vacuum",
        "TS2CryoInsVac",
        "TS2InnerCryoShell",
        "TS2OuterCryoShell",
        "TS2CryoInsVac",
        "TS3Vacuum",
        "TS3CryoInsVac",
        "TS3InnerCryoShell",
        "TS3OuterCryoShell",
        "TSudInterconnectu",
        "TSudInterconnectd",
        "TS4Vacuum",
        "TS4CryoInsVac",
        "TS4InnerCryoShell",
        "TS4OuterCryoShell",
        "TS5Vacuum",
        "TS5CryoInsVac",
        "TS5InnerCryoShell",
        "TS5OuterCryoShell"};
        for(auto& i: substring_ts){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
        static std::vector <std::string> substring_ts1_internals  {  
          "Coll11",
          "Coll12",
          "Coll13",
          "PbarAbsTS1Out",
          "TS1ThermalShieldEndPlateMLI1",
          "TS1ThermalShieldEndPlateMLI2",
          "TS1ThermalShieldEndPlateMid",
          "TS1InnerThermalShieldMLI1",
          "TS1InnerThermalShieldMid",
          "TS1InnerThermalShieldMLI2",
          "TS1OuterThermalShieldMLI1",
          "TS1OuterThermalShieldMid",
          "TS1OuterThermalShieldMLI2",
          "TS1_Coil1",
          "TS1_Coil2",
          "TS1_Coil3",
          "TS1_Coil4",
          "TS1CA"};
        for(auto& i: substring_ts1_internals){
          std::vector<double> shiftts1;
          shiftts1.push_back(x_ts1-x_ds3);
          shiftts1.push_back(y_ts1-y_ds3);
          shiftts1.push_back(z_ts1-z_ds3-z_trk);
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shiftts1, false, true, drawconfigf.getInt("BLColor"));
        }
         /*static std::vector <std::string> substring_ts3_internals  {  "Coll31","Coll32","pBarAbsSupport","PbarAbsDisk","PbarAbsWedge"};
        for(auto& i: substring_ts3_internals){
          std::vector<double> shiftts3;
          shiftts3.push_back(x_ts3-x_ds3);
          shiftts3.push_back(y_ts3-y_ds3);
          shiftts3.push_back(z_ts3-z_ds3-z_trk);
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shiftts3, false, true, drawconfigf.getInt("BLColor"));
        }*/
       
        static std::vector <std::string> substring_ts5_internals  {  "Coll51"};
        for(auto& i: substring_ts5_internals){
          std::vector<double> shiftts5;
          shiftts5.push_back(x_ts5-x_ds3);
          shiftts5.push_back(y_ts5-y_ds3);
          shiftts5.push_back(z_ts5-z_ds3-z_trk);
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shiftts5, false, true, drawconfigf.getInt("BLColor"));
        }
        
        
      }
      if(geomOpt.showDS){
        static std::vector <std::string> substring_ds3  {"DS3Vacuum"};
        for(auto& i: substring_ds3){
     
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
        static std::vector <std::string> substring_ds2  {"DS2Vacuum"};
        for(auto& i: substring_ds2){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
      }
  }

void MainWindow::GeomDrawerNominal(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, REX::REveElement* trackerholder, REX::REveElement* caloholder, REX::REveElement* crystalsholder, REX::REveElement* crvholder, REX::REveElement* targetholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets){

    
    double x_hall = 0;double x_world = 0;
    double y_hall = 0;double y_world = 0;
    double z_hall = 0;double z_world = 0;
    double x_crv = 0 ; double y_crv = 0 ; double z_crv = 0;
    x_trk = 0;y_trk = 0;z_trk = 0;
    double x_fp0 = 0; double y_fp0 = 0; double z_fp0 = 0;
    double x_fp1 = 0; double y_fp1 = 0; double z_fp1 = 0;
    for(unsigned int i = 0; i < offsets.size(); i++){
      
      if(offsets[i].first.find("World") != string::npos){
        x_world = offsets[i].second[0];
        y_world = offsets[i].second[1];
        z_world = offsets[i].second[2];
      }
      if(offsets[i].first.find("HallAir") != string::npos)
      {
        x_hall = x_world + offsets[i].second[0];
        y_hall = y_world + offsets[i].second[1];
        z_hall = z_world + offsets[i].second[2];
      }
      if(offsets[i].first.find("DS3Vacuum") != string::npos)      {
        x_ds3 = x_hall + offsets[i].second[0];
        y_ds3 = y_hall + offsets[i].second[1];
        z_ds3 = z_hall + offsets[i].second[2];
        z_crv = -1*offsets[i].second[2];
        x_crv = -1*offsets[i].second[0];
        y_crv = -1*offsets[i].second[1];
      }
      if(offsets[i].first.find("TrackerMother") != string::npos)       {
        x_trk = x_ds3 + offsets[i].second[0];
        y_trk = y_ds3 + offsets[i].second[1];
        z_trk = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CalorimeterMother") != string::npos)       {
        x_cal = x_ds3 + offsets[i].second[0];
        y_cal = y_ds3 + offsets[i].second[1];
        z_cal = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_0") != string::npos && offsets[i].first.find("CaloDisk_1") == string::npos){
        x_d0 = x_cal + offsets[i].second[0];
        y_d0 = y_cal + offsets[i].second[1];
        z_d0 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_1") != string::npos){
        x_d1 = x_cal + offsets[i].second[0];
        y_d1 = y_cal + offsets[i].second[1];
        z_d1 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("DS2Vacuum") != string::npos){
        x_ds2 = x_hall + offsets[i].second[0];
        y_ds2 = y_hall + offsets[i].second[1];
        z_ds2 = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("StoppingTargetMother") != string::npos) {
        x_st = x_ds2 + offsets[i].second[0];
        y_st = y_ds2 + offsets[i].second[1];
        z_st = z_ds2 + offsets[i].second[2];
      }
      if(offsets[i].first.find("protonabs1") != string::npos){
        x_ipa = x_st+ offsets[i].second[0];
        y_ipa = y_st + offsets[i].second[1];
        z_ipa = z_st + offsets[i].second[2];
      }
       if(offsets[i].first.find("CaloFP") != string::npos){
        x_fp0 = x_d0+ offsets[i].second[0];
        y_fp0 = y_d0 + offsets[i].second[1];
        z_fp0 = z_d0 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloFP") != string::npos){
        x_fp1 = x_d1 + offsets[i].second[0];
        y_fp1 = y_d1 + offsets[i].second[1];
        z_fp1 = z_d1 + offsets[i].second[2];
      }
    }
    std::vector<double> shift;

    shift.push_back(0);
    shift.push_back(0);
    shift.push_back(0);
    // set tracker to be central in the frame (shift to 0,0,0)
    if(geomOpt.showTracker){
      static std::vector <std::string> substring_tracker  {"TrackerPlaneEnvelope"};
      for(auto& i: substring_tracker){
        shift.at(0) = 0;
        shift.at(1) = 0;
        shift.at(2) = 0;
        showNodesByName(node,i,kFALSE, 0, trans, trackerholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("TRKColor"));
      }
    }
    
    if(geomOpt.showST and !geomOpt.extracted){
      static std::vector <std::string> substrings_stoppingtarget  {"TargetFoil"};
      shift.at(0) = x_st - x_trk;
      shift.at(1) = y_st - y_trk;
      shift.at(2) = z_st - z_trk;
      for(auto& i: substrings_stoppingtarget){
        showNodesByName(node,i,kFALSE, 0, trans, targetholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
      }
    }
    if(geomOpt.showCalo){
        static std::vector <std::string> substrings_disk  {"CaloDisk"};
        for(auto& i: substrings_disk){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, true, false, drawconfigf.getInt("CALColor") );
        }
      if(geomOpt.showCaloCrystals){
        static std::vector <std::string> substrings_crystals  {"CaloWrapper"};
        for(auto& i: substrings_crystals){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, true, false, drawconfigf.getInt("CALColor"));
        }
      }
      bool showPipes = false;
      if(showPipes){
        static std::vector <std::string> substrings_pipes1  {"CaloFP_pipe","CaloFP_pipe1","CaloFP_pipe2","CaloFP_pipe3","CaloFP_pipe4","CaloFP_pipe5"};
          for(auto& i: substrings_pipes1){
            shift.at(0) = x_fp1 - x_trk;
            shift.at(1) = y_fp1 - y_trk;
            shift.at(2) = z_fp1 - z_trk;
            showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, false,false, shift, true, false, drawconfigf.getInt("CALColor"));
            
          }
          static std::vector <std::string> substrings_pipes0  {"CaloFP_pipe","CaloFP_pipe1","CaloFP_pipe2","CaloFP_pipe3","CaloFP_pipe4","CaloFP_pipe5"};
          for(auto& i: substrings_pipes0){
            shift.at(0) = x_fp0 - x_trk; 
            shift.at(1) = y_fp0 - y_trk;
            shift.at(2) = z_fp0 - z_trk - 20/2; //crystal len/2
            showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, false, true, drawconfigf.getInt("CALColor"));
          }
      }
    }
    if(geomOpt.showST and !geomOpt.extracted){
      static std::vector <std::string> substrings_pa  {"protonabs1","protonabs3","protonabs4"};
      shift.at(0) = x_ipa - x_trk;
      shift.at(1) = y_ipa - y_trk;
      shift.at(2) = z_ds2 - z_trk;
      for(auto& i: substrings_pa){
        showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
      }
  }
  if(geomOpt.showCrv and !geomOpt.extracted){
    static std::vector <std::string> substrings_crv {"CRSmother_CRV"};
    shift.at(0) = x_crv;
    shift.at(1) = y_crv;
    shift.at(2) = z_crv - z_trk;
    for(auto& i: substrings_crv){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, false, drawconfigf.getInt("CrvColor"));
    }
  }
}

void MainWindow::GeomDrawerExtracted(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, REX::REveElement* trackerholder, REX::REveElement* caloholder, REX::REveElement* crystalsholder, REX::REveElement* crvholder, REX::REveElement* targetholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets, DrawOptions drawOpts){
    // Load tracker geometry for CRV Z-shift calculation in extracted mode
    std::string trackerfilename("Offline/Mu2eG4/geom/tracker_v7.txt");
    SimpleConfig trackerconfig(trackerfilename);

    // Get CRV Z-shift for extracted geometry alignment
    double tracker_half_length_cm = 0.0; 
    if(drawOpts.shiftTracker) tracker_half_length_cm = -125.0;
    double x_world = 0;
    double y_world = 0;
    double z_world = 0;
    double x_hall = 0;
    double y_hall = 0;
    double z_hall = 0;
    double x_crvex = 0; double y_crvex = 0;  double z_crvex = 0; 
    double x_crvt1 = 0; double y_crvt1 = 0;  double z_crvt1 = 0;     
    double x_crvt2 = 0; double y_crvt2 = 0;  double z_crvt2 = 0;
    // Note: cannot use GeomService due to positionining of this in separate thread and initialization order of geometry. Instead, get positions from offsets and apply appropriate shifts to align for display.
    for(unsigned int i = 0; i < offsets.size(); i++){
      if(offsets[i].first.find("World") != string::npos){
        x_world = offsets[i].second[0];
        y_world = offsets[i].second[1];
        z_world = offsets[i].second[2];
      }
      if(offsets[i].first.find("HallAir") != string::npos)
      {
        x_hall = x_world + offsets[i].second[0];
        y_hall = y_world + offsets[i].second[1];
        z_hall = z_world + offsets[i].second[2];
      }
      if(offsets[i].first.find("CRSmother_CRV_EX") != string::npos)       {
          x_crvex = x_hall + offsets[i].second[0];
          y_crvex = y_hall + offsets[i].second[1];
          z_crvex = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("CRSmother_CRV_T1") != string::npos)       {
          x_crvt1 = x_hall + offsets[i].second[0];
          y_crvt1 = y_hall + offsets[i].second[1];
          z_crvt1 = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("CRSmother_CRV_T2") != string::npos)       {
          x_crvt2 = x_hall + offsets[i].second[0];
          y_crvt2 = y_hall + offsets[i].second[1];
          z_crvt2 = z_hall + offsets[i].second[2];
      }
    }
    
    // Get detector positions from offsets
    double x_ds3 = 0;
    for(unsigned int i = 0; i < offsets.size(); i++){
      if(offsets[i].first.find("garageFakeDS3Vacuum") != string::npos)      {
        x_ds3 = x_hall + offsets[i].second[0];
        y_ds3 = y_hall + offsets[i].second[1];
        z_ds3 = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("TrackerMother") != string::npos)       {
        x_trk = x_ds3 + offsets[i].second[0];
        y_trk = y_ds3 + offsets[i].second[1];
        z_trk = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CalorimeterMother") != string::npos)       {
        x_cal = x_ds3 + offsets[i].second[0];
        y_cal = y_ds3 + offsets[i].second[1];
        z_cal = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_0") != string::npos && offsets[i].first.find("CaloDisk_1") == string::npos){
        x_d0 = x_cal + offsets[i].second[0];
        y_d0 = y_cal + offsets[i].second[1];
        z_d0 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_1") != string::npos){
        x_d1 = x_cal + offsets[i].second[0];
        y_d1 = y_cal + offsets[i].second[1];
        z_d1 = z_cal + offsets[i].second[2] + 22.5;
      }
    }
    std::vector<double> shift;

    shift.push_back(0);
    shift.push_back(0);
    shift.push_back(0);
    // set tracker to be central in the frame (shift to 0,0,0)
    if(geomOpt.showTracker){
      static std::vector <std::string> substring_tracker  {"TrackerPlaneEnvelope"};
      for(auto& i: substring_tracker){
        shift.at(0) = 0;
        shift.at(1) = 0;
        shift.at(2) = tracker_half_length_cm;
        showNodesByName(node,i,kFALSE, 0, trans, trackerholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("TRKColor"));
      }
    }
    // everything else needs to be shifted such that its relative to the tracker center at 0,0,0
  
    if(geomOpt.showCalo){
        static std::vector <std::string> substrings_disk0  {"CaloDisk"};
        for(auto& i: substrings_disk0){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk; //+ 50.0;
          showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, true, false, drawconfigf.getInt("CALColor") );
        }
        static std::vector <std::string> substrings_disk1  {"CaloDisk_1"};
        for(auto& i: substrings_disk1){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk; //+ 75.0;
          showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, true, false, drawconfigf.getInt("CALColor") );
        }
      if(geomOpt.showCaloCrystals){
        static std::vector <std::string> substrings_crystals0  {"CaloWrapper"};
        for(auto& i: substrings_crystals0){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk; //+ 50.0;
          showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, true, false, drawconfigf.getInt("CALColor"));
        }
        static std::vector <std::string> substrings_crystals1  {"CaloWrapper_1"};
        for(auto& i: substrings_crystals1){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk; //+ 75.0;
          showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, true, false, drawconfigf.getInt("CALColor"));
        }
      }
    }
  if(geomOpt.showCrv and geomOpt.extracted){

    static std::vector <std::string> substrings_ex {"CRSmotherLayer_CRV_EX"};
    shift.at(0) = x_crvex - x_trk  ;
    shift.at(1) = y_crvex - y_trk;
    shift.at(2) = z_crvex - z_trk; // + tracker_half_length_cm;
    for(auto& i: substrings_ex){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CrvColor"));
    }

    static std::vector <std::string> substrings_t1  {"CRSmotherLayer_CRV_T1"};
    shift.at(0) = x_crvt1 - x_trk;
    shift.at(1) = y_crvt1 - y_trk;
    shift.at(2) = z_crvt1 - z_trk; // + tracker_half_length_cm;
    for(auto& i: substrings_t1){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CrvColor"));
    }

    static std::vector <std::string> substrings_t2  {"CRSmotherLayer_CRV_T2"};
    shift.at(0) = x_crvt2 - x_trk;
    shift.at(1) = y_crvt2 - y_trk;
    shift.at(2) = z_crvt2 - z_trk; // + tracker_half_length_cm;

    for(auto& i: substrings_t2){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CrvColor"));
    }
  }
}

void MainWindow::projectEvents(REX::REveManager */*eveMng*/)
{
}

void MainWindow::createProjectionStuff(REX::REveManager */*eveMng*/)
{
}


void MainWindow::showEvents(REX::REveManager *eveMng, REX::REveScene* &eventScene, bool firstLoop, DataCollections &data, DrawOptions drawOpts, std::vector<int> particleIds, bool strawdisplay, GeomOptions geomOpts, KinKalOptions KKOpts, int run, int subRun, int event){
  if (!eventScene){
    std::cerr<<"ERROR: eventScene is not an REveScene!"<<std::endl;
    return;
  }
  eventScene->BeginAcceptingChanges();
  if(!firstLoop){
    eventScene->DestroyElements();
  }
  //...addReco:
  // for start and end of track
  double t1 = 0.;
  double t2 = 1696.;
  std::vector<const KalSeedPtrCollection*> track_list = std::get<1>(data.track_tuple);
  if(drawOpts.addTracks and track_list.size() !=0) {
    pass_data->FillKinKalTrajectory(eveMng, eventScene, data.track_tuple, KKOpts.addKalInter,  KKOpts.addTrkStrawHits, KKOpts.addTrkCaloHits, t1, t2, drawOpts.shiftTracker);
    //redrawCanvas(seedcol);
     if(drawOpts.addTrackerHist and track_list.size() !=0) {
      if (!fTrackerCalo2DViews) {
        fTrackerCalo2DViews = new TrackerCalo2DViews();
        fTrackerCalo2DViews->createHistogramView();
        fTrackerCalo2DViews->createStationView();
      }
      auto const& track_list = std::get<1>(data.track_tuple);
      const mu2e::KalSeedPtrCollection* seedcol = track_list[0];
      fTrackerCalo2DViews->drawTrackerStation(seedcol, art::EventID(run, subRun, event), drawOpts.useAlignedTracker, drawOpts.shiftTracker);
      fTrackerCalo2DViews->drawTrackerXYView(seedcol, run, subRun, event);
     }
  }
   if(drawOpts.addCrvTrack) {
     pass_data->AddCRVKalIntersection(eveMng, eventScene, data.track_tuple, KKOpts.addKalInter,  KKOpts.addTrkStrawHits, KKOpts.addTrkCaloHits, t1, t2, data.crvcoin_tuple, geomOpts.extracted, drawOpts.addCrvBars);
   }
  if(drawOpts.addComboHits) {
    std::vector<const ComboHitCollection*> combohit_list = std::get<1>(data.combohit_tuple);
    if(combohit_list.size() !=0 ) pass_data->AddComboHits(eveMng, data.combohit_tuple, eventScene, strawdisplay, drawOpts.addTrkErrBar);
  }
  if(drawOpts.addBkgClusters) {
    std::vector<const BkgClusterCollection*> bkgcluster_list = std::get<1>(data.bkgcluster_tuple);
    if(bkgcluster_list.size() !=0 ) pass_data->AddBkgClusters(eveMng, data.bkgcluster_tuple, eventScene);
  }
    if(drawOpts.addCrvRecoPulse){
      // removed AddCrvInfo call
    }

  /*if(drawOpts.addCrvClusters){
    std::vector<const CrvCoincidenceClusterCollection*> crvcoin_list = std::get<1>(data.crvcoin_tuple);
    if(crvcoin_list.size() !=0) pass_data->AddCrvClusters(eveMng, data.crvcoin_tuple, eventScene, geomOpts.extracted, drawOpts.addCrvBars);
    }*/

  if(drawOpts.addCaloDigis){
    std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(data.calodigi_tuple);
    if(calodigi_list.size() !=0 ) {
      pass_data->AddCaloDigis(eveMng, data.calodigi_tuple, eventScene);
    }
  }
  if(drawOpts.addClusters){
    std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(data.calocluster_tuple);
    if(calocluster_list.size() !=0 ) 
      pass_data->AddCaloClusters(eveMng, data.calocluster_tuple, eventScene, drawOpts.addCrystalDraw);
    if(drawOpts.addCaloHist and calocluster_list.size() !=0) {
      if (!fTrackerCalo2DViews)
        fTrackerCalo2DViews = new TrackerCalo2DViews();
      const CaloClusterCollection* clustercol = (calocluster_list.size() > 1) ? calocluster_list[1] : calocluster_list[0];
      auto const& track_list_calo = std::get<1>(data.track_tuple);
      const mu2e::KalSeedPtrCollection* seedcol_calo = track_list_calo.size() > 0 ? track_list_calo[0] : nullptr;
      fTrackerCalo2DViews->drawCalorimeterDisk(clustercol, seedcol_calo);
   }
  }

  std::vector<const HelixSeedCollection*> helix_list = std::get<1>(data.helix_tuple);
  if(drawOpts.addHelices and helix_list.size() !=0) {
    pass_data->AddHelixSeedCollection(eveMng, data.helix_tuple, eventScene);
  }


  if(drawOpts.addCosmicTracks){
    pass_data->AddCosmicTrackFit(eveMng, data.CosmicTrackSeedcol, eventScene);
  }
  if(drawOpts.addTimeClusters){
    std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(data.timecluster_tuple);
    if(timecluster_list.size() !=0) pass_data->AddTimeClusters(eveMng, data.timecluster_tuple, data.combohit_tuple, eventScene);
  }

  //... add MC:
  std::vector<const MCTrajectoryCollection*> mctrack_list = std::get<1>(data.mctrack_tuple);
  if(drawOpts.addMCTrajectories and mctrack_list.size() !=0){
    pass_mc->AddMCTrajectoryCollection(eveMng, data.mctrack_tuple, eventScene, particleIds, geomOpts.extracted);
  }

  std::vector<const SurfaceStepCollection*> surfstep_list = std::get<1>(data.surfstep_tuple);
  if(drawOpts.addSurfaceSteps and surfstep_list.size() !=0){
    pass_mc->AddSurfaceStepCollection(eveMng, data.surfstep_tuple, eventScene, particleIds, geomOpts.extracted);
  }
  std::vector<const SimParticleCollection*> sim_list = std::get<1>(data.sim_tuple);
  if(drawOpts.addSimParts and sim_list.size() !=0){
    pass_mc->AddSimParticleCollection(eveMng, data.sim_tuple, eventScene, particleIds, geomOpts.extracted);
  }

  // ... project these events onto 2D geometry:
  projectEvents(eveMng);

  // Run/Event overlay — added to event scene so it's recreated each event.
  // Mode 1 = NDC screen-space; position (0.005, 0.99) = top-left, TextAlign(13) = left-top anchor.
  auto* overlayScene = eveMng->SpawnNewScene("EventOverlay", "Event Overlay Scene");
  auto* ann = new REX::REveText("EventAnnotation");
  ann->SetText("Run: " + std::to_string(run) + " / Event: " + std::to_string(event));
  ann->SetMode(0);
  //ann->SetTextColor(kBlack);
  ann->SetMainColor(1);
  ann->SetPosition(REX::REveVector(0.0,2500.0, 0.0));
  ann->SetFontSize(150);
  //ann->SetFont("LiberationSans-Bold");
  //ann->SetTextAlign(13);
  //ann->SetDrawFrame(false);
  //eventScene->AddElement(ann);
  ann->StampObjProps();
  overlayScene->AddElement(ann);
  eveMng->GetDefaultViewer()->AddScene(overlayScene);
  std::string title = "Run: " + std::to_string(run) + " Event: " + std::to_string(event);
  eveMng->GetWorld()->SetTitle(title.c_str());
  //eveMng->GetEventScene()->AddElement(ann); 
  eventScene->EndAcceptingChanges();
}

void MainWindow::makeGeometryScene(REX::REveManager *eveMng, GeomOptions geomOpt, std::string gdmlname, DrawOptions drawOpts)
{
  std::vector<std::pair<std::string, std::vector<float>>> offsets; // a pair of offsets relative to the world
  TGeoManager *geom = TGeoManager::Import(gdmlname.c_str());
  TGeoVolume* topvol = geom->GetTopVolume();
  gGeoManager->SetTopVolume(topvol);
  gGeoManager->SetTopVisible(kFALSE);
  TGeoNode* topnode = gGeoManager->GetTopNode();
  createProjectionStuff(eveMng);
  std::string name(topnode->GetName());
  {
    // make holders for different parts of geometry
    //auto holder = new REX::REveElement("World");
    auto trackerholder = new REX::REveElement("Tracker");
    auto caloholder = new REX::REveElement("Calorimeter");
    auto crystalsholder = new REX::REveElement("CalorimeterCrystals");
    auto crvholder = new REX::REveElement("Crv");
    auto targetholder = new REX::REveElement("TargetElements");
    auto beamlineholder = new REX::REveElement("DSBeamlineElements");
    auto solenoidholder = new REX::REveElement("SolenoidElements");
    eveMng->GetGlobalScene()->AddElement(trackerholder);
    eveMng->GetGlobalScene()->AddElement(caloholder);
    eveMng->GetGlobalScene()->AddElement(crvholder);
    eveMng->GetGlobalScene()->AddElement(targetholder);
    eveMng->GetGlobalScene()->AddElement(beamlineholder);
    eveMng->GetGlobalScene()->AddElement(crystalsholder);

    REX::REveTrans trans;
    static std::vector <std::string> substrings {"*"};
    for(auto& i: substrings){
      getOffsets(topnode, i, trans, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), offsets);
    }

    if(!geomOpt.extracted) GeomDrawerNominal(topnode, trans, beamlineholder, trackerholder, caloholder, crystalsholder, crvholder, targetholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets);
    if(geomOpt.extracted) GeomDrawerExtracted(topnode, trans, beamlineholder, trackerholder, caloholder, crystalsholder, crvholder, targetholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets, drawOpts);
    if(geomOpt.showPS or geomOpt.showTS or geomOpt.showDS){
      eveMng->GetGlobalScene()->AddElement(solenoidholder);
      GeomDrawerSol(topnode, trans, solenoidholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets);
      
    }
  }

  try {
    eveMng->Show();
  } catch (...) {
    std::cout<< "Blocking window" << std::endl;
  }

}

#ifdef __CINT__
ClassImp(MainWindow)
#endif
