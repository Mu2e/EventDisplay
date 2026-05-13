#ifndef TrackerCalo2DViews_hh
#define TrackerCalo2DViews_hh

#include <TCanvas.h>
#include <TColor.h>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveScene.hxx>
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include <map>

namespace REX = ROOT::Experimental;

namespace mu2e {

class TrackerCalo2DViews {
public:
    TrackerCalo2DViews();
    virtual ~TrackerCalo2DViews();

    void createHistogramView();
  void drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol); //, const CaloDigiCollection* calodigicol);
  void drawTrackerXYView(const mu2e::KalSeedPtrCollection* seedcol);
  void drawCalorimeterDisk(const CaloClusterCollection* clustercol = nullptr);
  //void redrawCanvas(const mu2e::KalSeedPtrCollection* seedcol);

private:
    REX::REvePointSet* fCanvasHolder{nullptr};
    TCanvas* fCanvas{nullptr};
    TCanvas* fCaloCanvas{nullptr};
    TCanvas* fXYCanvas{nullptr};
    std::map<int, TCanvas*> fPlaneCanvases;
};

} // namespace mu2e

#endif
