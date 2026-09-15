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
#include "canvas/Persistency/Provenance/EventID.h"
#include <map>
#include <vector>

namespace REX = ROOT::Experimental;

namespace mu2e {

class TrackerCalo2DViews {
public:
    TrackerCalo2DViews();
    virtual ~TrackerCalo2DViews();

    void createHistogramView();
    void createStationView();
    void createCaloView();
    void drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol, const art::EventID& eventID, bool useAlignedTracker = false, bool shifttracker = false);
    void drawTrackerXYView(const mu2e::KalSeedPtrCollection* seedcol, int run = 0, int subRun = 0, int event = 0);
    void drawCalorimeterDisk(const CaloClusterCollection* clustercol = nullptr, const mu2e::KalSeedPtrCollection* seedcol = nullptr);

private:
    REX::REvePointSet* fXYCanvasHolder{nullptr};
    REX::REvePointSet* fCaloDisk0CanvasHolder{nullptr};
    REX::REvePointSet* fCaloDisk1CanvasHolder{nullptr};
    TCanvas* fXYCanvas{nullptr};
    TCanvas* fCaloCanvas{nullptr};
    TCanvas* fCaloCanvas1{nullptr};
    std::vector<TCanvas*> fStationCanvases;
    std::vector<REX::REvePointSet*> fStationCanvasHolders;
    int fEventCount{0};
};

} // namespace mu2e

#endif
