#ifndef TrackerCalo2DViews_hh
#define TrackerCalo2DViews_hh

#include <TCanvas.h>
#include <TColor.h>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveScene.hxx>
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

namespace REX = ROOT::Experimental;

namespace mu2e {

class TrackerCalo2DViews {
public:
    TrackerCalo2DViews();
    virtual ~TrackerCalo2DViews();

    void createHistogramView();
    void drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol);
    void redrawCanvas(const mu2e::KalSeedPtrCollection* seedcol);

private:
    REX::REvePointSet* fCanvasHolder{nullptr};
    TCanvas* fCanvas{nullptr};
};

} // namespace mu2e

#endif
