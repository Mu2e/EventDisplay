#include "EventDisplay/inc/TrackerCalo2DViews.hh"
#include <TPad.h>
#include <TH2F.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TBufferJSON.h>
#include <TBase64.h>
#include <iostream>
#include <map>
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"

namespace mu2e {

TrackerCalo2DViews::TrackerCalo2DViews() {}
TrackerCalo2DViews::~TrackerCalo2DViews() {}

void TrackerCalo2DViews::createHistogramView() {
    auto &evMng = *REX::gEve;
    auto histScene = evMng.SpawnNewScene("Histograms", "Histogram Scene");
    fCanvasHolder = new REX::REvePointSet("CanvasData");
    histScene->AddElement(fCanvasHolder);
    fCanvas = new TCanvas("Mu2eCanvas", "Tracker Plane Y v/s Z View", 800, 600);
    auto canvasViewer = evMng.SpawnNewViewer("Histogram viewer", "Mu2e Histogram");
    canvasViewer->AddScene(histScene);
}

void TrackerCalo2DViews::drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol) {
    if (!fCanvas || !fCanvasHolder) return;
    fCanvas->cd();
    fCanvas->Clear();

    fCanvas->Divide(2, 3, 0.005, 0.005);
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    int planeId = 22;
    
    std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitDataMap;
    if(seedcol != nullptr){
        for (auto const& kseedptr : *seedcol){
            const mu2e::KalSeed& kseed = *kseedptr; 
            for (auto const& hit : kseed.hits()){
                hitDataMap[hit.strawId()] = &hit;
            }
        }
    }
    
    const mu2e::Plane& plane = tracker->getPlane(planeId);
    double strawRadius = tracker->strawProperties()._strawOuterRadius;

    for (int panelId = 0; panelId < 6; ++panelId) {
        const mu2e::Panel& panel = plane.getPanel(panelId);
        std::cout<<"Plane ID = "<<planeId<<" panel = "<<panelId<<std::endl;
        fCanvas->cd(panelId+1);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        std::string frameTitle = "Plane " + std::to_string(planeId) + " Panel " + std::to_string(panelId);
        TH2F* frame = new TH2F(Form("h_%d", panelId), frameTitle.c_str(), 100, -20, 20, 100, -200, 200);
        frame->SetStats(0);
        frame->Draw();
        std::cout<<"nstraws = "<<panel.nStraws()<<std::endl;
        for (size_t iStraw=0; iStraw < panel.nStraws(); ++iStraw){
            const mu2e::Straw& straw = panel.getStraw(iStraw);
            CLHEP::Hep3Vector pos_l = panel.dsToPanel()*straw.getMidPoint();
            
            TEllipse *circ = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
            circ->SetLineColor(kGray+2);
            circ->SetFillStyle(0);
            circ->Draw();
            
            if(hitDataMap.count(straw.id())){
                const auto* hit = hitDataMap[straw.id()];
                TEllipse *hitcirc = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
                hitcirc->SetLineColor(kBlack);
                hitcirc->SetFillStyle(0);
                hitcirc->Draw();
                
                double rdrift = hit->driftRadius();
                TEllipse *rcirc = new TEllipse(pos_l.z(), pos_l.y(), rdrift, rdrift);
                rcirc->SetLineColor(kRed);
                
                mu2e::WireHitState whs(mu2e::WireHitState::State(hit->strawHitState()), mu2e::StrawHitUpdaters::algorithm(hit->strawHitAlgo()), hit->strawHitKkshFlag());
                if (whs.active() && whs.driftConstraint()) {
                    rcirc->SetFillColor(kLightCyan);
                    rcirc->SetFillStyle(1);
                } else {
                    rcirc->SetFillStyle(0);
                }
                rcirc->Draw();
                
                double sz = pos_l.z();
                double sy = pos_l.y();
                TGraph *g = new TGraph(1, &sz, &sy);
                g->SetMarkerColorAlpha(kWhite,0);
                g->SetName(Form("Straw %d: %.2f MeV",straw.id().getStraw(), hit->energyDep()));
                g->Draw("P SAME");
            }
        }
    }
}

void TrackerCalo2DViews::redrawCanvas(const mu2e::KalSeedPtrCollection* seedcol) {
    if (!fCanvas || !fCanvasHolder) return;
    drawTrackerStation(seedcol);
    fCanvas->Modified();
    fCanvas->Update();
    
    TString json = TBufferJSON::ToJSON(fCanvas);
    fCanvasHolder->SetTitle(TBase64::Encode(json).Data());
    fCanvasHolder->StampObjProps();
}

} // namespace mu2e
