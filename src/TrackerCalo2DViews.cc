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
    fCanvas = new TCanvas("Mu2eCanvas", "Mu2e stats", 800, 600);
    auto canvasViewer = evMng.SpawnNewViewer("Histogram viewer", "Mu2e Histogram");
    canvasViewer->AddScene(histScene);
}

void TrackerCalo2DViews::drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol) {
    if (!fCanvas || !fCanvasHolder) return;
    fCanvas->cd();
    fCanvas->Clear();
    
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    int planeId = 22;
    int panelId = 5;
    const mu2e::Plane& plane = tracker->getPlane(planeId);
    const mu2e::Panel& panel = plane.getPanel(panelId);
    
    TPad* pad = new TPad("p_0_0", "First Panel", 0.05, 0.05, 0.95, 0.95);
    pad->Draw();
    pad->cd();
    
    std::string frameTitle = "Plane" + std::to_string(planeId) + "Panel" + std::to_string(panelId) + "YZ view; Z [mm]; Y [mm]";
    TH2F* frame = new TH2F("frame", frameTitle.c_str(), 100, -20, 20, 100, -200, 200);
    frame->SetStats(0);
    frame->Draw();
    
    std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitDataMap;
    if(seedcol != nullptr){
        for (auto const& kseedptr : *seedcol){
            const mu2e::KalSeed& kseed = *kseedptr; 
            for (auto const& hit : kseed.hits()){
                hitDataMap[hit.strawId()] = &hit;
                std::cout << "Hit straw ID = " << hit.strawId() << std::endl;
            }
        }
    }
    
    double strawRadius = tracker->strawProperties()._strawOuterRadius;
    for (size_t iStraw=0; iStraw < panel.nStraws(); ++iStraw){
        const mu2e::Straw& straw = panel.getStraw(iStraw);
        CLHEP::Hep3Vector pos_l = panel.dsToPanel()*straw.getMidPoint();
        
        std::cout << "iStraw = " << iStraw << " pos_l x = " << pos_l.x() << " y = " << pos_l.y() << " z = " << pos_l.z() << " radius = " << strawRadius << " ID = " << straw.id() << std::endl;
        
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
            rcirc->SetFillStyle(0);
            rcirc->Draw();
            
            double sz = pos_l.z();
            double sy = pos_l.y();
            TGraph *g = new TGraph(1, &sz, &sy);
            g->SetMarkerColorAlpha(kWhite,0);
            g->SetName(Form("Straw %d: %.2f MeV",straw.id().getStraw(), hit->energyDep()));
            g->Draw("P SAME");
        }
    }
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.7, 0.92, Form("Plane %d : Panel %d", planeId, panelId));
    pad->Modified();
    pad->Update();
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
