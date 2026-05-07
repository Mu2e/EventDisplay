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
    // Note: We can draw without fCanvas/fCanvasHolder since you want offline viewing
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    
    std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitDataMap;
    std::set<int> uniquePlanes;

    if (seedcol != nullptr) {
        for (auto const& kseedptr : *seedcol) {
            const mu2e::KalSeed& kseed = *kseedptr; 
            for (auto const& hit : kseed.hits()) {
                mu2e::StrawId sid = hit.strawId();
                uniquePlanes.insert(sid.getPlane());
                hitDataMap[sid] = &hit;
            }
        }
    }

    double strawRadius = tracker->strawProperties()._strawOuterRadius;
    std::array<int, 6> padMap = {5, 4, 1, 2, 3, 6};

    for (const auto& planeId : uniquePlanes) {
        std::cout << "Processing Canvas for Plane " << planeId << std::endl;
        
        // Fix: Proper TCanvas name and title string formatting
        TString canvasName = Form("Canvas_Plane_%d", planeId);
        TString canvasTitle = Form("Mu2e Tracker Plane %d - Y vs Z View", planeId);
        
        TCanvas* planeCanvas = new TCanvas(canvasName, canvasTitle, 1000, 800);
        planeCanvas->Divide(2, 3, 0.005, 0.005);

        const mu2e::Plane& plane = tracker->getPlane(planeId);

        for (int panelId = 0; panelId < 6; ++panelId) {
            const mu2e::Panel& panel = plane.getPanel(panelId);
            
            // 1. Calculate Bounds for this specific panel
            double ymin = 1e9, ymax = -1e9, zmin = 1e9, zmax = -1e9;
            for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
                const mu2e::Straw& straw = panel.getStraw(iStraw);
                CLHEP::Hep3Vector pos_l = straw.getMidPoint();
                ymin = std::min(ymin, pos_l.y());
                ymax = std::max(ymax, pos_l.y());
                zmin = std::min(zmin, pos_l.z());
                zmax = std::max(zmax, pos_l.z());
            }

            // 2. Prepare the Pad
            planeCanvas->cd(padMap[panelId]);
            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.15);

            TString frameName = Form("h_plane%d_panel%d", planeId, panelId);
            TString frameTitle = Form("Plane %d: Panel %d;Z (mm);Y (mm)", planeId, panelId);
            
            TH2F* frame = new TH2F(frameName, frameTitle, 100, zmin - 20, zmax + 20, 100, ymin - 10, ymax + 10);
            frame->SetStats(0);
            frame->Draw();

            // 3. Draw Straws and Hits
            for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
                const mu2e::Straw& straw = panel.getStraw(iStraw);
                CLHEP::Hep3Vector pos_l = straw.getMidPoint();
                
                // Base Straw Geometry
                TEllipse *circ = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
                circ->SetLineColor(kGray + 1);
                circ->SetFillStyle(0);
                circ->Draw();
                
                // Check for Hits
                if (hitDataMap.count(straw.id())) {
                    const auto* hit = hitDataMap[straw.id()];
                    
                    // Hit Outline
                    TEllipse *hitcirc = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
                    hitcirc->SetLineColor(kBlack);
                    hitcirc->SetLineWidth(2);
                    hitcirc->SetFillStyle(0);
                    hitcirc->Draw();
                    
                    // Drift Circle
                    double rdrift = hit->driftRadius();
                    TEllipse *rcirc = new TEllipse(pos_l.z(), pos_l.y(), rdrift, rdrift);
                    
                    mu2e::WireHitState whs = hit->wireHitState();
                    if (whs.active() && whs.driftConstraint()) {
                        rcirc->SetFillColor(kAzure - 9);
                        rcirc->SetFillStyle(1001);
                        rcirc->SetLineColor(kBlue);
                    } else {
                        rcirc->SetFillStyle(0);
                        rcirc->SetLineColor(kRed);
                        rcirc->SetLineStyle(2); // Dashed for inactive/unconstrained
                    }
                    rcirc->Draw();

                    // Tooltip/Data Graph
                    double sz = pos_l.z();
                    double sy = pos_l.y();
                    TGraph *g = new TGraph(1, &sz, &sy);
                    g->SetMarkerStyle(1);
                    g->SetMarkerColorAlpha(kWhite, 0);
                    g->SetName(Form("Straw %d: %.2f MeV", straw.id().getStraw(), hit->energyDep()));
                    g->Draw("P SAME");
                }
            }
        }
        // Update the canvas for offline viewing
        planeCanvas->Update();
        // Optional: planeCanvas->SaveAs(Form("Plane_%d.png", planeId));
    }
}
  
void TrackerCalo2DViews::redrawCanvas(const mu2e::KalSeedPtrCollection* seedcol) {
  //if (!fCanvas || !fCanvasHolder) return;
    drawTrackerStation(seedcol);
    /* fCanvas->Modified();
    fCanvas->Update();
    TString json = TBufferJSON::ToJSON(fCanvas);
    fCanvasHolder->SetTitle(TBase64::Encode(json).Data());
    fCanvasHolder->StampObjProps();*/
}

} // namespace mu2e
