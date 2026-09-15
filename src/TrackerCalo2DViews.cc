#include "EventDisplay/inc/TrackerCalo2DViews.hh"
#include <TPad.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TBufferJSON.h>
#include <TBase64.h>
#include <algorithm>
#include <array>
#include <limits>
#include <iostream>
#include <map>
#include <set>
#include <TROOT.h>
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

namespace mu2e {

TrackerCalo2DViews::TrackerCalo2DViews() {}
TrackerCalo2DViews::~TrackerCalo2DViews() {}

void TrackerCalo2DViews::createHistogramView() {
    auto &evMng = *REX::gEve;
    auto histScene = evMng.SpawnNewScene("Histograms", "Histogram Scene");
    // "Calo Lego" is the element name expected by rootui5.eve7.view.Lego
    fXYCanvasHolder = new REX::REvePointSet("Calo Lego");
    histScene->AddElement(fXYCanvasHolder);
    // Viewer name must be exactly "Lego" so the browser routes to rootui5.eve7.view.Lego
    auto canvasViewer = evMng.SpawnNewViewer("Lego", "Tracker X-Y View");
    canvasViewer->AddScene(histScene);
}

void TrackerCalo2DViews::createStationView() {
    // Scenes/viewers/holders are created lazily in drawTrackerStation.
}

void TrackerCalo2DViews::createCaloView() {
    auto &evMng = *REX::gEve;

    auto disk0Scene = evMng.SpawnNewScene("CaloDisk0", "Calo Disk 0 Scene");
    fCaloDisk0CanvasHolder = new REX::REvePointSet("Calo Lego");
    disk0Scene->AddElement(fCaloDisk0CanvasHolder);
    auto disk0Viewer = evMng.SpawnNewViewer("Lego", "Calorimeter Disk 0");
    disk0Viewer->AddScene(disk0Scene);

    auto disk1Scene = evMng.SpawnNewScene("CaloDisk1", "Calo Disk 1 Scene");
    fCaloDisk1CanvasHolder = new REX::REvePointSet("Calo Lego");
    disk1Scene->AddElement(fCaloDisk1CanvasHolder);
    auto disk1Viewer = evMng.SpawnNewViewer("Lego", "Calorimeter Disk 1");
    disk1Viewer->AddScene(disk1Scene);
}

template<class KTRAJ>
static void drawTrajectory2D(const KTRAJ& trajectory, const mu2e::Plane& plane, std::map<int, TPad*>& panelPads, const std::set<int>& activePanels, bool shifttracker)
{
    double t1 = trajectory.range().begin();
    double t2 = trajectory.range().end();
    double step = 0.5;
    std::map<int, TGraph*> panelGraphs;
    for (int pid : activePanels) {
        panelGraphs[pid] = new TGraph();
    }
    for(double t = t1; t <= t2; t += step) {
        auto pos = trajectory.position3(t);
        double posz = pos.z();
        if(shifttracker)
          posz = pos.z()+1250.0;
        CLHEP::Hep3Vector global(pos.x(), pos.y(), posz);
        for (int pid : activePanels) {
            const mu2e::Panel& panel = plane.getPanel(pid);
            CLHEP::Hep3Vector local = panel.dsToPanel() * global;
            double w = local.z();
            double v = local.y();
            auto* g = panelGraphs[pid];
            g->SetPoint(g->GetN(), w, v);
        }
    }
    for (auto& [pid, graph] : panelGraphs) {
        TPad* pad = panelPads.count(pid) ? panelPads.at(pid) : nullptr;
        if (!pad) continue;
        pad->cd();
        graph->SetLineColor(kMagenta);
        graph->SetLineWidth(2);
        graph->Draw("L SAME");
    }
}

template<class KTRAJ>
static void drawTrajectoryXY(const KTRAJ& trajectory)
{
    double t1 = trajectory.range().begin();
    double t2 = trajectory.range().end();
    double step = 0.5;

    TGraph* graph = new TGraph();
    for (double t = t1; t <= t2; t += step) {
        auto pos = trajectory.position3(t);
        graph->SetPoint(graph->GetN(), pos.x(), pos.y());
    }

    graph->SetLineColor(kMagenta);
    graph->SetLineWidth(3);
    graph->Draw("L SAME");
}

  void TrackerCalo2DViews::drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol, const art::EventID& eventID, bool useAlignedTracker, bool shifttracker) {
    // Collect hit data and identify which (plane, panel) pairs have hits.
    std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitDataMap;
    std::set<std::pair<int,int>> seenPanels;
    std::vector<std::pair<int,int>> panelsWithHits; // (planeId, panelId)

    if (seedcol != nullptr) {
        for (auto const& kseedptr : *seedcol) {
            for (auto const& hit : kseedptr->hits()) {
                mu2e::StrawId sid = hit.strawId();
                hitDataMap[sid] = &hit;
                if (hit.wireHitState().active()) {
                  auto key = std::make_pair((int)sid.getPlane(), (int)sid.getPanel());
                  if (seenPanels.insert(key).second)
                    panelsWithHits.push_back(key);
                }
            }
        }
    }

    std::sort(panelsWithHits.begin(), panelsWithHits.end());

    if (panelsWithHits.empty()) return;

    ProditionsHandle<TrackerStatus> trackerStatus_h_;
    TrackerStatus const& trackerStatus = *trackerStatus_h_.getPtr(eventID);

    // Select tracker geometry source based on useAlignedTracker flag.
    // Aligned (proditions): uses survey-corrected straw positions for the specific run.
    // Nominal (GeomHandle): uses the design geometry, no run-by-run corrections.
    const mu2e::Tracker* tracker = nullptr;
    mu2e::GeomHandle<mu2e::Tracker> nominalTracker;
    ProditionsHandle<Tracker> alignedTracker_h_;
    if (useAlignedTracker) {
        tracker = alignedTracker_h_.getPtr(eventID).get();
    } else {
        tracker = &(*nominalTracker);
    }

    double strawRadius = tracker->strawProperties()._strawOuterRadius;

    const int kPanelsPerCanvas = 18;
    int nPanels   = (int)panelsWithHits.size();
    int nCanvases = (nPanels + kPanelsPerCanvas - 1) / kPanelsPerCanvas;

    // Create canvases and Eve holders on demand; reuse existing ones.
    while ((int)fStationCanvases.size() < nCanvases) {
        int idx = (int)fStationCanvases.size();
        bool wasBatch = gROOT->IsBatch();
        gROOT->SetBatch(kTRUE);
        auto* c = new TCanvas(Form("TrackerStation%d", idx),
                              Form("Tracker Station View %d", idx), 2400, 1200);
        c->SetBatch(kTRUE);
        gROOT->SetBatch(wasBatch);
        fStationCanvases.push_back(c);
    }
    while ((int)fStationCanvasHolders.size() < nCanvases) {
        int idx = (int)fStationCanvasHolders.size();
        auto& evMng = *REX::gEve;
        auto* holder = new REX::REvePointSet("Calo Lego");
        auto* scene  = evMng.SpawnNewScene(Form("TrackerStation%d", idx),
                                            Form("Tracker Station Scene %d", idx));
        scene->AddElement(holder);
        auto* viewer = evMng.SpawnNewViewer("Lego",
                                            Form("Tracker Station View %d", idx));
        viewer->AddScene(scene);
        fStationCanvasHolders.push_back(holder);
    }

    // Maps needed for trajectory drawing (grouped by plane).
    std::map<std::pair<int,int>, TPad*> panelPadMap;
    std::map<int, std::set<int>> activePanelsPerPlane;

    // Draw panels, 18 per canvas in a 6x3 grid.
    for (int ci = 0; ci < nCanvases; ++ci) {
        TCanvas* canvas = fStationCanvases[ci];
        canvas->SetName(Form("TrackerStation%d_%d_%d_%d", ci,
                             eventID.run(), eventID.subRun(), eventID.event()));
        canvas->cd();
        canvas->Clear();
        canvas->Divide(6, 3, 0.005, 0.005);

        int start = ci * kPanelsPerCanvas;
        int end   = std::min(start + kPanelsPerCanvas, nPanels);

        for (int i = start; i < end; ++i) {
            auto [planeId, panelId] = panelsWithHits[i];
            activePanelsPerPlane[planeId].insert(panelId);

            canvas->cd(i - start + 1);
            panelPadMap[{planeId, panelId}] = dynamic_cast<TPad*>(gPad);
            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.15);
            gPad->SetFixedAspectRatio();

            const mu2e::Plane& plane = tracker->getPlane(planeId);
            const mu2e::Panel& panel = plane.getPanel(panelId);

            double yMin = std::numeric_limits<double>::max();
            double yMax = std::numeric_limits<double>::lowest();
            for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
                const mu2e::Straw& straw = panel.getStraw(iStraw);
                if (!hitDataMap.count(straw.id())) continue;
                CLHEP::Hep3Vector pos_l = panel.dsToPanel() * straw.strawPosition(hitDataMap.at(straw.id())->refPOCA_Upos());
                yMin = std::min(yMin, pos_l.y());
                yMax = std::max(yMax, pos_l.y());
            }
            if (yMin > yMax) { yMin = -170.0; yMax = 170.0; }

            TH2F* frame = new TH2F(
                Form("h_pl%d_pa%d", planeId, panelId),
                Form("Plane %d Panel %d;W (mm);V (mm)", planeId, panelId),
                100, -20, 20, 100, yMin - 20.0, yMax + 20.0);
            frame->SetStats(0);
            frame->Draw();

            for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
                const mu2e::Straw& straw = panel.getStraw(iStraw);
                float rupos = hitDataMap.count(straw.id()) ? hitDataMap.at(straw.id())->refPOCA_Upos() : 0.0f;
                CLHEP::Hep3Vector pos_l = panel.dsToPanel() * straw.strawPosition(rupos);
                CLHEP::Hep3Vector pos_hit = panel.dsToPanel() * straw.wirePosition(rupos);

                TEllipse* circ = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
                if (trackerStatus.noSignal(straw.id()) || trackerStatus.suppress(straw.id())) {
                    circ->SetFillColor(kGray);
                    circ->SetFillStyle(1001);
                    circ->SetLineColor(kGray);
                } else {
                    circ->SetLineColor(kGray + 1);
                    circ->SetFillStyle(0);
                }
                circ->Draw();

                if (!hitDataMap.count(straw.id())) continue;

                const auto* hit    = hitDataMap[straw.id()];
                mu2e::WireHitState whs = hit->wireHitState();
                double rdrift      = hit->driftRadius();
                TEllipse* rcirc   = new TEllipse(pos_hit.z(), pos_hit.y(), rdrift,      rdrift);
                TEllipse* hitcirc = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);

                if (!whs.active()) {
                    rcirc->SetFillStyle(0);
                    rcirc->SetLineColor(kBlack);
                    rcirc->SetLineStyle(2);
                    hitcirc->SetLineColor(kBlack);
                    hitcirc->SetLineWidth(0);
                    hitcirc->SetFillStyle(0);
                } else if (!whs.driftConstraint()) {
                    hitcirc->SetFillColor(kAzure - 9);
                    hitcirc->SetFillStyle(1001);
                    hitcirc->SetLineColor(kBlue);
                    rcirc->SetLineColor(kWhite);
                } else {
                    rcirc->SetFillStyle(0);
                    rcirc->SetLineColor(kRed);
                    rcirc->SetLineWidth(hit->radialErr() * 10);
                    hitcirc->SetLineColor(kBlack);
                    hitcirc->SetLineWidth(0);
                    hitcirc->SetFillStyle(0);
                }
                rcirc->Draw();
                hitcirc->Draw();

                double sz = pos_l.z(), sy = pos_l.y();
                TGraph* g = new TGraph(1, &sz, &sy);
                g->SetMarkerStyle(1);
                g->SetMarkerColorAlpha(kWhite, 0);
                g->SetName(Form("hit_%d_%d_%d",
                                straw.id().getPlane(), straw.id().getPanel(), straw.id().getStraw()));
                g->Draw("P SAME");
            }
        }
    }

    // Draw trajectories, grouped by plane so dsToPanel transforms are correct.
    if (seedcol != nullptr) {
        for (auto const& kseedptr : *seedcol) {
            const mu2e::KalSeed& kseed = *kseedptr;
            for (auto& [planeId, panelSet] : activePanelsPerPlane) {
                const mu2e::Plane& plane = tracker->getPlane(planeId);
                std::map<int, TPad*> thisPlanesPads;
                for (int pid : panelSet)
                    thisPlanesPads[pid] = panelPadMap.at({planeId, pid});

                if (kseed.loopHelixFit()) {
                    auto traj = kseed.loopHelixFitTrajectory();
                    if (traj) drawTrajectory2D(*traj, plane, thisPlanesPads, panelSet, shifttracker);
                } else if (kseed.centralHelixFit()) {
                    auto traj = kseed.centralHelixFitTrajectory();
                    if (traj) drawTrajectory2D(*traj, plane, thisPlanesPads, panelSet, shifttracker);
                } else if (kseed.kinematicLineFit()) {
                    auto traj = kseed.kinematicLineFitTrajectory();
                    if (traj) drawTrajectory2D(*traj, plane, thisPlanesPads, panelSet, shifttracker);
                }
            }
        }
    }

    // Serialize each canvas to its Eve holder.
    for (int ci = 0; ci < nCanvases; ++ci) {
        TCanvas* canvas = fStationCanvases[ci];
        canvas->Modified();
        canvas->Update();
        REX::REvePointSet* holder = fStationCanvasHolders[ci];
        if (holder) {
            TString json = TBufferJSON::ToJSON(canvas);
            holder->SetTitle(TBase64::Encode(json).Data());
            holder->SetMainColor(kWhite);
            holder->StampObjProps();
        }
    }
}

void TrackerCalo2DViews::drawTrackerXYView(const mu2e::KalSeedPtrCollection* seedcol, int run, int subRun, int event) {

    mu2e::GeomHandle<mu2e::Tracker> tracker;

    double rInner = tracker->g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();
    double rOuter = tracker->g4Tracker()->getInnerTrackerEnvelopeParams().outerRadius();

    if (!fXYCanvas){
      bool wasBatch = gROOT->IsBatch();
      gROOT->SetBatch(kTRUE);
      fXYCanvas = new TCanvas("TrackerXY", "Tracker X-Y View", 800, 800);
      fXYCanvas->SetBatch(kTRUE);
      gROOT->SetBatch(wasBatch);
    }
    fXYCanvas->SetName(Form("TrackerXY_%d_%d_%d", run, subRun, event));
    fXYCanvas->cd();
    fXYCanvas->Clear();
    gPad->SetFixedAspectRatio();

    double rMax = rOuter * 1.1;
    TH2F* frame = new TH2F("TrackerXY_frame", "Tracker X-Y View;X (mm);Y (mm)",
                            100, -rMax, rMax, 100, -rMax, rMax);
    frame->SetStats(0);
    frame->Draw();

    TEllipse* innerCircle = new TEllipse(0.0, 0.0, rInner, rInner);
    innerCircle->SetLineColor(kBlue + 1);
    innerCircle->SetLineWidth(2);
    innerCircle->SetFillStyle(0);
    innerCircle->Draw();

    TEllipse* outerCircle = new TEllipse(0.0, 0.0, rOuter, rOuter);
    outerCircle->SetLineColor(kBlue + 1);
    outerCircle->SetLineWidth(2);
    outerCircle->SetFillStyle(0);
    outerCircle->Draw();

    TLatex* evLabel = new TLatex(0.5, 0.91, Form("%d:%d:%d", run, subRun, event));
    evLabel->SetNDC(true);
    evLabel->SetTextAlign(22);
    evLabel->SetTextSize(0.035);
    evLabel->SetTextColor(kBlack);
    evLabel->Draw();

    if (seedcol == nullptr) {
        fXYCanvas->Modified();
        fXYCanvas->Update();
        return;
    }

    // Collect unique hit straws across all seeds (last seed wins for duplicate straws)
    std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitMap;
    for (auto const& kseedptr : *seedcol) {
        for (auto const& hit : kseedptr->hits()) {
            hitMap[hit.strawId()] = &hit;
        }
    }

    for (auto const& [sid, hit] : hitMap) {
        const mu2e::Straw& straw = tracker->getStraw(sid);
        const CLHEP::Hep3Vector mid = straw.strawPosition(hit->refPOCA_Upos());
        const CLHEP::Hep3Vector& dir = straw.getDirection();
        double hlen = straw.halfLength();

        // Full straw extent projected to XY, drawn in light grey
        TLine* strawLine = new TLine(
            mid.x() - hlen * dir.x(), mid.y() - hlen * dir.y(),
            mid.x() + hlen * dir.x(), mid.y() + hlen * dir.y());
        strawLine->SetLineColor(kGray);
        strawLine->SetLineWidth(1);
        strawLine->Draw();

        mu2e::WireHitState whs(mu2e::WireHitState::State(hit->_ambig), mu2e::StrawHitUpdaters::algorithm(hit->_algo), hit->_kkshflag);
        bool active = whs.active();
        if(active) {
            CLHEP::Hep3Vector hitPos3D = straw.wirePosition(hit->_wdist);
            double hx = hitPos3D.x();
            double hy = hitPos3D.y();
            float werr = hit->_werr;
            // Longitudinal error bar: ±werr along wire direction, projected to XY
            TLine* errBar = new TLine(
                                      hx - werr * dir.x(), hy - werr * dir.y(),
                                      hx + werr * dir.x(), hy + werr * dir.y());
            errBar->SetLineColor(kBlack);
            errBar->SetLineWidth(2);
            errBar->Draw();
            // Hit marker as TGraph so ROOT shows the title in the status bar on hover
            double hx_d = hx, hy_d = hy;
            TGraph* hitPoint = new TGraph(1, &hx_d, &hy_d);
            hitPoint->SetMarkerStyle(20);
            hitPoint->SetMarkerSize(0.95);
            hitPoint->SetMarkerColor(kRed);
            hitPoint->SetName(Form("hit_%d_%d_%d", sid.getPlane(), sid.getPanel(), sid.getStraw()));
            hitPoint->SetTitle(Form("Plane %d  Panel %d  Straw %d  rdrift=%.3f mm",
                                sid.getPlane(), sid.getPanel(), sid.getStraw(),
                                hit->driftRadius()));
            hitPoint->Draw("P SAME");
        }
    }

    // Draw track trajectory in XY for each seed
    for (auto const& kseedptr : *seedcol) {
        const mu2e::KalSeed& kseed = *kseedptr;
        if (kseed.loopHelixFit()) {
            auto traj = kseed.loopHelixFitTrajectory();
            if (!traj) continue;
            drawTrajectoryXY(*traj);
        } else if (kseed.centralHelixFit()) {
            auto traj = kseed.centralHelixFitTrajectory();
            if (!traj) continue;
            drawTrajectoryXY(*traj);
        } else if (kseed.kinematicLineFit()) {
            auto traj = kseed.kinematicLineFitTrajectory();
            if (!traj) continue;
            drawTrajectoryXY(*traj);
        }
    }

    fXYCanvas->Modified();
    fXYCanvas->Update();

    if (fXYCanvasHolder) {
        TString json = TBufferJSON::ToJSON(fXYCanvas);
        fXYCanvasHolder->SetTitle(TBase64::Encode(json).Data());
        fXYCanvasHolder->SetMainColor(kWhite);
        fXYCanvasHolder->StampObjProps();
    }
}

void TrackerCalo2DViews::drawCalorimeterDisk(const CaloClusterCollection* clustercol, const mu2e::KalSeedPtrCollection* seedcol) {
    ++fEventCount;
    if(!fCaloDisk0CanvasHolder)
      createCaloView();
    mu2e::GeomHandle<mu2e::DiskCalorimeter> calo;

    // Collect per-hit data from cluster hit vectors
    struct HitInfo { int crystalID; float time; float eDep; double cx; double cy; double dx; double dy; int diskID; };
    std::vector<HitInfo> allHits;
    if (clustercol != nullptr) {
        for (const auto& cluster : *clustercol) {
            for (const auto& hitPtr : cluster.caloHitsPtrVector()) {
                const mu2e::CaloHit& hit = *hitPtr;
                const mu2e::Crystal& cr  = calo->crystal(hit.crystalID());
                allHits.push_back({hit.crystalID(), hit.time(), hit.energyDep(),
                                   cr.localPosition().x(), cr.localPosition().y(),
                                   cr.size().x() / 2.0, cr.size().y() / 2.0,
                                   (int)cr.diskID()});
            }
        }
    }

    auto drawDisk = [&](int diskID, TCanvas*& canvas, const char* canvasName,
                         const char* canvasTitle, REX::REvePointSet* holder) {
        const mu2e::Disk& disk = calo->disk(diskID);

        if (!canvas) {
            bool wasBatch = gROOT->IsBatch();
            gROOT->SetBatch(kTRUE);
            canvas = new TCanvas(canvasName, canvasTitle, 1400, 1200);
            canvas->SetBatch(kTRUE);
            gROOT->SetBatch(wasBatch);
        }
        canvas->SetName(canvasName);
        canvas->cd();
        canvas->Clear();
        canvas->SetRightMargin(0.15);

        double xmin =  1e9, xmax = -1e9;
        double ymin =  1e9, ymax = -1e9;
        for (size_t icr = 0; icr < disk.nCrystals(); ++icr) {
            const mu2e::Crystal& crystal = disk.crystal(icr);
            CLHEP::Hep3Vector pos  = crystal.localPosition();
            CLHEP::Hep3Vector size = crystal.size();
            double dx = size.x() / 2.0;
            double dy = size.y() / 2.0;
            xmin = std::min(xmin, pos.x() - dx);
            xmax = std::max(xmax, pos.x() + dx);
            ymin = std::min(ymin, pos.y() - dy);
            ymax = std::max(ymax, pos.y() + dy);
        }

        TH2Poly* energyHist = new TH2Poly(
            Form("calo_disk%d", diskID),
            Form("Disk %d;X (mm);Y (mm)", diskID),
            xmin, xmax, ymin, ymax);
        energyHist->SetDirectory(0);
        energyHist->SetStats(0);
        gStyle->SetPalette(kBird);
        energyHist->GetZaxis()->SetTitleOffset(1.5);
        energyHist->GetZaxis()->SetTitle("edep (MeV)");

        for (size_t icr = 0; icr < disk.nCrystals(); ++icr) {
            const mu2e::Crystal& crystal = disk.crystal(icr);
            CLHEP::Hep3Vector pos  = crystal.localPosition();
            CLHEP::Hep3Vector size = crystal.size();
            double x  = pos.x();
            double y  = pos.y();
            double dx = size.x() / 2.0;
            double dy = size.y() / 2.0;
            energyHist->AddBin(x - dx, y - dy, x + dx, y + dy);
        }

        for (const auto& h : allHits) {
            if (h.diskID != diskID) continue;
            energyHist->Fill(h.cx, h.cy, h.eDep);
        }

        energyHist->Draw("COLZ");

        // Always draw the full crystal geometry grid
        for (size_t icr = 0; icr < disk.nCrystals(); ++icr) {
            const mu2e::Crystal& crystal = disk.crystal(icr);
            CLHEP::Hep3Vector pos  = crystal.localPosition();
            CLHEP::Hep3Vector size = crystal.size();
            double x  = pos.x();
            double y  = pos.y();
            double dx = size.x() / 2.0;
            double dy = size.y() / 2.0;
            TBox* box = new TBox(x - dx, y - dy, x + dx, y + dy);
            box->SetFillStyle(0);
            box->SetLineColor(kGray + 1);
            box->SetLineWidth(1);
            box->Draw();
        }

        // Overlay hit markers for this event
        for (const auto& h : allHits) {
            if (h.diskID != diskID) continue;
            TGraph* g = new TGraph(1, &h.cx, &h.cy);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.5);
            g->SetMarkerColorAlpha(kWhite, 0);
            g->SetName(Form("Crystal %d  time=%.2f ns  eDep=%.2f MeV", h.crystalID, h.time, h.eDep));
            g->Draw("P SAME");
        }

        // Overlay track projection for this event
        if (seedcol != nullptr) {
            for (auto const& kseedptr : *seedcol) {
                const mu2e::KalSeed& kseed = *kseedptr;
                if (kseed.loopHelixFit()) {
                    auto traj = kseed.loopHelixFitTrajectory();
                    if (traj) drawTrajectoryXY(*traj);
                } else if (kseed.centralHelixFit()) {
                    auto traj = kseed.centralHelixFitTrajectory();
                    if (traj) drawTrajectoryXY(*traj);
                } else if (kseed.kinematicLineFit()) {
                    auto traj = kseed.kinematicLineFitTrajectory();
                    if (traj) drawTrajectoryXY(*traj);
                }
            }
        }

        canvas->Modified();
        canvas->Update();
        if (holder) {
            TString json = TBufferJSON::ToJSON(canvas);
            holder->SetTitle(TBase64::Encode(json).Data());
            holder->SetMainColor(kWhite);
            holder->StampObjProps();
        }
    };

    drawDisk(0, fCaloCanvas,  Form("calo_disk0_%d", fEventCount), "Disk 0", fCaloDisk0CanvasHolder);
    drawDisk(1, fCaloCanvas1, Form("calo_disk1_%d", fEventCount), "Disk 1", fCaloDisk1CanvasHolder);
}

} // namespace mu2e
