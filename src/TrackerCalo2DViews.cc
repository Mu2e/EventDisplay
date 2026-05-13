#include "EventDisplay/inc/TrackerCalo2DViews.hh"
#include <TPad.h>
#include <TH2F.h>
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
#include <iostream>
#include <map>
#include <set>
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

template<class KTRAJ>
static void drawTrajectory2D(const KTRAJ& trajectory, const mu2e::Plane& plane, std::map<int, TPad*>& panelPads, const std::set<int>& activePanels)
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
        CLHEP::Hep3Vector global(pos.x(), pos.y(), pos.z());
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
        graph->SetLineColor(kRed);
        graph->SetLineWidth(3);
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

    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);
    graph->Draw("L SAME");
}

  void TrackerCalo2DViews::drawTrackerStation(const mu2e::KalSeedPtrCollection* seedcol){ //, const CaloDigiCollection* calodigicol) {
    // Note: We can draw without fCanvas/fCanvasHolder since you want offline viewing
     std::map<mu2e::StrawId, const mu2e::TrkStrawHitSeed*> hitDataMap;
     std::set<int> uniquePlanes;
     std::map<int, std::set<int>> activePanelsPerPlane;
    if (seedcol != nullptr) {
      for (auto const& kseedptr : *seedcol) {
        const mu2e::KalSeed& kseed = *kseedptr;
        for (auto const& hit : kseed.hits()) {
          mu2e::StrawId sid = hit.strawId();
          uniquePlanes.insert(sid.getPlane());
          hitDataMap[sid] = &hit;
	  activePanelsPerPlane[sid.getPlane()].insert(sid.getPanel());
        }
      }
    }
      mu2e::GeomHandle<mu2e::Tracker> tracker;
     
      double strawRadius = tracker->strawProperties()._strawOuterRadius;
      std::vector<int> planeIdA = {0,3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35};
      for (const auto& planeId : uniquePlanes) {
        std::cout << "Processing Canvas for Plane " << planeId << std::endl;
        std::array<int, 6> padMap;
        if(std::find(planeIdA.begin(), planeIdA.end(), planeId) != planeIdA.end())
          padMap = {2, 3, 6, 5, 4, 1};
        else
          padMap = {5, 4, 1, 2, 3, 6};
        TString canvasName = Form("Canvas_Plane_%d", planeId);
        TString canvasTitle = Form("Plane %d - V vs W View", planeId);
        auto cit = fPlaneCanvases.find(planeId);
        if (cit != fPlaneCanvases.end()) {
            delete cit->second;
            fPlaneCanvases.erase(cit);
        }
        TCanvas* planeCanvas = new TCanvas(canvasName, canvasTitle, 1000, 800);
        fPlaneCanvases[planeId] = planeCanvas;
        planeCanvas->Divide(2, 3, 0.005, 0.005);
        const mu2e::Plane& plane = tracker->getPlane(planeId);
	std::map<int, TPad*> panelPads;
        for (int panelId = 0; panelId < 6; ++panelId) {
          const mu2e::Panel& panel = plane.getPanel(panelId);
          // 1. Calculate Bounds for this specific panel
          /*double ymin = 1e9, ymax = -1e9, zmin = 1e9, zmax = -1e9;
          for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
             const mu2e::Straw& straw = panel.getStraw(iStraw);
             CLHEP::Hep3Vector pos_l = plane.dsToPlane()*straw.getMidPoint();
             ymin = std::min(ymin, pos_l.y());
             ymax = std::max(ymax, pos_l.y());
             zmin = std::min(zmin, pos_l.z());
             zmax = std::max(zmax, pos_l.z());
           }*/
           //std::cout<<"Plane = "<<planeId<<" panel = "<<panelId<<" ymin = "<<ymin<<" max = "<<ymax<<" zmin = "<<zmin<<" max = "<<zmax<<std::endl;
           // 2. Prepare the Pad
           planeCanvas->cd(padMap[panelId]);
	   panelPads[panelId] = dynamic_cast<TPad*>(gPad);
           gPad->SetBottomMargin(0.15);
           gPad->SetLeftMargin(0.15);
           gPad->SetFixedAspectRatio();
           TString frameName = Form("h_plane%d_panel%d", planeId, panelId);
           TString frameTitle = Form("Plane %d: Panel %d;W (mm);V (mm)", planeId, panelId);
           //double dz = zmax - zmin;
           //double dy = ymax - ymin;
           //double margin = std::max(dz, dy) * 0.1;
           //TH2F* frame = new TH2F(frameName, frameTitle, 100, zmin - margin, zmax + margin, 100, ymin - margin, ymax + margin);
           TH2F* frame = new TH2F(frameName, frameTitle, 100, -20, 20, 100, -170, 170);
           frame->SetStats(0);
           frame->Draw();
           // 3. Draw Straws and Hits
           for (size_t iStraw = 0; iStraw < panel.nStraws(); ++iStraw) {
             const mu2e::Straw& straw = panel.getStraw(iStraw);
             CLHEP::Hep3Vector pos_l = panel.dsToPanel()*straw.getMidPoint();
             // Base Straw Geometry
             TEllipse *circ = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);
             circ->SetLineColor(kGray + 1);
             circ->SetFillStyle(0);
             circ->Draw();
             // Check for Hits
             if (hitDataMap.count(straw.id())) {
               const auto* hit = hitDataMap[straw.id()];
               mu2e::WireHitState whs = hit->wireHitState();
               double rdrift = hit->driftRadius();
               TEllipse *rcirc = new TEllipse(pos_l.z(), pos_l.y(), rdrift, rdrift);
               TEllipse *hitcirc = new TEllipse(pos_l.z(), pos_l.y(), strawRadius, strawRadius);   
               if (!whs.active()) {
                 rcirc->SetFillStyle(0);
                 rcirc->SetLineColor(kBlack);
                 rcirc->SetLineStyle(2);
                 hitcirc->SetLineColor(kBlack);
                 hitcirc->SetLineWidth(0);
                 hitcirc->SetFillStyle(0);
               }
               else{
                 if (!whs.driftConstraint()) {
                   hitcirc->SetFillColor(kAzure - 9);
                   hitcirc->SetFillStyle(1001);
                   hitcirc->SetLineColor(kBlue);
                   rcirc->SetLineColor(kWhite);
                 } else {
                   rcirc->SetFillStyle(0);
                   rcirc->SetLineColor(kRed);
                   rcirc->SetLineWidth(hit->radialErr()*10);
                   hitcirc->SetLineColor(kBlack);
                   hitcirc->SetLineWidth(0);
                   hitcirc->SetFillStyle(0);
                 }
               }
               rcirc->Draw();
               hitcirc->Draw();
               // Tooltip/Data Graph
               double sz = pos_l.z();
               double sy = pos_l.y();
               TGraph *g = new TGraph(1, &sz, &sy);
               g->SetMarkerStyle(1);
               g->SetMarkerColorAlpha(kWhite, 0);
               g->SetName(Form("hit_%d_%d_%d", straw.id().getPlane(), straw.id().getPanel(), straw.id().getStraw()));
               g->Draw("P SAME");
             }
           }
        }
        
	if(seedcol != nullptr) {
          for (auto const& kseedptr : *seedcol) {
            const mu2e::KalSeed& kseed = *kseedptr;
            if(kseed.loopHelixFit()) {
              auto traj = kseed.loopHelixFitTrajectory();
              if (!traj) continue;
              drawTrajectory2D(*traj, plane, panelPads, activePanelsPerPlane[planeId]);
            } else if(kseed.centralHelixFit()) {
              auto traj = kseed.centralHelixFitTrajectory();
              if (!traj) continue;
              drawTrajectory2D(*traj, plane, panelPads, activePanelsPerPlane[planeId]);
            } else if(kseed.kinematicLineFit()) {
              auto traj = kseed.kinematicLineFitTrajectory();
              if (!traj) continue;
              drawTrajectory2D(*traj, plane, panelPads, activePanelsPerPlane[planeId]);
            }
          }
	}
        /*planeCanvas->cd();
        TLegend *leg = new TLegend(0.72, 0.82, 0.97, 0.96);
        leg->SetTextSize(0.022);
        leg->SetBorderSize(0);
        //leg->SetHeader("Hit key", "C");
        TEllipse *dummyInactive = new TEllipse();
        dummyInactive->SetFillStyle(0);
        dummyInactive->SetLineColor(kBlack);
        dummyInactive->SetLineStyle(2);
        TEllipse *dummyNoDrift = new TEllipse();
        dummyNoDrift->SetFillColor(kAzure - 9);
        dummyNoDrift->SetFillStyle(1001);
        dummyNoDrift->SetLineColor(kBlue);
        TEllipse *dummyDrift = new TEllipse();
        dummyDrift->SetFillColor(kAzure - 9);
        dummyDrift->SetFillStyle(1001);
        dummyDrift->SetLineColor(kRed);
        leg->AddEntry(dummyInactive, "Inactive (dashed circle)", "l");
        leg->AddEntry(dummyNoDrift,  "Active, no drift (full straw)", "f");
        leg->AddEntry(dummyDrift,    "Active, drift constraint", "f");
        leg->Draw();*/
        planeCanvas->Update();
      }
}

void TrackerCalo2DViews::drawTrackerXYView(const mu2e::KalSeedPtrCollection* seedcol) {
    mu2e::GeomHandle<mu2e::Tracker> tracker;

    double rInner = tracker->g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();
    double rOuter = tracker->g4Tracker()->getInnerTrackerEnvelopeParams().outerRadius();

    if (!fXYCanvas)
        fXYCanvas = new TCanvas("TrackerXY", "Tracker X-Y View", 800, 800);
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
        const CLHEP::Hep3Vector& mid = straw.getMidPoint();
        const CLHEP::Hep3Vector& dir = straw.getDirection();
        double hlen = straw.halfLength();

        // Full straw extent projected to XY, drawn in light grey
        TLine* strawLine = new TLine(
            mid.x() - hlen * dir.x(), mid.y() - hlen * dir.y(),
            mid.x() + hlen * dir.x(), mid.y() + hlen * dir.y());
        strawLine->SetLineColor(kGray);
        strawLine->SetLineWidth(1);
        strawLine->Draw();

        // Use the fitted POCA position (refPOCA_Upos) — bounded within the straw.
        // wireDist() is the raw TDC measurement and can exceed halfLength().
        CLHEP::Hep3Vector hitPos3D = straw.wirePosition(hit->refPOCA_Upos());
        double hx = hitPos3D.x();
        double hy = hitPos3D.y();
        float werr = hit->wireRes();

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
}

void TrackerCalo2DViews::drawCalorimeterDisk(const CaloClusterCollection* clustercol) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> calo;
    const mu2e::Disk& disk = calo->disk(0);

    if (!fCaloCanvas)
        fCaloCanvas = new TCanvas("calo_disk", "Calorimeter Disk 0", 1400, 1200);
    fCaloCanvas->cd();
    fCaloCanvas->Clear();
    fCaloCanvas->SetRightMargin(0.15);

    // Compute crystal bounds to set histogram range
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
    double margin = 20.0;
    xmin -= margin; xmax += margin;
    ymin -= margin; ymax += margin;

    // 2D histogram: X-Y are disk-local crystal coordinates, fill weight = energyDep (MeV)
    // 200x200 bins gives ~7 mm/bin, finer than a crystal (~33 mm) so COG fills the right cell
    TH2F* energyHist = new TH2F("calo_disk_energy",
                                 "Calorimeter Disk 0 - Cluster Energy Deposition;X (mm);Y (mm)",
                                 200, xmin, xmax, 200, ymin, ymax);
    energyHist->SetDirectory(0);
    energyHist->SetStats(0);
    gStyle->SetPalette(kBird);
    energyHist->GetZaxis()->SetTitle("Energy Deposition (MeV)");

    // Fill from disk-0 clusters using their COG projected to disk-local frame
    if (clustercol != nullptr) {
        for (const auto& cluster : *clustercol) {
            if (cluster.diskID() != 0) continue;
            CLHEP::Hep3Vector cog = cluster.cog3Vector();
            CLHEP::Hep3Vector localPos = calo->geomUtil().mu2eToDisk(0, cog);
            energyHist->Fill(localPos.x(), localPos.y(), cluster.energyDep());
        }
    }

    energyHist->Draw("COLZ");

    // Draw crystal outlines on top of the energy map
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

    fCaloCanvas->Modified();
    fCaloCanvas->Update();
}

  /*void TrackerCalo2DViews::redrawCanvas(const mu2e::KalSeedPtrCollection* seedcol) {
    if (!fCanvas || !fCanvasHolder) return;
    drawTrackerStation(seedcol);
    fCanvas->Modified();
    fCanvas->Update();
    TString json = TBufferJSON::ToJSON(fCanvas);
    fCanvasHolder->SetTitle(TBase64::Encode(json).Data());
    fCanvasHolder->StampObjProps();
  }*/

} // namespace mu2e
