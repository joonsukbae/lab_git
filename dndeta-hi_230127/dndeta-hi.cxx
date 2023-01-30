// // Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "bestCollisionTable.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Index.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "TDatabasePDG.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis::hf_cuts_bplus_to_d0_pi;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiTracks = soa::Filtered<ExTracks>;
using Particles = soa::Filtered<aod::McParticles>;
using Particle = Particles::iterator;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using MyTracks = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::pidTPCPr>;

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kECend
};
enum {
  kTrigbegin = 0,
  kMBAND = 1,
  kTrigend
};
enum {
  kSpeciesbegin = 0,
  kK0short = 1,
  kLambda,
  kAntilambda,
  kSpeciesend
};

AxisSpec ZAxis = {60, -30, 30, "Z (cm)", "zaxis"};
AxisSpec DeltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
AxisSpec DCAAxis = {601, -3.01, 3.01, "", "DCA axis"};
AxisSpec EtaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
AxisSpec PhiAxis = {629, 0, 2 * M_PI, "Rad", "phi axis"};
AxisSpec PtAxis = {2401, -0.005, 24.005, "#it{p}_{T} (GeV/c)", "P_{T} axis"};
AxisSpec EvtClassAxis = {kECend - 1, kECbegin + 0.5, kECend - 0.5, "", "event class"};
AxisSpec TrigClassAxis = {kTrigend - 1, kTrigbegin + 0.5, kTrigend - 0.5, "", "trigger class"};
std::vector<double> centBinning = {0, 10., 20., 30., 40., 50., 60., 70., 80., 100};
AxisSpec CentAxis = {centBinning, "", "centrality"};
AxisSpec SpeciesAxis = {kSpeciesend - 1, kSpeciesbegin + 0.5, kSpeciesend - 0.5, "", "species class"};
AxisSpec massAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

struct MultiplicityCounter {
  Service<O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isMC{"isMC", false, "check if MC"};

  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  AxisSpec MultAxis = {multBinning, "N"};

  HistogramRegistry registry{
    "registry",
    {{"Events/Selection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}},                                                                          //
     {"hrecdndeta", "evntclass; triggerclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis, EtaAxis}}}, //
     {"hgendndeta", "evntclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis, EtaAxis}}},                              //
     {"hreczvtx", "evntclass; triggerclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis}}},                 //
     {"hgenzvtx", "evntclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis}}},                                              //
     {"PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}},                                                    //
     {"DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}},                                                                     //
     {"DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}},                                                                       //
     {"Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis}}},                                       //
     {"Centrality", " ; centrality_FT0C (%) ", {HistType::kTH1F, {CentAxis}}},                                                                            //
     {"Centrality_MBAND", " ; centrality_MBAND_FT0C (%) ", {HistType::kTH1F, {CentAxis}}},                                                                //
     {"hrecpt", " eventclass; centrality; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, PtAxis, PtAxis}}},                           //
     {"hgenpt", " eventclass; centrality; pt;  ", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, PtAxis}}},                                             //
     {"hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {SpeciesAxis, EvtClassAxis, massAxis}}}}};

  std::vector<int> usedTracksIds;

  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection()[kIsBBT0A] &
                        bc.selection()[kIsBBT0C]) != 0) {
        registry.fill(HIST("Events/Selection"), 5.);
        cols.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("Events/Selection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("Events/Selection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStat, "Collect event sample stats", false);

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                                    (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                                    true) &&
                                             ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA);
  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));
  std::vector<Double_t> tracketas;
  template <typename CollisionTypes, typename bcsTypes, typename ft0sTypes, typename fv0asTypes, typename tracksTypes, typename fullV0sTypes, typename atracksTypes>
  void runCounting(CollisionTypes const& collision, bcsTypes const& bcs, ft0sTypes const& ft0s, fv0asTypes const& fv0as, tracksTypes const& tracks, fullV0sTypes const& fullV0s, atracksTypes const& atracks, float cent)
  {
    const auto& foundBC = collision.template foundBC_as<BCsRun3>();
    float multT0A = 0;
    float multT0C = 0;
    float multV0A = 0;

    registry.fill(HIST("Centrality"), cent);
    registry.fill(HIST("Events/Selection"), 1.);
    auto z = collision.posZ();

    if (!useEvSel || collision.sel8()) {
      if (std::abs(z) < 10) {
        if (foundBC.has_ft0()) {
          for (auto amplitude : foundBC.ft0().amplitudeA()) {
            multT0A += amplitude;
          }
          for (auto amplitude : foundBC.ft0().amplitudeC()) {
            multT0C += amplitude;
          }
        } else {
          multT0A = multT0C = -999;
        }

        if (foundBC.has_fv0a()) {
          for (auto amplitude : foundBC.fv0a().amplitude()) {
            multV0A += amplitude;
          }
        } else {
          multV0A = -999;
        }
        registry.fill(HIST("Multiplicity"), multV0A, multT0A, multT0C);
        registry.fill(HIST("Centrality_MBAND"), cent);

        registry.fill(HIST("Events/Selection"), 2.);

        registry.fill(HIST("hreczvtx"), Double_t(kDATA), Double_t(kMBAND), cent, z);
        usedTracksIds.clear();

        tracketas.clear();
        for (auto& track : atracks) {
          auto otrack = track.template track_as<FiTracks>();
          tracketas.push_back(otrack.eta());
          registry.fill(HIST("PhiEta"), Double_t(kDATA), otrack.phi(), otrack.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), otrack.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), otrack.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, otrack.pt());
        }
        for (auto& track : tracks) {
          if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
            continue;
          }
          registry.fill(HIST("PhiEta"), Double_t(kDATA), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), track.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), track.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, track.pt());
          tracketas.push_back(track.eta());
        }

        for (auto eta : tracketas) {
          registry.fill(HIST("hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), cent, z, eta);
        }
        for (auto& v0 : fullV0s) {
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short());
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kLambda), v0.mLambda());
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), v0.mAntiLambda());
        }
      }
    }
  }
  void processCountingWithCent(
    MyCollisionsCent::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    aod::V0Datas const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    auto cent = collision.centFT0C();
    runCounting(collision, bcs, ft0s, fv0as, tracks, fullV0s, atracks, cent);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithCent, "Count tracks", false);

  void processCountingWithoutCent(
    MyCollisions::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    aod::V0Datas const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    auto cent = 50.;
    runCounting(collision, bcs, ft0s, fv0as, tracks, fullV0s, atracks, cent);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithoutCent, "Count tracks", false);

  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Preslice<FiTracks> perCol = aod::track::collisionId;
  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;

  void processMCCounting(
    soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&,
    aod::V0Datas const& fullV0s,
    Particles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    auto cent = 50.;

    for (auto& collision : collisions) {

      auto z = collision.posZ();
      if (useEvSel && !collision.sel8()) {
        if (std::abs(z) > 10) {
          continue;
        }
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      for (auto& v0 : fullV0s) {
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kK0short), v0.mK0Short());
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kLambda), v0.mLambda());
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kAntilambda), v0.mAntiLambda());
      }

      registry.fill(HIST("hreczvtx"), Double_t(kINEL), Double_t(kMBAND), cent, z);
      auto mcCollision = collision.mcCollision();
      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), cent, z, ttrack.mcParticle_as<Particles>().eta());
          registry.fill(HIST("hrecpt"), Double_t(kINEL), cent, ttrack.mcParticle_as<Particles>().pt(), ttrack.pt());

          registry.fill(HIST("PhiEta"), Double_t(kINEL), ttrack.phi(), ttrack.eta());
          registry.fill(HIST("DCAXY"), Double_t(kINEL), ttrack.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kINEL), ttrack.dcaZ());
        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (track.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), cent, z, track.mcParticle_as<Particles>().eta());
          registry.fill(HIST("hrecpt"), Double_t(kINEL), cent, track.mcParticle_as<Particles>().pt(), track.pt());

          registry.fill(HIST("PhiEta"), Double_t(kINEL), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), Double_t(kINEL), track.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kINEL), track.dcaZ());
        } else {
          // when secondary
        }
      }
      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          if (std::abs(p->Charge()) >= 3) {
            registry.fill(HIST("hgenpt"), Double_t(kINEL), cent, particle.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCCounting, "MC Count tracks", false);

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    auto cent = 50.;
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto genz = mcCollision.posZ();
    registry.fill(HIST("hgenzvtx"), Double_t(kINEL), cent, genz);
    for (auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          registry.fill(HIST("hgendndeta"), Double_t(kINEL), cent, genz, particle.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);
};

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

#include "Framework/runDataProcessing.h"

/// D0 analysis task
struct HfTaskD0 {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> recoFlag2Prong = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hYGenPrompt", "MC particles (matched, prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hPtProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hPtProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hDecLengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYSig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYSig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hPtProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hDecLengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtVsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng0", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng1", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthFinerBinning", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyFinerBinning", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0FinerBinning", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1FinerBinning", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0FinerBinning", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSFinerBinning", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtFinerBinning", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{500, -0., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAFinerBinning", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates)
  {
    for (auto& candidate : selectedD0Candidates) {
      if (!(candidate.hfflag() & 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
        continue;
      }

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), invMassD0ToPiK(candidate), candidate.pt());
        registry.fill(HIST("hMassFinerBinning"), invMassD0ToPiK(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), invMassD0barToKPi(candidate), candidate.pt());
        registry.fill(HIST("hMassFinerBinning"), invMassD0barToKPi(candidate), candidate.pt());
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      registry.fill(HIST("hNormalisedDecLength"), candidate.decayLengthNormalised(), candidate.pt());
      registry.fill(HIST("hNormalisedDecLengthxy"), candidate.decayLengthXYNormalised(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0ErrProng0"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hd0ErrProng1"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTS"), cosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCt"), ctD0(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("hDecLengthFinerBinning"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthxyFinerBinning"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0FinerBinning"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1FinerBinning"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0FinerBinning"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTSFinerBinning"), cosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCtFinerBinning"), ctD0(candidate), candidate.pt());
      registry.fill(HIST("hCPAFinerBinning"), candidate.cpa(), candidate.pt());
    }
  }

  void processMc(soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>& candidates,
                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    // MC rec.
    // Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : recoFlag2Prong) {
      if (!(candidate.hfflag() & 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>(), pdg::Code::kD0, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        auto ptRec = candidate.pt();
        auto yRec = yD0(candidate);
        if (candidate.isRecoHfFlag() >= selectionFlagHf) {
          registry.fill(HIST("hPtVsYRecSigRecoHFFlag"), ptRec, yRec);
        }
        if (candidate.isRecoTopol() >= selectionTopol) {
          registry.fill(HIST("hPtVsYRecSigRecoTopol"), ptRec, yRec);
        }
        if (candidate.isRecoCand() >= selectionCand) {
          registry.fill(HIST("hPtVsYRecSigRecoCand"), ptRec, yRec);
        }
        if (candidate.isRecoPid() >= selectionPid) {
          registry.fill(HIST("hPtVsYRecSig_RecoPID"), ptRec, yRec);
        }
        if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0) {
          registry.fill(HIST("hPtVsYRecSigReco"), ptRec, yRec); // rec. level pT
          registry.fill(HIST("hPtRecSig"), ptRec);
        }
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0) {
            registry.fill(HIST("hPtVsYRecSigPromptReco"), ptRec, yRec); // rec. level pT, prompt
            registry.fill(HIST("hPtRecSigPrompt"), ptRec);
          }
        } else {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0) {
            registry.fill(HIST("hPtVsYRecSigNonPromptReco"), ptRec, yRec); // rec. level pT, non-prompt
            registry.fill(HIST("hPtRecSigNonPrompt"), ptRec);
          }
        }
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
      }
      auto massD0 = invMassD0ToPiK(candidate);
      auto massD0bar = invMassD0barToKPi(candidate);
      auto ptCandidate = candidate.pt();
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto rapidityCandidate = yD0(candidate);
      auto declengthCandidate = candidate.decayLength();
      auto declengthxyCandidate = candidate.decayLengthXY();
      auto normaliseddeclengthCandidate = candidate.decayLengthNormalised();
      auto normaliseddeclengthxyCandidate = candidate.decayLengthXYNormalised();
      auto d0Prong0 = candidate.impactParameter0();
      auto d0Prong1 = candidate.impactParameter1();
      auto d0d0Candidate = candidate.impactParameterProduct();
      auto ctsCandidate = cosThetaStarD0(candidate);
      auto ctCandidate = ctD0(candidate);
      auto cpaCandidate = candidate.cpa();
      auto cpaxyCandidate = candidate.cpaXY();
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hPtProng0Sig"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Sig"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthSig"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYSig"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthSig"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYSig"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Sig"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Sig"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Sig"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSSig"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtSig"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPASig"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxySig"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hPtProng0Bkg"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Bkg"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthBkg"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYBkg"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthBkg"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYBkg"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Bkg"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Bkg"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Bkg"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSBkg"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtBkg"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPABkg"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxyBkg"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
          }
        }
      }
      if (candidate.isSelD0bar() >= selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == (1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          }
        }
      }
    }
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
        if (yCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > yCandMax) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("hPtGen"), ptGen);
        registry.fill(HIST("hPtVsYGen"), ptGen, yGen);
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hPtGenPrompt"), ptGen);
          registry.fill(HIST("hYGenPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenPrompt"), ptGen, yGen);
        } else {
          registry.fill(HIST("hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("hYGenNonPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenNonPrompt"), ptGen, yGen);
        }
        registry.fill(HIST("hEtaGen"), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(HfTaskD0, processMc, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc), adaptAnalysisTask<HfTaskD0>(cfgc)};
}
