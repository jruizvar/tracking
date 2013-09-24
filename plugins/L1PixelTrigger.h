#ifndef L1PixelTrigger_h
#define L1PixelTrigger_h

// system include files
#include <cmath>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// Tracker data formats
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// --- L1 Egamma dataFormats
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

// SiPixelRecHit dataFormat
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// --- GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// --- SimTracks & SimVertex
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

// --- Tracking Particles
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// --- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// ROOT includes
#include "TTree.h"

class L1PixelTrigger  : public edm::EDAnalyzer {
 public:

      explicit L1PixelTrigger(const edm::ParameterSet&);
      ~L1PixelTrigger();

 private:

      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      int                                 run;
      int                               event;
      float                                bz;
      int                             bunch_n;
      std::vector<int>                 pileup;
      float                       beamspot_x0;
      float                       beamspot_y0;
      float                       beamspot_z0;
      float                  beamspot_x0Error;
      float                  beamspot_y0Error;
      float                  beamspot_z0Error;
      float                       beam_widthX;
      float                       beam_widthY;
      float                       beam_sigmaZ;
      float                  beam_widthXError;
      float                  beam_widthYError;
      float                  beam_sigmaZError;
      int                               trk_n;
      std::vector<int>             trk_status;
      std::vector<int>                 trk_id;
      std::vector<int>             trk_charge;
      std::vector<float>           trk_energy;
      std::vector<float>               trk_et;
      std::vector<float>             trk_mass;
      std::vector<float>               trk_mt;
      std::vector<float>               trk_pt;
      std::vector<float>              trk_phi;
      std::vector<float>              trk_eta;
      std::vector<float>               trk_vx;
      std::vector<float>               trk_vy;
      std::vector<float>               trk_vz;
      std::vector<float>             trk_perp;
      std::vector<int>            trk_g4trk_n;
      std::vector<int>              trk_hit_n;
      std::vector<int>            trk_layer_n;
      int                               vtx_n;
      std::vector<int>           vtx_involume;
      std::vector<int>            vtx_g4vtx_n;
      std::vector<int>           vtx_source_n;
      std::vector<int>         vtx_daughter_n;
      std::vector<float>                vtx_x;
      std::vector<float>                vtx_y;
      std::vector<float>                vtx_z;
      std::vector<float>             vtx_perp;
      std::vector<float>              vtx_eta;
      std::vector<float>              vtx_phi;
      int                             g4trk_n;
      std::vector<float>             g4trk_pt;
      std::vector<float>            g4trk_eta;
      std::vector<float>            g4trk_phi;
      std::vector<int>               g4trk_id;
      std::vector<int>             g4trk_type;
      float                           g4vtx_x;
      float                           g4vtx_y;
      float                           g4vtx_z;
      float                        g4vtx_perp;
      float                         g4vtx_eta;
      float                         g4vtx_phi;
      int                           genpart_n;
      std::vector<float>            genpart_e;
      std::vector<float>           genpart_et;
      std::vector<float>           genpart_pt;
      std::vector<float>          genpart_eta;
      std::vector<float>          genpart_phi;
      std::vector<int>         genpart_charge;
      std::vector<int>             genpart_id;
      int                             genel_n;
      std::vector<float>              genel_e;
      std::vector<float>             genel_et;
      std::vector<float>             genel_pt;
      std::vector<float>            genel_eta;
      std::vector<float>            genel_phi;
      std::vector<int>           genel_charge;
      std::vector<int>               genel_id;
      int                            egamma_n;
      std::vector<float>             egamma_e;
      std::vector<float>            egamma_et;
      std::vector<float>           egamma_eta;
      std::vector<float>           egamma_phi;
      std::vector<float>            egamma_gx;
      std::vector<float>            egamma_gy;
      std::vector<float>            egamma_gz;
      std::vector<float>          egamma_perp;
      std::vector<int>          egamma_charge;
      int                            recHit_n;
      std::vector<float>            recHit_lx;
      std::vector<float>            recHit_ly;
      std::vector<float>            recHit_gx;
      std::vector<float>            recHit_gy;
      std::vector<float>            recHit_gz;
      std::vector<float>          recHit_perp;
      std::vector<float>           recHit_eta;
      std::vector<float>           recHit_phi;
      std::vector<int>           recHit_subid;
      std::vector<int>           recHit_layer;
      std::vector<int>            recHit_disk;
      std::vector<int>          recHit_spread;
      std::vector<int>         recHit_spreadx;
      std::vector<int>         recHit_spready;
      void InitializeVectors();
      // ----------Function taken from atani's code ---------------------------
      GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
      TTree* t;
};

#endif
