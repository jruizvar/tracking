#ifndef Pxecal_h
#define Pxecal_h

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
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// --- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// ROOT includes
#include "TTree.h"

class Pxecal  : public edm::EDAnalyzer {
 public:

      explicit Pxecal(const edm::ParameterSet&);
      ~Pxecal();

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
      int                          trktruth_n;
      std::vector<int>            trktruth_id;
      std::vector<int>        trktruth_charge;
      std::vector<float>      trktruth_energy;
      std::vector<float>          trktruth_et;
      std::vector<float>        trktruth_mass;
      std::vector<float>          trktruth_mt;
      std::vector<float>          trktruth_pt;
      std::vector<float>         trktruth_phi;
      std::vector<float>         trktruth_eta;
      std::vector<float>          trktruth_vx;
      std::vector<float>          trktruth_vy;
      std::vector<float>          trktruth_vz;
      std::vector<int>       trktruth_g4trk_n;
      std::vector<int>      trktruth_trkhit_n;
      std::vector<int>    trktruth_trklayer_n;
      int                            simtrk_n;
      std::vector<float>            simtrk_pt;
      std::vector<float>           simtrk_eta;
      std::vector<float>           simtrk_phi;
      std::vector<int>              simtrk_id;
      std::vector<int>            simtrk_type;
      std::vector<float>            simtrk_vx;
      std::vector<float>            simtrk_vy;
      std::vector<float>            simtrk_vz;
      float                            sim_vx;
      float                            sim_vy;
      float                            sim_vz;
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
      std::vector<int>          egamma_charge;
      int                            recHit_n;
      std::vector<float>            recHit_lx;
      std::vector<float>            recHit_ly;
      std::vector<float>            recHit_gx;
      std::vector<float>            recHit_gy;
      std::vector<float>            recHit_gz;
      std::vector<float>           recHit_rho;
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
