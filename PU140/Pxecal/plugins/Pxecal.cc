// -*- C++ -*-
//
// Package:    Pxecal
// Class:      Pxecal
// 
/**\class Pxecal Pxecal.h 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jose Cupertino Ruiz Vargas,32 4-C14,+41227674949,
//         Created:  Wed Aug  7 11:57:33 CEST 2013
// $Id$
//
//

#include "Pxecal.h"

Pxecal::Pxecal(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("t","t");
  t->Branch("Run",               &run);
  t->Branch("Event",             &event);
  t->Branch("BZfield",           &bz);
  t->Branch("BunchN",            &bunch_n);
  t->Branch("Pileup",            &pileup);
  t->Branch("BeamSpotX0",        &beamspot_x0);
  t->Branch("BeamSpotY0",        &beamspot_y0);
  t->Branch("BeamSpotZ0",        &beamspot_z0);
  t->Branch("BeamSpotX0Error",   &beamspot_x0Error);
  t->Branch("BeamSpotY0Error",   &beamspot_y0Error);
  t->Branch("BeamSpotZ0Error",   &beamspot_z0Error);
  t->Branch("BeamWidthX",        &beam_widthX);
  t->Branch("BeamWidthY",        &beam_widthY);
  t->Branch("BeamSigmaZ",        &beam_sigmaZ);
  t->Branch("BeamWidthXError",   &beam_widthXError);
  t->Branch("BeamWidthYError",   &beam_widthYError);
  t->Branch("BeamSigmaZError",   &beam_sigmaZError);
  t->Branch("TrkTruthN",         &trktruth_n);
  t->Branch("TrkTruthId",        &trktruth_id);
  t->Branch("TrkTruthCharge",    &trktruth_charge);
  t->Branch("TrkTruthEnergy",    &trktruth_energy);
  t->Branch("TrkTruthEt",        &trktruth_et);
  t->Branch("TrkTruthMass",      &trktruth_mass);
  t->Branch("TrkTruthMt",        &trktruth_mt);
  t->Branch("TrkTruthPt",        &trktruth_pt);
  t->Branch("TrkTruthPhi",       &trktruth_phi);
  t->Branch("TrkTruthEta",       &trktruth_eta);
  t->Branch("TrkTruthVx",        &trktruth_vx);
  t->Branch("TrkTruthVy",        &trktruth_vy);
  t->Branch("TrkTruthVz",        &trktruth_vz);
  t->Branch("TrkTruthG4trkN",    &trktruth_g4trk_n);
  t->Branch("TrkTruthTrkHitN",   &trktruth_trkhit_n);
  t->Branch("TrkTruthTrkLayerN", &trktruth_trklayer_n);
  t->Branch("SimTrkN",           &simtrk_n);
  t->Branch("SimTrkPt",          &simtrk_pt);
  t->Branch("SimTrkEta",         &simtrk_eta);
  t->Branch("SimTrkPhi",         &simtrk_phi);
  t->Branch("SimTrkId",          &simtrk_id);
  t->Branch("SimTrkType",        &simtrk_type);
  t->Branch("SimTrkVx",          &simtrk_vx);
  t->Branch("SimTrkVy",          &simtrk_vy);
  t->Branch("SimTrkVz",          &simtrk_vz);
  t->Branch("GenVx",             &sim_vx);
  t->Branch("GenVy",             &sim_vy);
  t->Branch("GenVz",             &sim_vz);
  t->Branch("GenPartN",          &genpart_n);
  t->Branch("GenPartE",          &genpart_e);
  t->Branch("GenPartEt",         &genpart_et);
  t->Branch("GenPartPt",         &genpart_pt);
  t->Branch("GenPartEta",        &genpart_eta);
  t->Branch("GenPartPhi",        &genpart_phi);
  t->Branch("GenPartCharge",     &genpart_charge);
  t->Branch("GenPartId",         &genpart_id);
  t->Branch("GenElN",            &genel_n);
  t->Branch("GenElE",            &genel_e);
  t->Branch("GenElEt",           &genel_et);
  t->Branch("GenElPt",           &genel_pt);
  t->Branch("GenElEta",          &genel_eta);
  t->Branch("GenElPhi",          &genel_phi);
  t->Branch("GenElCharge",       &genel_charge);
  t->Branch("GenElId",           &genel_id);
  t->Branch("EgN",               &egamma_n);
  t->Branch("EgE",               &egamma_e);
  t->Branch("EgEt",              &egamma_et);
  t->Branch("EgEta",             &egamma_eta);
  t->Branch("EgPhi",             &egamma_phi);
  t->Branch("EgGx",              &egamma_gx);
  t->Branch("EgGy",              &egamma_gy);
  t->Branch("EgGz",              &egamma_gz);
  t->Branch("EgCharge",          &egamma_charge);
  t->Branch("ClN",               &recHit_n);
  t->Branch("ClLx",              &recHit_lx);
  t->Branch("ClLy",              &recHit_ly);
  t->Branch("ClGx",              &recHit_gx);
  t->Branch("ClGy",              &recHit_gy);
  t->Branch("ClGz",              &recHit_gz);
  t->Branch("ClRho",             &recHit_rho);
  t->Branch("ClSubid",           &recHit_subid);
  t->Branch("ClLayer",           &recHit_layer);
  t->Branch("ClDisk",            &recHit_disk);
  t->Branch("ClSize",            &recHit_spread);
  t->Branch("ClSizeX",           &recHit_spreadx);
  t->Branch("ClSizeY",           &recHit_spready);
}

Pxecal::~Pxecal()
{
}

// ------------ method called for each event  ------------
void Pxecal::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  InitializeVectors();
  run = e.id().run();
  event = e.id().event();
  ///////////////////////////////////////////////////////////
  // Pileup  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<PileupSummaryInfo> > puinfo;
  e.getByLabel( "addPileupInfo", puinfo );
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  bunch_n=0;
  for(PVI = puinfo->begin(); PVI != puinfo->end(); ++PVI) {
     pileup.push_back(PVI->getPU_NumInteractions());
     bunch_n++;
  }
  ///////////////////////////////////////////////////////////
  // Gen Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< reco::GenParticleCollection > genParticles;
  e.getByLabel( "genParticles", genParticles );
  genpart_n = 0;
  genel_n = 0;
  for ( size_t i = 0; i < genParticles->size(); ++i ){
      const reco::GenParticle & genParticle = genParticles->at(i);
      genpart_e.push_back( genParticle.energy() );
      genpart_et.push_back( genParticle.et() );
      genpart_pt.push_back( genParticle.pt() );
      genpart_eta.push_back( genParticle.eta() );
      genpart_phi.push_back( genParticle.phi() );
      genpart_charge.push_back( genParticle.charge() );
      genpart_id.push_back( genParticle.pdgId() );
      genpart_n++;
      if ( abs(genParticles->at(i).pdgId()) == 11  ) {
         const reco::GenParticle & genElectron = genParticles->at(i);
         genel_e.push_back( genElectron.energy() );
         genel_et.push_back( genElectron.et() );
         genel_pt.push_back( genElectron.pt() );
         genel_eta.push_back( genElectron.eta() );
         genel_phi.push_back( genElectron.phi() );
         genel_charge.push_back( genElectron.charge() );
         genel_id.push_back( genElectron.pdgId() );
         genel_n++;
      }
  }
  ///////////////////////////////////////////////////////////
  // SimTracks & SimVertices
  //////////////////////////////////////////////////////////
  edm::Handle< edm::SimTrackContainer >   simTrackHandle;
  edm::Handle< edm::SimVertexContainer >  simVertex;
  e.getByLabel( "g4SimHits", simTrackHandle );
  e.getByLabel( "g4SimHits", simVertex );
  edm::SimTrackContainer::const_iterator iterSimTracks;
  simtrk_n = 0;
  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {
      simtrk_n++;
      simtrk_pt.push_back( iterSimTracks->momentum().pt() );   
      simtrk_eta.push_back( iterSimTracks->momentum().eta() );   
      simtrk_phi.push_back( iterSimTracks->momentum().phi() );   
      simtrk_id.push_back( iterSimTracks->trackId() );   
      simtrk_type.push_back( iterSimTracks->type() );   
      int index = iterSimTracks->vertIndex();
      simtrk_vx.push_back( simVertex->at(index).position().x() );   
      simtrk_vy.push_back( simVertex->at(index).position().y() );   
      simtrk_vz.push_back( simVertex->at(index).position().z() );  
      // index==0 gets the primary vertex;
      if ( index==0 ){
         sim_vx = simVertex->at(index).position().x();   
         sim_vy = simVertex->at(index).position().y();   
         sim_vz = simVertex->at(index).position().z();   
      }
  }
  ///////////////////////////////////////////////////////////
  // Tracking Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<TrackingParticle> > mergedPH;
  e.getByLabel( "mix","MergedTrackTruth", mergedPH );
  trktruth_n = 0;
  for ( std::vector<TrackingParticle>::const_iterator iTrack = mergedPH->begin(); iTrack != mergedPH->end(); ++iTrack ){
      trktruth_n++;
      trktruth_id.push_back( iTrack->pdgId() ); 
      trktruth_charge.push_back( iTrack->charge() ); 
      trktruth_energy.push_back( iTrack->energy() );
      trktruth_et.push_back( iTrack->et() );
      trktruth_mass.push_back( iTrack->mass() );
      trktruth_mt.push_back( iTrack->mt() );
      trktruth_pt.push_back( iTrack->pt() );
      trktruth_phi.push_back( iTrack->phi() );
      trktruth_eta.push_back( iTrack->eta() );
      trktruth_vx.push_back( iTrack->vx() );
      trktruth_vy.push_back( iTrack->vy() );
      trktruth_vz.push_back( iTrack->vz() );
      trktruth_g4trk_n.push_back( iTrack->g4Tracks().size() );
      trktruth_trkhit_n.push_back( iTrack->numberOfTrackerHits() );
      trktruth_trklayer_n.push_back( iTrack->numberOfTrackerLayers() );
  }
  ///////////////////////////////////////////////////////////
  // L1 Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > Egamma;
  e.getByLabel( "SLHCL1ExtraParticles","EGamma", Egamma );
  egamma_n = 0;
  for (l1extra::L1EmParticleCollection::const_iterator it = Egamma->begin(); it!=Egamma->end(); ++it){
      egamma_e.push_back(it->energy());
      egamma_et.push_back(it->et());
      egamma_eta.push_back(it->eta());
      egamma_phi.push_back(it->phi());
      egamma_charge.push_back(it->charge());
      GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
      egamma_gx.push_back(pos.x());
      egamma_gy.push_back(pos.y());
      egamma_gz.push_back(pos.z());
      egamma_n++;
  }
  ///////////////////////////////////////////////////////////
  // Beam Spot  
  //////////////////////////////////////////////////////////
  edm::Handle<reco::BeamSpot> thebeamSpot;
  e.getByLabel( "BeamSpotFromSim", "BeamSpot", thebeamSpot );
  beamspot_x0      = thebeamSpot->x0();
  beamspot_y0      = thebeamSpot->y0();
  beamspot_z0      = thebeamSpot->z0();
  beamspot_x0Error = thebeamSpot->x0Error();
  beamspot_y0Error = thebeamSpot->y0Error();
  beamspot_z0Error = thebeamSpot->z0Error();
  beam_widthX      = thebeamSpot->BeamWidthX();
  beam_widthY      = thebeamSpot->BeamWidthY();
  beam_sigmaZ      = thebeamSpot->sigmaZ();
  beam_widthXError = thebeamSpot->BeamWidthXError();
  beam_widthYError = thebeamSpot->BeamWidthYError();
  beam_sigmaZError = thebeamSpot->sigmaZ0Error();
  ///////////////////////////////////////////////////////////
  // Magnetic Field  
  //////////////////////////////////////////////////////////
  edm::ESHandle<MagneticField> theMagField;
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  bz = theMagField.product()->inTesla(GlobalPoint(0,0,0)).z();
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get(geom);
  const TrackerGeometry* theTracker = geom.product();
  //////////////////////////////////////////////////////////
  // Topology 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerTopology> topo;
  es.get<IdealGeometryRecord>().get(topo);
  //////////////////////////////////////////////////////////
  // RecHits
  //////////////////////////////////////////////////////////
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( "siPixelRecHits", recHitColl );
  SiPixelRecHitCollection::const_iterator recHitIdIterator    = (recHitColl.product())->begin();
  SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd = (recHitColl.product())->end();
  recHit_n = 0;
  for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++ ) {
      SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
      DetId detId = DetId(detset.detId()); // Get the Detid object
      const GeomDet* geomDet( theTracker->idToDet(detId) );
      int subid = detId.subdetId();
      SiPixelRecHitCollection::DetSet::const_iterator iterRecHit    = detset.begin();
      SiPixelRecHitCollection::DetSet::const_iterator iterRecHitEnd = detset.end();
      // Loop over rechits for this detid
      for ( ; iterRecHit != iterRecHitEnd; ++iterRecHit) {
          LocalPoint lp = iterRecHit->localPosition();
          GlobalPoint GP = geomDet->surface().toGlobal(lp);
          double rho = sqrt(GP.x()*GP.x() + GP.y()*GP.y());
          if ( rho < 20 ){ // rejects outer tracker
             recHit_lx.push_back(lp.x());
             recHit_ly.push_back(lp.y());
             recHit_gx.push_back(GP.x());
             recHit_gy.push_back(GP.y());
             recHit_gz.push_back(GP.z());
             recHit_rho.push_back(rho);
             recHit_subid.push_back(subid);
             if ( subid==PixelSubdetector::PixelBarrel ){
                recHit_layer.push_back(topo->pxbLayer(detId.rawId()));
                recHit_disk.push_back(-99);
             }
             if ( subid==PixelSubdetector::PixelEndcap ){
                recHit_layer.push_back(-99);
                recHit_disk.push_back(topo->pxfDisk(detId()));
             }
             // Cluster size
             SiPixelRecHit::ClusterRef const& Cluster =  iterRecHit->cluster();
             recHit_spread.push_back(Cluster->size());
             recHit_spreadx.push_back(Cluster->sizeX());
             recHit_spready.push_back(Cluster->sizeY());
             recHit_n++;
           }
      } // close loop over rechits for this detid
  } 
  // Fill tree
  t->Fill();

}// close Pxecal::analyze

// Taken from atanis' code
GlobalPoint Pxecal::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0;
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 )
    {
      double ecalZ = 314*fabs(eta)/eta;

      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {
      double rperp = 129.0;
      double zface =  sqrt( cos( theta ) * cos( theta ) /
                            ( 1 - cos( theta ) * cos( theta ) ) *
                            rperp * rperp ) * abs( eta ) / eta;
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}

void Pxecal::InitializeVectors()
{
                pileup.clear();
           trktruth_id.clear();
       trktruth_charge.clear(); 
       trktruth_energy.clear();
           trktruth_et.clear();
         trktruth_mass.clear();
           trktruth_mt.clear();
           trktruth_pt.clear();
          trktruth_phi.clear();
          trktruth_eta.clear();
           trktruth_vx.clear();
           trktruth_vy.clear();
           trktruth_vz.clear();
      trktruth_g4trk_n.clear();
     trktruth_trkhit_n.clear();
   trktruth_trklayer_n.clear();
             simtrk_pt.clear();  
            simtrk_eta.clear();
            simtrk_phi.clear();
             simtrk_id.clear();
           simtrk_type.clear();
             simtrk_vx.clear();
             simtrk_vy.clear();
             simtrk_vz.clear();
             genpart_e.clear();
            genpart_et.clear();
            genpart_pt.clear();
           genpart_eta.clear();
           genpart_phi.clear();
        genpart_charge.clear();
            genpart_id.clear();
               genel_e.clear();
              genel_et.clear();
              genel_pt.clear();
             genel_eta.clear();
             genel_phi.clear();
          genel_charge.clear();
              genel_id.clear();
              egamma_e.clear();
             egamma_et.clear();
            egamma_eta.clear();
            egamma_phi.clear();
             egamma_gx.clear();
             egamma_gy.clear();
             egamma_gz.clear();
         egamma_charge.clear();
             recHit_lx.clear();
             recHit_ly.clear();
             recHit_gx.clear();
             recHit_gy.clear();
             recHit_gz.clear();
            recHit_rho.clear();
          recHit_subid.clear();
          recHit_layer.clear();
           recHit_disk.clear();
         recHit_spread.clear();
        recHit_spreadx.clear();
        recHit_spready.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pxecal);
