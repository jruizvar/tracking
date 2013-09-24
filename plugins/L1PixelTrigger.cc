// -*- C++ -*-
//
// Package:    L1PixelTrigger
// Class:      L1PixelTrigger
// 
/**\class L1PixelTrigger L1PixelTrigger.h 

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

#include "L1PixelTrigger.h"

L1PixelTrigger::L1PixelTrigger(const edm::ParameterSet& iConfig)

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
  t->Branch("TrkN",              &trk_n);
  t->Branch("TrkStatus",         &trk_status);
  t->Branch("TrkId",             &trk_id);
  t->Branch("TrkCharge",         &trk_charge);
  t->Branch("TrkEnergy",         &trk_energy);
  t->Branch("TrkEt",             &trk_et);
  t->Branch("TrkMass",           &trk_mass);
  t->Branch("TrkMt",             &trk_mt);
  t->Branch("TrkPt",             &trk_pt);
  t->Branch("TrkPhi",            &trk_phi);
  t->Branch("TrkEta",            &trk_eta);
  t->Branch("TrkVx",             &trk_vx);
  t->Branch("TrkVy",             &trk_vy);
  t->Branch("TrkVz",             &trk_vz);
  t->Branch("TrkPerp",           &trk_perp);
  t->Branch("TrkG4trkN",         &trk_g4trk_n);
  t->Branch("TrkHitN",           &trk_hit_n);
  t->Branch("TrkLayerN",         &trk_layer_n);
  t->Branch("VtxN",              &vtx_n);
  t->Branch("VtxInVolume",       &vtx_involume);
  t->Branch("VtxG4vtxN",         &vtx_g4vtx_n);
  t->Branch("VtxSourceN",        &vtx_source_n);
  t->Branch("VtxDaughterN",      &vtx_daughter_n);
  t->Branch("VtxX",              &vtx_x);
  t->Branch("VtxY",              &vtx_y);
  t->Branch("VtxZ",              &vtx_z);
  t->Branch("VtxPerp",           &vtx_perp);
  t->Branch("VtxEta",            &vtx_eta);
  t->Branch("VtxPhi",            &vtx_phi);
  t->Branch("G4TrkN",            &g4trk_n);
  t->Branch("G4TrkPt",           &g4trk_pt);
  t->Branch("G4TrkEta",          &g4trk_eta);
  t->Branch("G4TrkPhi",          &g4trk_phi);
  t->Branch("G4TrkId",           &g4trk_id);
  t->Branch("G4TrkType",         &g4trk_type);
  t->Branch("G4VtxX",            &g4vtx_x);
  t->Branch("G4VtxY",            &g4vtx_y);
  t->Branch("G4VtxZ",            &g4vtx_z);
  t->Branch("G4VtxPerp",         &g4vtx_perp);
  t->Branch("G4VtxEta",          &g4vtx_eta);
  t->Branch("G4VtxPhi",          &g4vtx_phi);
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
  t->Branch("EgPerp",            &egamma_perp);
  t->Branch("EgCharge",          &egamma_charge);
  t->Branch("ClN",               &recHit_n);
  t->Branch("ClLx",              &recHit_lx);
  t->Branch("ClLy",              &recHit_ly);
  t->Branch("ClGx",              &recHit_gx);
  t->Branch("ClGy",              &recHit_gy);
  t->Branch("ClGz",              &recHit_gz);
  t->Branch("ClPerp",            &recHit_perp);
  t->Branch("ClEta",             &recHit_eta);
  t->Branch("ClPhi",             &recHit_phi);
  t->Branch("ClSubid",           &recHit_subid);
  t->Branch("ClLayer",           &recHit_layer);
  t->Branch("ClDisk",            &recHit_disk);
  t->Branch("ClSize",            &recHit_spread);
  t->Branch("ClSizeX",           &recHit_spreadx);
  t->Branch("ClSizeY",           &recHit_spready);
}

L1PixelTrigger::~L1PixelTrigger()
{
}

// ------------ method called for each event  ------------
void L1PixelTrigger::analyze(const edm::Event& e, const edm::EventSetup& es)
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
      if ( genParticles->at(i).pdgId() == 11  ) {
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
  g4trk_n = 0;
  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {
      g4trk_n++;
      g4trk_pt.push_back( iterSimTracks->momentum().pt() );   
      g4trk_eta.push_back( iterSimTracks->momentum().eta() );   
      g4trk_phi.push_back( iterSimTracks->momentum().phi() );   
      g4trk_id.push_back( iterSimTracks->trackId() );   
      g4trk_type.push_back( iterSimTracks->type() );   
      int index = iterSimTracks->vertIndex();
      // index==0 gets the primary vertex;
      if ( index==0 ){
         GlobalPoint pos(simVertex->at(index).position().x(),simVertex->at(index).position().y(),simVertex->at(index).position().z()); 
         g4vtx_x = pos.x();   
         g4vtx_y = pos.y();   
         g4vtx_z = pos.z();   
         g4vtx_perp = pos.perp();   
         g4vtx_eta = pos.eta();   
         g4vtx_phi = pos.phi();   
      }
  }
  ///////////////////////////////////////////////////////////
  // Tracking Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<TrackingParticle> > mergedPH;
  e.getByLabel( "mix","MergedTrackTruth", mergedPH );
  trk_n = 0;
  for ( std::vector<TrackingParticle>::const_iterator iTrack = mergedPH->begin(); iTrack != mergedPH->end(); ++iTrack ){
      trk_n++;
      trk_status.push_back( iTrack->longLived() );
      trk_id.push_back( iTrack->pdgId() ); 
      trk_charge.push_back( iTrack->charge() ); 
      trk_energy.push_back( iTrack->energy() );
      trk_et.push_back( iTrack->et() );
      trk_mass.push_back( iTrack->mass() );
      trk_mt.push_back( iTrack->mt() );
      trk_pt.push_back( iTrack->pt() );
      trk_phi.push_back( iTrack->phi() );
      trk_eta.push_back( iTrack->eta() );
      GlobalPoint pos(iTrack->vx(),iTrack->vy(),iTrack->vz()); 
      trk_vx.push_back( pos.x() );
      trk_vy.push_back( pos.y() );
      trk_vz.push_back( pos.z() );
      trk_perp.push_back( pos.perp() );
      trk_g4trk_n.push_back( iTrack->g4Tracks().size() );
      trk_hit_n.push_back( iTrack->numberOfTrackerHits() );
      trk_layer_n.push_back( iTrack->numberOfTrackerLayers() );
  }
  ///////////////////////////////////////////////////////////
  // Tracking Vertices  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<TrackingVertex> > mergedVH;
  e.getByLabel( "mix","MergedTrackTruth", mergedVH );
  vtx_n = 0;
  for ( std::vector<TrackingVertex>::const_iterator iVertex = mergedVH->begin(); iVertex != mergedVH->end(); ++iVertex ){
      vtx_n++;
      vtx_involume.push_back( iVertex->inVolume() );
      vtx_g4vtx_n.push_back( iVertex->nG4Vertices() );
      vtx_source_n.push_back( iVertex->nSourceTracks() );
      vtx_daughter_n.push_back( iVertex->nDaughterTracks() );
      GlobalPoint pos(iVertex->position().x(),iVertex->position().y(),iVertex->position().z()); 
      vtx_x.push_back( pos.x() );
      vtx_y.push_back( pos.y() );
      vtx_z.push_back( pos.z() );
      vtx_perp.push_back( pos.perp() );
      vtx_eta.push_back( pos.eta() );
      vtx_phi.push_back( pos.phi() );
  }
  ///////////////////////////////////////////////////////////
  // L1 Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > Egamma;
  e.getByLabel( "SLHCL1ExtraParticles","IsoEGamma", Egamma );
  egamma_n = 0;
  for (l1extra::L1EmParticleCollection::const_iterator it = Egamma->begin(); it!=Egamma->end(); ++it){
      egamma_e.push_back(it->energy());
      egamma_et.push_back(it->et());
      egamma_eta.push_back(it->eta());
      egamma_phi.push_back(it->phi());
      egamma_charge.push_back(it->charge());
      GlobalPoint pos = getCalorimeterPosition(it->phi(), it->eta(), it->energy());
      egamma_gx.push_back(pos.x());
      egamma_gy.push_back(pos.y());
      egamma_gz.push_back(pos.z());
      egamma_perp.push_back(pos.perp());
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
          double rho = GP.perp();
          if ( rho < 20 ){ // rejects outer tracker
             recHit_lx.push_back(lp.x());
             recHit_ly.push_back(lp.y());
             recHit_gx.push_back(GP.x());
             recHit_gy.push_back(GP.y());
             recHit_gz.push_back(GP.z());
             recHit_perp.push_back(rho);
             recHit_eta.push_back(GP.eta());
             recHit_phi.push_back(GP.phi());
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

}// close L1PixelTrigger::analyze

// Taken from atanis' code
GlobalPoint L1PixelTrigger::getCalorimeterPosition(double phi, double eta, double e) {
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

void L1PixelTrigger::InitializeVectors()
{
                pileup.clear();
            trk_status.clear();
                trk_id.clear();
            trk_charge.clear(); 
            trk_energy.clear();
                trk_et.clear();
              trk_mass.clear();
                trk_mt.clear();
                trk_pt.clear();
               trk_phi.clear();
               trk_eta.clear();
                trk_vx.clear();
                trk_vy.clear();
                trk_vz.clear();
              trk_perp.clear();
           trk_g4trk_n.clear();
             trk_hit_n.clear();
           trk_layer_n.clear();
          vtx_involume.clear();
           vtx_g4vtx_n.clear();
          vtx_source_n.clear();
        vtx_daughter_n.clear();
                 vtx_x.clear();
                 vtx_y.clear();
                 vtx_z.clear();
              vtx_perp.clear();
               vtx_eta.clear();
               vtx_phi.clear();
              g4trk_pt.clear();  
             g4trk_eta.clear();
             g4trk_phi.clear();
              g4trk_id.clear();
            g4trk_type.clear();
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
           egamma_perp.clear();
         egamma_charge.clear();
             recHit_lx.clear();
             recHit_ly.clear();
             recHit_gx.clear();
             recHit_gy.clear();
             recHit_gz.clear();
           recHit_perp.clear();
            recHit_eta.clear();
            recHit_phi.clear();
          recHit_subid.clear();
          recHit_layer.clear();
           recHit_disk.clear();
         recHit_spread.clear();
        recHit_spreadx.clear();
        recHit_spready.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelTrigger);
