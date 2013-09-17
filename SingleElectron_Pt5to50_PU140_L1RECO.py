# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: L1RECO -s L1TrackTrigger,RECO:pixeltrackerlocalreco --geometry ExtendedPhase2TkBE --conditions POSTLS261_V2::All --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE.customise,SLHCUpgradeSimulations/Configuration/customise_mixing.customise_NoCrossing --eventcontent FEVTDEBUG --datatier L1RECO --filein file:SingleElectron_Pt5to50_PU140_DIGI.root --fileout file:SingleElectron_Pt5to50_PU140_L1RECO.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBEReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:SingleElectron_Pt5to50_PU140_DIGI.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('L1RECO nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:dummy.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('L1RECO')
    )
)

# Additional output definition
#################################################################################################
# Calo trigger information
#################################################################################################
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.L1CaloTowerProducer.ECALDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.L1Calo = cms.Path(process.SLHCCaloTrigger)

#################################################################################################
# Beam Spot
#################################################################################################
process.BeamSpotFromSim=cms.EDProducer("BeamSpotFromSimProducer")
process.BeamSpot  = cms.Path(process.BeamSpotFromSim)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')

# Path and EndPath definitions
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.reconstruction_step = cms.Path(process.pixeltrackerlocalreco)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.demo = cms.EDAnalyzer('Pxecal')
process.p = cms.Path(process.demo)
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleElectron_PU140_ntuple.root') )

# Schedule definition
process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.reconstruction_step,process.endjob_step,process.L1Calo,process.BeamSpot,process.p)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.customise_mixing
from SLHCUpgradeSimulations.Configuration.customise_mixing import customise_NoCrossing 

#call to customisation function customise_NoCrossing imported from SLHCUpgradeSimulations.Configuration.customise_mixing
process = customise_NoCrossing(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE import customise 

#call to customisation function customise imported from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE
process = customise(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
