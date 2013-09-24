# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: DIGI -s DIGI,L1,DIGI2RAW --geometry ExtendedPhase2TkBE --conditions POSTLS261_V2::All --pileup AVE_140_BX_25ns --datamix NODATAMIXER --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE.customise,SLHCUpgradeSimulations/Configuration/customise_mixing.customise_NoCrossing --eventcontent FEVTDEBUG --datatier DIGI --filein root://osg-se.sprace.org.br//store/mc/Summer13/SingleElectronFlatPt5To50/GEN-SIM/UpgrdPhase2BE_POSTLS261_V2-v3/00000/225B081C-F6C4-E211-B833-00261894387A.root --pileup_input /store/mc/Summer13/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/UpgrdPhase2BE_POSTLS261_V2-v1/10000/000C99F4-F4BE-E211-9627-003048344B08.root --fileout file:SingleElectron_Pt5to50_PU140_DIGI.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('DIGI2RAW')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_E8TeV_AVE_16_BX_25ns_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBEReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('root://osg-se.sprace.org.br//store/mc/Summer13/Neutrino_Pt2to20_gun/GEN-SIM/UpgrdPhase2BE_POSTLS261_V2-v3/00000/3EEBE103-E8C4-E211-BD61-002618943976.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('DIGI nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:Neutrino_Pt2to20_PU140_DIGI.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI')
    )
)

# Additional output definition
from SimGeneral.MixingModule.trackingTruthProducerSelection_cfi import *
process.mix.digitizers.mergedtruth = trackingParticles
process.mix.digitizers.mergedtruth.volumeRadius = cms.double(200.0)
process.mix.digitizers.mergedtruth.volumeZ = cms.double(400.0)
process.mix.digitizers.mergedtruth.ignoreTracksOutsideVolume = cms.bool(True)
process.mix.digitizers.mergedtruth.select.chargedOnlyTP = cms.bool(False)
process.mix.digitizers.mergedtruth.select.signalOnlyTP = cms.bool(False) 
process.mix.digitizers.mergedtruth.select.stableOnlyTP = cms.bool(False)
process.mix.digitizers.mergedtruth.select.tipTP = cms.double(200.0)
process.mix.digitizers.mergedtruth.select.lipTP = cms.double(400.0)

# Other statements
from MinBiasSource import MinBias
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = MinBias 

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.FEVTDEBUGoutput_step)

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
