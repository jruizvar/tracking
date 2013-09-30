tracking
========

Full simulation with CMSSW_62X 

DIGI

`ssh lxplus5`
`cmsrel CMSSW_6_2_0_SLHC1`
`cd CMSSW_6_2_0_SLHC1/src`
`cmsenv`
`cmsRun DIGI_DIGI_L1_DIGI2RAW_PU.py`

NTUPLE

`ssh lxplus5`
`cmsrel CMSSW_6_2_0_SLHC1`
`cd CMSSW_6_2_0_SLHC1/src`
`cvs co -r jimb3June2013 SimDataFormats/SLHC`
`cvs co -r sh9Jul2013 SLHCUpgradeSimulations/L1CaloTrigger`
`cd SLHCUpgradeSimulations`
`mkedanlzr L1PixelTrigger`

(Replace the folder L1PixelTrigger/plugins) 

`cd CMSSW_6_2_0_SLHC1/src`
`scram b -j 8`
`cmsRun Neutrino_PU140_L1RECO.py`
