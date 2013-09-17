tracking
========

Full simulation with CMSSW_62X 

DIGI

`ssh lxplus5`
`cmsrel CMSSW_6_2_0_SLHC1`
`cd CMSSW_6_2_0_SLHC1/src`
`source digi_step.sh`

NTUPLE

`ssh lxplus5`
`cmsrel CMSSW_6_2_0_SLHC1`
`cd CMSSW_6_2_0_SLHC1/src`
`cvs co -r jimb3June2013 SimDataFormats/SLHC`
`cvs co -r sh9Jul2013 SLHCUpgradeSimulations/L1CaloTrigger`
`scram b -j 8`
