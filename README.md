# EcalClusteringML

Where:

    /afs/cern.ch/user/a/amassiro/work/ECAL/SIM/EcalClusteringML/
    
Starting point:

    https://github.com/bmarzocc/RecoSimStudies/tree/CMSSW_11_2_2_patch1

    
Install:

    scram project CMSSW_11_2_2_patch1
    cd CMSSW_11_2_2_patch1/src/
    cmsenv
    git cms-init
    git cms-merge-topic bmarzocc:CaloParticles_CMSSW_11_2_2_patch1 #if you produce RAW samples with caloParticles
    
    git clone https://github.com/bmarzocc/RecoSimStudies
    cd RecoSimStudies
    git checkout CMSSW_11_2_2_patch1
    cd -
    scram b -j 5
    
    
Produce GEN-SIM (QCD):

    cd RecoSimStudies/Dumpers/test/
    cmsRun QCD_Pt-15to7000_TuneCUETP8M1_Flat_14TeV-pythia8_cfi_GEN_SIM.py maxEvents=10 #QCD
    
Produce DIGI-RAW (Standard Mixing):

    cd RecoSimStudies/Dumpers/test/
    cmsRun step2_DIGI_L1_DIGI2RAW_HLT_PU_Run3_2021.py

Produce RECO:

    cd RecoSimStudies/Dumpers/test/
    cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py

Run general dumper (per crystal, CaloParticle, PFcluster, superCluster infos) on a RECO sample (produced in the previous steps):

    cd RecoSimStudies/Dumpers/
    cmsRun python/RecoSimDumper_cfg.py
    
    
New dumper
====
    
    mkdir ECALValidation
    cd ECALValidation/
    git clone git@github.com:amassiro/EcalClusteringML.git
    
    cmsenv
    scramv1 b -j 20
    
    
New dumper with only:

    - sim energy (and all sim hits?)
    - digis

for EE and EB
    
    
    
    
    
    
    
    