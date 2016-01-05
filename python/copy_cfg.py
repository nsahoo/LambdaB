import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring (
#"/store/mc/Summer12_DR53X/BuToKstarMuMu_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/006C6D81-8B4D-E311-BEBC-E0CB4E1A1194.root")                   
#"/store/mc/Summer11LegDR/LambdaBToLambdaMuMu_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/149B747C-D2CA-E311-B850-001E6739726F.root")
#"/store/mc/Summer12_DR53X/LambdaBToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/0000/003E1158-0AE3-E111-9A35-00A0D1EE8EE0.root")
"/store/mc/Summer11LegDR/LambdaBToLambdaJpsi_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/04B659AC-DFCA-E311-931B-001E67397F2B.root")
                             
)

# tell the process to only run over 500 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32 (100)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string ("LbToLzMuMu_MC_7TeV_aod_1000.root")                       
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)

