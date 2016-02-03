################
#  variables
################
runMC               = True
run2012not2011      = False  # True for 2012 and False for 2011

print "\n@@@ CMSSW run configuration flags @@@"

if (runMC == True):
    print " running MC : ", runMC
else:
    print " not running MC, please xCheck !!"

if (run2012not2011 == True):
    print "---> 2012 MC running"
else:
    print "---> 2011 MC running"

print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"

##################
# cmssw configs
##################

import FWCore.ParameterSet.Config as cms
from lambdab_cfi import process 

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",

#         fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/LambdaBToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/0000/003E1158-0AE3-E111-9A35-00A0D1EE8EE0.root')
        fileNames = cms.untracked.vstring(
#'/store/mc/Summer11LegDR/LambdaBToLambdaMuMu_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/149B747C-D2CA-E311-B850-001E6739726F.root')
'/store/mc/Summer11LegDR/LambdaBToLambdaJpsi_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/04B659AC-DFCA-E311-931B-001E67397F2B.root',
'/store/mc/Summer11LegDR/LambdaBToLambdaJpsi_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/0C764327-AECA-E311-8EE4-001E67397F8F.root',
'/store/mc/Summer11LegDR/LambdaBToLambdaJpsi_EtaPtFilter_7TeV-pythia6-evtgen/AODSIM/PU_S13_START53_LV6-v1/00000/0CC61275-D7CA-E311-B177-002590200894.root')
#'/store/mc/Summer12DR53X/LambdaBToJPsiLambda_lambdaFilter_8TeV-pythia6-evtgen/GEN-SIM-RECO/PU_RD2_START53_V19F-v2/30000/001ABDFD-CA41-E511-A0FD-008CFA008C04.root')
                        )

if (run2012not2011 == True):
    process.GlobalTag.globaltag = cms.string('START53_V19F::All') # run dependent MC: START53_V19F , Jpsi X MC: START53_V7G
    print "\n=> Global Tag : START53_V19F::All\n"
else:
    process.GlobalTag.globaltag = cms.string('START53_LV6A1::All') # GT for 7TeV legacy samples reprocessed in cmssw_5_3_X
    print "\n=> Global Tag : START53_LV6A1::All\n"

# do trigger matching for muons
triggerProcessName = 'HLT'


if (run2012not2011 == True):
    process.cleanMuonTriggerMatchHLT  = cms.EDProducer(
            # match by DeltaR only (best match by DeltaR)
            'PATTriggerMatcherDRLessByR',
                src                   = cms.InputTag('cleanPatMuons'),
                # default producer label as defined in
                # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                matched               = cms.InputTag('patTrigger'),
                matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced*",0,0)'),
                maxDeltaR             = cms.double(0.1),
                # only one match per trigger object
                resolveAmbiguities    = cms.bool(True),
                # take best match found per reco object (by DeltaR here, see above)
                resolveByMatchQuality = cms.bool(False))
else:
    process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
                src = cms.InputTag('cleanPatMuons'),
                # default producer label as defined in                                                                                                               
                # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py 
                matched = cms.InputTag('patTrigger'),
                matchedCuts = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
                maxDeltaR = cms.double(0.1),
                # only one match per trigger object
                resolveAmbiguities = cms.bool(True),
                # take best match found per reco object (by DeltaR here, see above) 
                resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
                src = cms.InputTag('cleanPatMuons'),
                matched = cms.InputTag('patTrigger'),
                matchedCuts = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
                maxDeltaR = cms.double(0.1),
                resolveAmbiguities = cms.bool(True),
                resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
                src = cms.InputTag('cleanPatMuons'),
                matched = cms.InputTag('patTrigger'),
                matchedCuts = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
                maxDeltaR = cms.double(0.1),
                resolveAmbiguities = cms.bool(True),
                resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
                src = cms.InputTag('cleanPatMuons'),
                matched = cms.InputTag('patTrigger'),
                matchedCuts = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
                maxDeltaR = cms.double(0.1),
                resolveAmbiguities = cms.bool(True),
                resolveByMatchQuality = cms.bool(False))
    
    

from PhysicsTools.PatAlgos.tools.trigTools import *

if (run2012not2011 == True):
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')
    
    g_TriggerNames_LastFilterNames = [
        ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass')  #5E32, v8.1-v8.3
        ]
    
else:
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')
    
    g_TriggerNames_LastFilterNames = [
            ('HLT_Dimuon7_LowMass_Displaced',
             'hltDisplacedmumuFilterLowMass'), #1E33, 1.4E33

            ('HLT_DoubleMu4_LowMass_Displaced',
             'hltDisplacedmumuFilterLowMass'), #2E33

            ('HLT_DoubleMu4p5_LowMass_Displaced',
             'hltDisplacedmumuFilterDoubleMu4p5LowMass'), #3E33, 5E33

            ('HLT_DoubleMu5_LowMass_Displaced',
             'hltDisplacedmumuFilterDoubleMu5LowMass'), #3E33, 5E33
            ]
    


    
g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
