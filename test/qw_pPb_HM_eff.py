import FWCore.ParameterSet.Config as cms

process = cms.Process("QVector")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.GlobalTag.globaltag = 'GR_P_V43::All'

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

#fN = cms.untracked.vstring();
#for line in open('flist').read().splitlines():
#	fN.append('file:'+line);
#
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:../../../../pPb_HM_1000_1_Bgt.root")
#        fileNames = cms.untracked.vstring("file:pPb_HM.root")
#        fileNames = cms.untracked.vstring("file:pA_run210638_evt40050923.root")
#        fileNames = cms.untracked.vstring("file:pA_run210638_evt38874293.root")
)

#import FWCore.PythonUtilities.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#JSONfile = 'Cert_210498-211631_HI_PromptReco_Collisions13_JSON_v2.txt'
#myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
#process.source.lumisToProcess.extend(myLumis)
#
#
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltHM100 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM100.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
#       "HLT_PAPixelTracks_Multiplicity130_v*",
#       "HLT_PAPixelTracks_Multiplicity160_v*",
#       "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM100.andOr = cms.bool(True)
process.hltHM100.throw = cms.bool(False)

process.hltHM130 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM130.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
#       "HLT_PAPixelTracks_Multiplicity160_v*",
##      "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM130.andOr = cms.bool(True)
process.hltHM130.throw = cms.bool(False)


process.hltHM160 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM160.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
#       "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM160.andOr = cms.bool(True)
process.hltHM160.throw = cms.bool(False)




process.hltHM190 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM190.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
        "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM190.andOr = cms.bool(True)
process.hltHM190.throw = cms.bool(False)

process.hltHM220 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM220.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
        "HLT_PAPixelTracks_Multiplicity190_v*",
        "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM220.andOr = cms.bool(True)
process.hltHM220.throw = cms.bool(False)

process.QVector100 = cms.EDAnalyzer('QWQVector'
	, centrality = cms.InputTag("centralityBin", "HFtowers")
	, trackTag = cms.untracked.InputTag('generalTracks')
	, vertexSrc = cms.untracked.InputTag('offlinePrimaryVertices', "")
	, pterrorpt = cms.untracked.double(0.1)
	, dzdzerror = cms.untracked.double(3.0)
	, d0d0error = cms.untracked.double(3.0)
	, minvz = cms.untracked.double(-1.0)
	, maxvz = cms.untracked.double(15.0)
	, minEta = cms.untracked.double(-2.4)
	, maxEta = cms.untracked.double(2.4)
	, minPt = cms.untracked.double(0.3)
	, maxPt = cms.untracked.double(3.0)
	, minCent = cms.untracked.int32(120)
	, maxCent = cms.untracked.int32(150)
	, epSrc = cms.untracked.InputTag("hiEvtPlane")
	, fweight = cms.untracked.InputTag('TrackCorrections_HIJING_538_OFFICIAL_Mar24.root')
	, bEff = cms.untracked.bool(True)
	, algoParameters = cms.vint32()
	, bGen = cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('qvector.root')
)


process.QVector130 = process.QVector100.clone()
process.QVector160 = process.QVector100.clone()
process.QVector190 = process.QVector100.clone()
process.QVector220 = process.QVector100.clone()

process.QVector130.minCent = cms.untracked.int32(150)
process.QVector130.maxCent = cms.untracked.int32(185)
process.QVector160.minCent = cms.untracked.int32(185)
process.QVector160.maxCent = cms.untracked.int32(220)
process.QVector190.minCent = cms.untracked.int32(220)
process.QVector190.maxCent = cms.untracked.int32(260)
process.QVector220.minCent = cms.untracked.int32(260)
process.QVector220.maxCent = cms.untracked.int32(420)

process.p100 = cms.Path(process.hltHM100*process.QVector100)
process.p130 = cms.Path(process.hltHM130*process.QVector130)
process.p160 = cms.Path(process.hltHM160*process.QVector160)
process.p190 = cms.Path(process.hltHM190*process.QVector190)
process.p220 = cms.Path(process.hltHM220*process.QVector220)

process.schedule = cms.Schedule(
        process.p100,
        process.p130,
        process.p160,
        process.p190,
        process.p220,
)
