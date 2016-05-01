import FWCore.ParameterSet.Config as cms

process = cms.Process("QVector")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v13', '')

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

#fN = cms.untracked.vstring();
#for line in open('flist').read().splitlines():
#	fN.append('file:'+line);
#
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/q/qwang/work/cleanroomRun2/Ana/CMSSW_7_5_8_patch2/src/QWAna/QWCumuV3/test/HIMinBias_28.root")
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

process.hltMB = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMB.HLTPaths = [
	"HLT_HIL1MinimumBiasHF2AND*",
	"HLT_HIL1MinimumBiasHF1AND*",
]
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)


process.QVector = cms.EDAnalyzer('QWQVector'
	, centrality = cms.InputTag("centralityBin", "HFtowers")
	, trackTag = cms.untracked.InputTag('hiGeneralTracks')
	, vertexSrc = cms.untracked.InputTag('hiSelectedVertex', "")
	, pterrorpt = cms.untracked.double(0.1)
	, dzdzerror = cms.untracked.double(3.0)
	, d0d0error = cms.untracked.double(3.0)
	, minvz = cms.untracked.double(-1.0)
	, maxvz = cms.untracked.double(15.0)
	, minEta = cms.untracked.double(-2.4)
	, maxEta = cms.untracked.double(2.4)
	, minPt = cms.untracked.double(1.0)
	, maxPt = cms.untracked.double(3.0)
	, minCent = cms.untracked.int32(-1)
	, maxCent = cms.untracked.int32(500)
	, epSrc = cms.untracked.InputTag("hiEvtPlane")
#	, fweight_ = cms.untracked.InputTag('PbPb_dijet_TT_5TeV_v2.root')
#	, bEff_ = cms.untracked.bool(False)
	, algoParameters = cms.vint32(4,5,6,7)
	, bGen = cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('qvector.root')
)

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.clusterCompatibilityFilter.clusterPars = cms.vdouble(0.0,0.006)

process.eventSelection = cms.Sequence(
        process.hfCoincFilter3
        + process.primaryVertexFilter
        + process.clusterCompatibilityFilter
)

process.path= cms.Path(process.hltMB*process.eventSelection*process.centralityBin*process.QVector)

process.schedule = cms.Schedule(
	process.path,
)
