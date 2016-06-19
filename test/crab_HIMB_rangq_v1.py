from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

config = config()

config.General.requestName = 'HIMB3_QVector_randq_v1'
config.General.workArea = 'CrabArea'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'qw_PbPb15_HIMB_randq.py'
#config.JobType.inputFiles = ['PbPb_dijet_TT_5TeV_v2.root']
config.Data.inputDataset =  '/HIMinimumBias3/qwang-HIMinBias_v2-ce439b1c24fa1bf3a491f2ccb0fd72a9/USER'
#config.Data.inputDataset = '/HIMinimumBias3/HIRun2015-PromptReco-v1/AOD'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/group/phys_heavyions/qwang/QVector/'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/Cert_262548-263757_PromptReco_HICollisions15_JSON_v2.txt'
config.Data.publication = False
config.Data.useParent = False
config.Site.storageSite = 'T2_CH_CERN'
#try:
#        crabCommand('submit', config = config)
#except HTTPException as hte:
#        print "Failed submitting task: %s" % (hte.headers)
#except ClientException as cle:
#        print "Failed submitting task: %s" % (cle)


# HIMB1
config.General.requestName = 'HIMB1_QVector_randq_v1_hlt'
config.Data.inputDataset =  '/HIMinimumBias1/HIRun2015-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
try:
        crabCommand('submit', config = config)
except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

# HIMB2
config.General.requestName = 'HIMB2_QVector_randq_v1'
config.Data.inputDataset =  '/HIMinimumBias2/HIRun2015-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
#try:
#        crabCommand('submit', config = config)
#except HTTPException as hte:
#        print "Failed submitting task: %s" % (hte.headers)
#except ClientException as cle:
#        print "Failed submitting task: %s" % (cle)

# HIMB4
config.General.requestName = 'HIMB4_QVector_randq_v1'
config.Data.inputDataset =  '/HIMinimumBias4/HIRun2015-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
#try:
#        crabCommand('submit', config = config)
#except HTTPException as hte:
#        print "Failed submitting task: %s" % (hte.headers)
#except ClientException as cle:
#        print "Failed submitting task: %s" % (cle)

