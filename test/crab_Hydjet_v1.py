from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

config = config()

config.General.requestName = 'Hydjet_QVector_v5'
config.General.workArea = 'CrabArea'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'qw_PbPb15_Hydjet.py'
#config.JobType.inputFiles = ['PbPb_dijet_TT_5TeV_v2.root']
config.Data.inputDataset =  '/Hydjet_Quenched_MinBias_5020GeV_750/HINPbPbWinter16DR-NoPU_75X_mcRun2_HeavyIon_v1_75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/HIMinimumBias3/HIRun2015-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/group/phys_heavyions/qwang/QVector/'
config.Data.publication = False
config.Data.useParent = False
config.Site.storageSite = 'T2_CH_CERN'
try:
        crabCommand('submit', config = config)
except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
except ClientException as cle:
        print "Failed submitting task: %s" % (cle)


