import re,sys
from CRABClient.UserUtilities import config
config = config()

tag = "ntuple_aod"
project = "unprefirable_egamma"

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_egamma.py'
config.JobType.outputFiles = ['output.root']

# disable ASO
# https://github.com/dmwm/CRABServer/blob/master/test/templates/config/MinBias_PrivateMC_EventBased_ExtraParams.py#L45
# config.section_("Debug")
# config.Debug.extraJDL = ['+CRAB_StageoutPolicy="remote"'] # %%stageout%%

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/Crab3DataHandling
#config.Data.inputDataset = '/EGamma/Run2018D-12Nov2019_UL2018-v8/AOD'
#config.Data.inputDataset = '/EGamma/Run2018C-12Nov2019_UL2018-v2/AOD'
#config.Data.inputDataset = '/EGamma/Run2018B-12Nov2019_UL2018-v2/AOD'
#config.Data.inputDataset = '/EGamma/Run2018A-12Nov2019_UL2018-v2/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2017B-09Aug2019_UL2017-v1/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2017C-09Aug2019_UL2017-v1/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2017D-09Aug2019_UL2017-v1/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2017E-09Aug2019_UL2017-v1/AOD'
config.Data.inputDataset = '/SingleElectron/Run2017F-09Aug2019_UL2017_EcalRecovery-v1/AOD'
##/SingleMuon/Run2018A-17Sep2018-v2/AOD
##/SingleMuon/Run2018B-17Sep2018-v1/AOD
##/SingleMuon/Run2018C-17Sep2018-v1/AOD
##/SingleMuon/Run2018D-22Jan2019-v2/AOD

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.publication = False
#config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/group/phys_muon/jschulte'

config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = [ 'T2_CH_CERN' ]

#    'T2_CH_CERN','T2_US_MIT','T2_BE_IIHE','T2_US_Vanderbilt','T2_US_Caltech',
#    'T2_US_UCSD','T2_US_Nebraska','T2_US_Purdue','T2_EE_Estonia',
#    'T2_US_Wisconsin','T2_CH_CSCS','T2_US_Florida','T2_DE_DESY','T1_UK_RAL']


# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
match = re.search('^\/[^\/]+\/([^\/]+)\/',config.Data.inputDataset)
if match:
    config.General.requestName = '%s_%s_%s_v2' % (tag,match.group(1),project)
else:
    sys.exit(1)
# config.General.requestName = 'V08-00-03_Run2016B-23Sep2016-v3_muon-cleaning'
config.General.transferLogs = True

