#!/bin/env python
import multiprocessing,re,subprocess,commands
import os,re,hashlib

unmerged_path = "/tmp/jschulte/prefiring/"

ncpus = 16

def find_root_files(path):
    return commands.getoutput("find %s -type f -name '*root'"%path).split("\n")

eras = {
     'Run2016B':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016B-07Aug17_ver2-v1_unprefirable/190714_073234'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016C':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016C-07Aug17-v1_unprefirable/190714_073319'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016D':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016D-07Aug17-v1_unprefirable/190714_073354'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016E':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016E-07Aug17-v1_unprefirable/190716_044107'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016F':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016F-07Aug17-v1_unprefirable/190714_073451'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016G':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016G-07Aug17-v1_unprefirable/190726_094125'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2016H':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016H-07Aug17-v1_unprefirable/190714_073545'),
         'tag':'80X_dataRun2_2016SeptRepro_v7'
         },
     'Run2017B':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017B-17Nov2017-v1_unprefirable/190820_154737'),
         'tag':'94X_dataRun2_v11'
         },    
     'Run2017C':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017C-17Nov2017-v1_unprefirable/190820_154812'),
         'tag':'94X_dataRun2_v11'
         },    
     'Run2017D':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017D-17Nov2017-v1_unprefirable/190819_154116'),
         'tag':'94X_dataRun2_v11'
         },    
     'Run2017E':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017E-17Nov2017-v1_unprefirable/190728_065825'),
         'tag':'94X_dataRun2_v11'
         },    
     'Run2017F':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017F-17Nov2017-v1_unprefirable/190820_154833'),
         'tag':'94X_dataRun2_v11'
         },    
     'Run2018A':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018A-17Sep2018-v2_unprefirable/190725_130436'),
         'tag':'102X_dataRun2_v11'
         },    
     'Run2018B':{
         'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018B-17Sep2018-v1_unprefirable/190715_070339'),
         'tag':'102X_dataRun2_v11'
         },    
    'Run2018C':{
        'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018C-17Sep2018-v1_unprefirable/190728_064656'),
        'tag':'102X_dataRun2_v11'
        },    
    'Run2018D':{
        'files':find_root_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018D-22Jan2019-v2_unprefirable/190726_103909'),
        'tag':'102X_dataRun2_v11'
        },    
    }


jobs = dict()
                
def update_job_status(job):
    processed = False
    prefix = "%s/%s/%s" % (unmerged_path,job['era'],job['hash_name'])
    if os.path.exists("%s.root" % prefix) and os.path.exists("%s.done" % prefix):
        job['processed'] = True
    else:
        job['processed'] = False

def process_job(job):
    if job['processed']: return True
    pid = multiprocessing.current_process().pid
    prefix = "%s/%s/%s" % (unmerged_path,job['era'],job['hash_name'])
    files = ",".join(job['inputFiles'])
    command = "cmsRun prefireStudyAODtoNtuple.py inputFiles='%s' globaltag='%s' outputFile='%s.root' >& '%s.log'" % (files,job['tag'],prefix,prefix)
    exit_code = subprocess.call(command,shell=True)
    if exit_code==0:
        subprocess.call("touch %s.done"%prefix,shell=True)

for era,info in eras.items():
    print "Checkings jobs for %s" % era
    if os.path.exists("%s.root"%era):
        print "Merged output found. Nothing to process. Skip it."
        continue
    jobs[era] = []
    tmp_era_path = "%s/%s" % (unmerged_path,era)
    if not os.path.exists(tmp_era_path):
        os.makedirs(tmp_era_path)
    for file in info['files']:
        job = {
            'hash_name':hashlib.md5(file).hexdigest(),
            'tag':info['tag'],
            'inputFiles':[file],
            'era':era
            }
        update_job_status(job)
        jobs[era].append(job)
    print "Njobs (%s) : %u " % (era,len(jobs[era]))

pool = multiprocessing.Pool(ncpus)
for era in jobs:
    print "Processing jobs for era: %s" % era
    nProcessed = 0
    for job in jobs[era]:
        if job['processed']: nProcessed+=1
    print "\tCurrent status: processed %u out of %u" % (nProcessed,len(jobs[era]))
    pool.map(process_job, jobs[era])
    nProcessed = 0
    for job in jobs[era]:
        update_job_status(job)
        if job['processed']: nProcessed+=1
    print "\tCycle is completed: processed %u out of %u" % (nProcessed,len(jobs[era]))
    if nProcessed==len(jobs[era]):
        output_file_name = "%s.root" % era
        if not os.path.exists(output_file_name):
            command = "hadd %s %s/%s/*.root" % (output_file_name,unmerged_path,era)
            print command
            subprocess.call(command,shell=True)
