#for file in TRK_newTracking MUO_newTracking MUO_updateTkMu MUO_updateOpenMu MUO_updateNoVtx MUO_newIO MUO_newOI BTV_noCalo_roiPF_DeepCSV BTV_noCalo_roiPF_DeepJet BTV_roiCalo_roiPF_DeepCSV BTV_roiCalo_roiPF_DeepJet BTV_roiCalo_globalPF_DeepCSV BTV_roiCalo_globalPF_DeepJet BTV_globalCalo_globalPF_DeepCSV BTV_globalCalo_globalPF_DeepJet MUO_useGEM MUO_newReco combined ;

options=" --globaltag auto:run2_data --data --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input --input file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics7/RAW/v1/000/323/790/00000/B543D251-40F1-CB46-A6A1-046CF3D78D6D.root --era Run2_2018 --output none --max-events 100 --setup /users/sdonato/CMSSW_12_2_0_GRunV5/combined/V6"


#for file in BTV_noCalo_roiPF_DeepCSV BTV_noCalo_roiPF_DeepJet BTV_roiCalo_globalPF_DeepCSV MUO_newIO MUO_newOI MUO_newReco MUO_newTracking MUO_updateNoVtx MUO_updateOpenMu MUO_updateTkMu MUO_useGEM TRK_newTracking combined nothing;
for file in combination;
do 
    config="/users/sdonato/CMSSW_12_2_0_GRunV5/"$file
    echo $config
    hltGetConfiguration $config $options > $file.py && edmConfigDump $file.py > $file\_dump.py && cmsRun $file\_dump.py >& $file\_dump.py.log &
done;

BTV_noCalo_roiPF_DeepCSV BTV_noCalo_roiPF_DeepJet BTV_roiCalo_globalPF_DeepCSV MUO_newIO MUO_newOI MUO_newReco MUO_newTracking MUO_updateNoVtx MUO_updateOpenMu MUO_updateTkMu MUO_useGEM TRK_newTracking combined nothing


