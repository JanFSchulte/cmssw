# MuonAnalysis-MuonAnalyzer

Package to run tag/probe ntuples for Muon POG on both AOD and miniAOD format.

Requires CMSSW_12_4_8 or higher.

## Setup

For Run2023:

```bash
$ cmsrel CMSSW_13_0_5_patch2 
$ cd CMSSW_13_0_5_patch2/src
$ cmsenv
$ git cms-init
$ git clone https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer.git MuonAnalysis/MuonAnalyzer -b 12_4_X
$ scram b -j 8
```

For Run2022:

```bash
$ cmsrel CMSSW_12_4_8 
$ cd CMSSW_12_4_8/src
$ cmsenv
$ git cms-init
$ git clone https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer.git MuonAnalysis/MuonAnalyzer -b 12_4_X
$ scram b -j 8
```

## Usage
```bash
$ cmsRun MuonAnalysis/MuonAnalyzer/test/run_muonAnalyzer_cfg.py
```

## Contributing

Before creating a Merge Request, make sure the code is properly formatted with:

```bash
$ find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\)' | xargs clang-format -i
```

Run this from the base directory. It should go through all source code files in the repository and format them according to CMSSW's style (the `Google` style).


### Gitlab CI setup for Developers (required to submit merge requests!!!)

1. Add predefined variables

    i. from your fork go to `Muonanalysis-Muonanalyzer > settings > CI/CD`

    ii. expand the `Variables` section

    iii. click the `Add variable` button

        - add key: USER_NAME and value: <your CERN user name>

    iv. click the `Add variable` button

        - select the flag `Mask variable`
        - add key: USER_PASS and value: <your CERN password encoded to base64>
          e.g encode "thisword": printf "thisword" | base64
