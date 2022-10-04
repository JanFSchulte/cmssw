from __future__ import print_function
import ROOT as r 
import sys 
from collections import defaultdict

r.gROOT.SetBatch(True)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print(bcolors.OKCYAN, "========== Comparing", bcolors.ENDC, sys.argv[1], bcolors.OKCYAN, 'and', bcolors.ENDC, sys.argv[2], bcolors.OKCYAN, "==========", bcolors.ENDC)

reference = r.TFile.Open(sys.argv[1])
target   =  r.TFile.Open(sys.argv[2])
suffix   = sys.argv[3]

ref_tree=reference.Get("muon/Events")
tar_tree=target   .Get("muon/Events")

print(bcolors.OKCYAN, '[INFO]: Checking if the branches contained are the same', bcolors.ENDC)

ref_branches = [ br.GetName() for br in ref_tree.GetListOfBranches()]
tar_branches = [ br.GetName() for br in tar_tree.GetListOfBranches()]

fail=False

for ref_br in ref_branches:
    if ref_br not in tar_branches:
        print(bcolors.WARNING, "[FAIL]: Branch", bcolors.ENDC, bcolors.OKGREEN, ref_br, bcolors.ENDC, bcolors.WARNING, 'is removed by the PR', bcolors.ENDC)
        fail=True

for tar_br in tar_branches:
    if tar_br not in ref_branches:
        print(bcolors.WARNING, "[FAIL]: Branch", bcolors.ENDC, bcolors.OKGREEN, tar_br, bcolors.ENDC, bcolors.WARNING, 'is added by the PR', bcolors.ENDC)
        fail=True

if fail: 
    print(bcolors.WARNING, "[FAIL]: This PR changes the branches in the tree. We will continue with the checks just in case there are other differences.", bcolors.ENDC)
else:
    print(bcolors.OKCYAN, "[INFO]: This PR contains the exact same branches as in the reference", bcolors.ENDC)


print(bcolors.OKCYAN, "[INFO]: Now we will compare event by event", bcolors.ENDC)
## Continue with the checks so other differences may be spotted
allBranches=[tar_br for tar_br in tar_branches if tar_br in ref_branches]

ref_events=[]
for ev in ref_tree:
    ref_events.append( dict( (br, getattr(ev, br)) for br in allBranches))

tar_events=[]
for ev in tar_tree:
    tar_events.append( dict( (br, getattr(ev, br)) for br in allBranches))

residuals=defaultdict(list)

def matchEvents(ref_ev, tar_ev):
    if ref_ev['event'] != tar_ev['event']: return False
    if ref_ev['lumi'] != tar_ev['lumi']: return False
    if ref_ev['run'] != tar_ev['run']: return False
    if ref_ev['run'] != tar_ev['run']: return False
    if abs(ref_ev['tag_pt'] - tar_ev['tag_pt'] ) > 0.01: return False
    if abs(ref_ev['tag_eta'] - tar_ev['tag_eta'] ) > 0.001: return False
    if abs(ref_ev['probe_pt'] - tar_ev['probe_pt'] ) > 0.01: return False
    if abs(ref_ev['probe_eta'] - tar_ev['probe_eta'] ) > 0.001: return False
    return True

print(bcolors.OKCYAN, '[INFO]: Checking if events in the reference are also in the target ', bcolors.ENDC)

branchesToSkip=[]
for ref_ev in ref_events:
    isMatched=False
    for tar_ev in tar_events:
        if not matchEvents(ref_ev, tar_ev): continue
        isMatched=True
        for br in allBranches:
            if br in branchesToSkip: continue # so its no so verbose
            try:
                if "cppyy.gbl.std.vector" in str(type(tar_ev[br])) or "cppyy.gbl.std.vector" in str(type(ref_ev[br])):
                    branchesToSkip.append(br)
                    continue
                
                residuals[br].append( tar_ev[br] - ref_ev[br] )
            except TypeError:
                print(bcolors.FAIL, '[WARNING]: Branches of type ', type(tar_ev[br]),'like', br, 'not supported for event by event comparison', bcolors.ENDC)
                branchesToSkip.append(br)
        break
    if not isMatched:
        print(bcolors.WARNING, '[FAIL]: Tag-probe pair in event %d:%d:%d tag_pt=%4.2f,tag_eta=%4.2f,probe_pt=%4.2f,eta_pt=%4.2f removed by the PR'%(ref_ev['event'], ref_ev['lumi'],ref_ev['run'], ref_ev['tag_pt'], ref_ev['tag_eta'], ref_ev['probe_pt'], ref_ev['probe_eta']), bcolors.ENDC)
        fail=True

print(bcolors.OKCYAN, '[INFO]: Checking if variable definition is the same ', bcolors.ENDC)

differentVariables=[]
c=r.TCanvas()
for i,var in enumerate(residuals):
    if not len(residuals[var]): continue
    mini=min(residuals[var])
    maxi=max(residuals[var])
    xmin=mini-(maxi-mini)/100
    xmax=maxi+(maxi-mini)/100
    h=r.TH1F(var, '', 100, xmin, xmax)
    h.GetXaxis().SetTitle("Target - reference")
    for what in residuals[var]:
        h.Fill(what)
    if abs(mini) > 1e-7 or abs(maxi) > 1e-7:
        h.GetXaxis().SetAxisColor(r.kRed)
        h.GetYaxis().SetAxisColor(r.kRed)
        h.SetLineColor(r.kRed)
        h.SetLineWidth(2)
        fail=True
        differentVariables.append(var)
    h.Draw()
    c.SetLogy(True)
    if not i: 
        c.Print("../../residuals%s.pdf("%suffix)
    if i==len(residuals)-1:
        c.Print("../../residuals%s.pdf)"%suffix)
    else:
        c.Print("../../residuals%s.pdf"%suffix)

for var in differentVariables:
    print(bcolors.WARNING, "[FAIL]: Variable %s has significant differences between reference and target"%var, bcolors.ENDC)
    


for tar_ev in tar_events:
    isMatched=False
    for ref_ev in ref_events:
        if not matchEvents(ref_ev, tar_ev): continue
        isMatched=True
    if not isMatched:
        print(bcolors.WARNING, '[FAIL]: Tag-probe pair in event %d:%d:%d tag_pt=%4.2f,tag_eta=%4.2f,probe_pt=%4.2f,eta_pt=%4.2f added by the PR'%(tar_ev['event'], tar_ev['lumi'],tar_ev['run'], tar_ev['tag_pt'], tar_ev['tag_eta'], tar_ev['probe_pt'], tar_ev['probe_eta']), bcolors.ENDC)
        fail=True

if fail: 
    print(bcolors.FAIL, '[FAIL]: Some of the checks above have failed. Please make sure the differences are intended and send the muon conveners the new reference together with this log so they can update the reference file. ', bcolors.ENDC)
else:
    print(bcolors.OKCYAN, '[INFO]: Tests have passed. This either means that a) your PR doesnt change anything, b) you are changing something that is not caught by our CI (action item?) or c) you have understood the differences and have updated the reference file ', bcolors.ENDC)
    
assert(not fail)
