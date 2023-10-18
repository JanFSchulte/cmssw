//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// flat tree branches/var declaration

#ifndef NTUPLECONTENT_H
#define NTUPLECONTENT_H
#include <string>
#include <vector>
#include "TString.h"
#include "TTree.h"
#include <variant>

using AllTypes = std::variant<bool, int, float, double, long,
unsigned int, long unsigned int, long long unsigned int,
std::vector<float>, std::vector<int>,std::vector<bool>>;

template<typename... Ts>
struct getVal {
    std::variant<Ts...> & var;
    template <typename T>
    operator T() { return std::get<T>(var); }
};

class NtupleContent {
public:
  NtupleContent();
  virtual ~NtupleContent();
  void SetTree(TTree *t1);
  void CreateBranches(const std::vector<std::string> &, const std::vector<std::string> &);
  void CreateExtraTrgBranches(const std::vector<std::string> &, bool);
  void ClearBranches();

  AllTypes & operator()(std::string key) {
    return branches[key].value;
  }

  typedef struct BranchInfo {
    BranchInfo() {}
    BranchInfo(AllTypes in_val) :
      default_value(in_val), value(in_val) {}
    BranchInfo& operator= (AllTypes val) { value = val; return * this; }
    AllTypes & operator()() { return value; }
    template<typename T, typename... Ts>
    friend std::ostream & operator<<(std::ostream& os, const std::variant<T, Ts...>& v) {
      std::visit([&os](auto&& arg) {
                   os << arg;
                 }, v);
      return os;
    }

    AllTypes default_value;
    AllTypes value;
  } BranchInfo;

  std::map<TString, BranchInfo> branches;


  // Triggers
  static const int NTRIGGERMAX = 100;
  bool trigger[NTRIGGERMAX];
  std::vector<TString> trg_filter;
  std::vector<float> trg_pt;
  std::vector<float> trg_eta;
  std::vector<float> trg_phi;

  std::vector<TString> prb_filter;
  std::vector<float> prb_pt;
  std::vector<float> prb_eta;
  std::vector<float> prb_phi;

  // Trigger matches
  bool tag_trg[NTRIGGERMAX];
  float tag_trg_pt[NTRIGGERMAX];
  float tag_trg_eta[NTRIGGERMAX];
  float tag_trg_phi[NTRIGGERMAX];
  float tag_trg_dr[NTRIGGERMAX];
  bool probe_trg[NTRIGGERMAX];
  float probe_trg_pt[NTRIGGERMAX];
  float probe_trg_eta[NTRIGGERMAX];
  float probe_trg_phi[NTRIGGERMAX];
  float probe_trg_dr[NTRIGGERMAX];

  // Standard selectors in reco::muon::Selector
  bool probe_selectors[100];

private:
  TTree *t1;
};
#endif
