#pragma once
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <memory>

#include "Data.h"
#include "Constants.h"

// duneanaobj
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "AnalysisBinningManager.h"
#include "FluxManager.h"
#include "Oscillator/OscillatorFactory.h"

template <typename T>
class Reader {
 private:
  std::string _fname = "";
  TFile *_file = nullptr;
  TChain *_Chain = nullptr;
  TChain *_global_chain = nullptr;
  TTree *_genie_tree = nullptr;
  double _POT = 0;
  Data<T> _data;
  int _nentries;
  int _entry = -1;
  //caf::StandardRecordProxy* _sr = nullptr;
  caf::StandardRecord* _sr = nullptr;

  FluxManager* FlxMgr;
  AnalysisBinningManager<T>* AnalysisBinning;
  OscillatorBase* OscillBase;

  void Open(std::string fname, std::string subfolder);
  void SetupTree();
  void SetupTreeHierarchical();
  void GetPOT();
  void UpdateData(std::string fname);
  
 public:
  Reader(std::string fname, std::string subfolder = "");
  Reader(TChain *Chain);
  ~Reader();
  double POT(){return _POT;};
  bool GetEntry(int i = -1);
  const Data<T>& GetData();
  TChain* GetGlobalTree();
  TTree* GetGenieTree();
  TChain* GetTree();
  TFile* GetFile();
  int GetNentries();
  void SetAnalysisBinning(AnalysisBinningManager<T>* AnalysisBinning_) {AnalysisBinning = AnalysisBinning_;}
  void SetOscillator(OscillatorBase* OscillBase_) {OscillBase = OscillBase_;}
  void SetFluxManager(FluxManager* FlxMgr_) {FlxMgr = FlxMgr_;}
  const T ReturnKinematicParameter(int Par);
  const T GetEventWeight();
};
