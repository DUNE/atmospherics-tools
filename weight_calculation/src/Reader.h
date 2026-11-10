#pragma once
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <memory>
#include "Data.h"
// duneanaobj
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

// Input formats understood by the Reader:
//  - CAF : hierarchical CAF (cafTree/meta/genieEvt), read through SRProxy
//  - Flat: flat ntuple produced from LArSoft files by larsoft/extract_larsoft.C
//          (events/meta trees, see that macro for the branch content)
enum class InputFormat { CAF, Flat };

// Branch buffers for the flat input format
struct FlatInput
{
    int run = 0;
    int subrun = 0;
    int eid = 0;
    int nuPDG = 0;
    double Ev = 0;
    double NuMomX = 0;
    double NuMomY = 0;
    double NuMomZ = 0;
    double genie_weight = 0;
};

template <typename T>
class Reader
{
private:
    InputFormat _format = InputFormat::CAF;
    FlatInput _flat;
    std::string _fname = "";
    TFile *_file = nullptr;
    TTree *_tree = nullptr;
    TTree *_global_tree = nullptr;
    TTree *_genie_tree = nullptr;
    double _POT = 0;
    Data<T> _data;
    int _nentries;
    int _entry = -1;
    caf::StandardRecordProxy* _sr = nullptr;

    void Open(std::string fname, std::string subfolder);
    void SetupTree();
    void SetupTreeHierarchical();
    void GetPOT();
    void UpdateData();

public:
    Reader(std::string fname, std::string subfolder = "", InputFormat format = InputFormat::CAF);
    Reader(TTree *tree);
    ~Reader();
    double POT(){return _POT;};
    bool GetEntry(int i = -1);
    const Data<T>& GetData();
    TTree* GetGlobalTree();
    TTree* GetGenieTree();
    TTree* GetTree();
    TFile* GetFile();
    int GetNentries();
    InputFormat GetFormat(){return _format;};
};