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

template <typename T>
class Reader
{
private:
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
    Reader(std::string fname, std::string subfolder = "");
    Reader(TTree *tree);
    ~Reader();
    double POT(){return _POT;};
    bool GetEntry(int i = -1);
    const Data<T>& GetData();
    TTree* GetGlobalTree();
    TTree* GetGenieTree();
    TTree* GetTree();
    TFile* GetFile();
};