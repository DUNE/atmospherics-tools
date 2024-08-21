#pragma once
#include "Data.h"
#include <string>
#include "TFile.h"
#include "TTree.h"

class Writer
{
protected:
    
    TTree *_tree;
    Data<double> _data;
    std::vector<TTree*> _additional_trees;

    void Create(std::string fname);
    void SetupTree();
public:
    TFile *_file;
    Writer();
    Writer(std::string fname);
    ~Writer();
    void Fill(const Data<double>& data);
    void Write();
    void WriteFile();
    TTree* GetTree();
    void AddNorm();
    void AddTree(TTree *tree);
    TFile* GetFile();
};
