#include <iostream>
#include "SRWriter.h"

SRWriter::SRWriter(std::string name, caf::StandardRecord *sr)
{
    _sr = sr;
    this->Create(name);
    this->SetupTree();
}

SRWriter::~SRWriter()
{
}

void SRWriter::Create(std::string name)
{
    std::cout << "Opening output file: " << name << std::endl;
    _file = TFile::Open(name.c_str(), "RECREATE");
    if(_file == nullptr){
        std::cout << "Could not open file: " << name << std::endl;
        abort();
    }
    _tree = new TTree("cafTree", "cafTree");
    _tree->SetDirectory(_file);
}

void SRWriter::SetupTree(){
    _tree->Branch("rec", &_sr);
}

void SRWriter::Fill(){
    _tree->Fill();
}
void SRWriter::SetGenieTree(TTree *tree){
    genieTree = tree;
    genieTree->SetDirectory(_file);
    _additional_trees.push_back(genieTree);
}

TTree* SRWriter::GetGenieTree(){
    return genieTree;
}