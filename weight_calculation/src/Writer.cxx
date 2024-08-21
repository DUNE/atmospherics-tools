#include "Writer.h"
#include "TH1D.h"
#include <iostream>

Writer::Writer()
{
}

Writer::Writer(std::string fname)
{
    this->Create(fname);
    this->SetupTree();
}

Writer::~Writer()
{
    if(_file)
        _file->Close();
}

void Writer::AddNorm(){
    _file->cd();
    TH1D *norm = new TH1D("norm","norm",10,0,10);
    norm->SetDirectory(_file);
    // Set norm bin 1 content (neutrino type normalisation)
    double val = 1e21;
    norm->SetBinContent(2,val);
    norm->SetBinContent(1,6e-5);
    norm->Write();
    delete norm;
}

void Writer::Create(std::string fname)
{
    std::cout << "Opening output file: " << fname << std::endl;
    _file = TFile::Open(fname.c_str(), "RECREATE");
    if(_file == nullptr){
        std::cout << "Could not open file: " << fname << std::endl;
        abort();
    }
    _tree = new TTree("cafTree", "cafTree");
    _tree->SetDirectory(_file);
}

void Writer::Fill(const Data<double>& data){
    _data = data;
    _tree->Fill();
}

void Writer::SetupTree(){
    _tree->Branch("Ev", &_data.ev);
    _tree->Branch("Ev_reco_numu", &_data.erec);
    _tree->Branch("Ev_reco_nue", &_data.erec_nue);
    _tree->Branch("mode",&_data.mode);
    _tree->Branch("NuMomX", &_data.NuMomX);
    _tree->Branch("NuMomY", &_data.NuMomY);
    _tree->Branch("NuMomZ", &_data.NuMomZ);
    _tree->Branch("LepMomX", &_data.LepMomX);
    _tree->Branch("LepMomY", &_data.LepMomY);
    _tree->Branch("LepMomZ", &_data.LepMomZ);
    _tree->Branch("cvnnumu", &_data.cvnnumu);
    _tree->Branch("cvnnue", &_data.cvnnue);
    _tree->Branch("isCC", &_data.isCC);
    _tree->Branch("nuPDGunosc", &_data.nuPDGunosc);
    _tree->Branch("nuPDG", &_data.nuPDG);
    _tree->Branch("run", &_data.run);
    _tree->Branch("isFHC", &_data.isFHC);
    _tree->Branch("BeRPA_A_cvwgt", &_data.BeRPA_cvwgt);
    _tree->Branch("vtx_x", &_data.vtx_x);
    _tree->Branch("vtx_y", &_data.vtx_y);
    _tree->Branch("vtx_z", &_data.vtx_z);  
    _tree->Branch("LepPDG",&_data.ipnu); 
    _tree->Branch("LepNuAngle",&_data.LepTheta); 
    _tree->Branch("Q2",&_data.Q2);
    _tree->Branch("weight",&_data.weight);
    _tree->Branch("flux_nue",&_data.flux_nue);
    _tree->Branch("flux_numu",&_data.flux_numu);
    _tree->Branch("xsec",&_data.xsec);
    _tree->Branch("nue_w",&_data.nue_w);
    _tree->Branch("numu_w",&_data.numu_w);
    _tree->Branch("osc_from_e_w",&_data.osc_from_e_w);
    _tree->Branch("osc_from_mu_w",&_data.osc_from_mu_w);
    _tree->Branch("final_oscillated_w",&_data.final_oscillated_w);
    _tree->Branch("genie_weight",&_data.genie_weight);

    //Current location of systs weights
    _tree->Branch("xsSyst_wgt",&_data.weightVec);
}

void Writer::Write(){
    // _tree->SetDirectory(_file);
    std::cout << "Nentries: " << _tree->GetEntries() << std::endl;
    _tree->Write();

    for(TTree *tree : _additional_trees){
        tree->Write();
    }
}

void Writer::WriteFile(){
    // _tree->SetDirectory(_file);
    _file->Write();
}

TTree* Writer::GetTree(){
    _file->cd();
    TTree* tree = static_cast<TTree*>(_tree->CloneTree());
    // tree->SetDirectory(0);
    return tree;
}

void Writer::AddTree(TTree *tree){
    if(tree != nullptr){
        // TTree* new_tree = tree->CloneTree();
        tree->SetDirectory(_file);
        _additional_trees.push_back(tree);
    }
    else{
        std::cout << "ERROR: Cannot add a null tree to Writer! Fix your code!" << std::endl;
        abort();
    }
}

TFile* Writer::GetFile(){
    return _file;
}