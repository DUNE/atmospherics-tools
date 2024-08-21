#include "Reader.h"

template<typename T>
Reader<T>::Reader(std::string fname, std::string subfolder)
{
   this->Open(fname, subfolder);
   _data = Data<T>();
   this->SetupTree();
   this->GetPOT();
   _entry = -1;
   _nentries = _tree->GetEntries();
}

template<typename T>
Reader<T>::Reader(TTree *tree)
: _tree(tree)
{
   _data = Data<T>();
   this->SetupTree();
   _entry = -1;
   _nentries = _tree->GetEntries();
}

template<typename T>
Reader<T>::~Reader()
{
    if(_file){
        _file->Close();
    }
    if(_sr){
        delete _sr;
    }
}

template<typename T>
void Reader<T>::Open(std::string fname, std::string subfolder)
{
    std::cout << "Opening file: " << fname << std::endl;
    _file = TFile::Open(fname.c_str(), "READ");
    if(_file == nullptr){
        std::cout << "Could not open file: " << fname << std::endl;
        abort();
    }

    std::string tree_name = "cafTree";
    std::string globtree_name = "meta";
    std::string genietree_name = "genieEvt";

    if(!subfolder.empty()){
        tree_name = subfolder + "/" + tree_name;
        globtree_name = subfolder + "/" + globtree_name;
        genietree_name = subfolder + "/" + genietree_name;
    }

    TTree *tree = nullptr;
    
    _file->GetObject(tree_name.c_str(), tree);
    if(tree == nullptr){
        std::cout << "No tree named " << tree_name << " in file" << std::endl;
        abort();
    }
    _tree = tree;

    //Loading the gloabl tree
    
    TTree *global_tree;
    _file->GetObject(globtree_name.c_str(), global_tree);
    if(global_tree == nullptr){
        std::cout << "WARNING: No tree named " << globtree_name << " in the input file " << fname << std::endl;
    }
    _global_tree = global_tree;

    TTree *genie_tree;
    _file->GetObject(genietree_name.c_str(), genie_tree);
    if(genie_tree == nullptr){
        std::cout << "WARNING: No tree named " << genietree_name << " in the input file " << fname << std::endl;
    }
    _genie_tree = genie_tree;
}

template<typename T>
TFile* Reader<T>::GetFile(){
    return _file;
}

template<typename T>
void Reader<T>::SetupTree(){
    _sr = new caf::StandardRecordProxy(_tree, "rec");
    // TString start_str = "";
    // _tree->SetBranchAddress(start_str + "Ev", &_data.ev);
    // _tree->SetBranchAddress(start_str + "Ev_reco_numu", &_data.erec);
    // _tree->SetBranchAddress(start_str + "Ev_reco_nue", &_data.erec_nue);
    // _tree->SetBranchAddress(start_str + "mode",&_data.mode);
    // _tree->SetBranchAddress(start_str + "NuMomX", &_data.NuMomX);
    // _tree->SetBranchAddress(start_str + "NuMomY", &_data.NuMomY);
    // _tree->SetBranchAddress(start_str + "NuMomZ", &_data.NuMomZ);
    // _tree->SetBranchAddress(start_str + "LepMomX", &_data.LepMomX);
    // _tree->SetBranchAddress(start_str + "LepMomY", &_data.LepMomY);
    // _tree->SetBranchAddress(start_str + "LepMomZ", &_data.LepMomZ);
    // _tree->SetBranchAddress(start_str + "cvnnumu", &_data.cvnnumu);
    // _tree->SetBranchAddress(start_str + "cvnnue", &_data.cvnnue);
    // _tree->SetBranchAddress(start_str + "isCC", &_data.isCC);
    // _tree->SetBranchAddress(start_str + "nuPDGunosc", &_data.nuPDGunosc);
    // _tree->SetBranchAddress(start_str + "nuPDG", &_data.nuPDG);
    // _tree->SetBranchAddress(start_str + "run", &_data.run);
    // _tree->SetBranchAddress(start_str + "isFHC", &_data.isFHC);
    // _tree->SetBranchAddress(start_str + "BeRPA_A_cvwgt", &_data.BeRPA_cvwgt);
    // _tree->SetBranchAddress(start_str + "vtx_x", &_data.vtx_x);
    // _tree->SetBranchAddress(start_str + "vtx_y", &_data.vtx_y);
    // _tree->SetBranchAddress(start_str + "vtx_z", &_data.vtx_z);  
    // _tree->SetBranchAddress(start_str + "LepPDG",&_data.ipnu); 
    // _tree->SetBranchAddress(start_str + "LepNuAngle",&_data.LepTheta); 
    // _tree->SetBranchAddress(start_str + "Q2",&_data.Q2);
    
    // //Current location of systs weights
    // _tree->SetBranchAddress("xsSyst_wgt",&_data.weightVec);

    // if(true){
    //     _tree->SetBranchAddress(start_str + "pot",&_data.weight); //The weight is stored here in a hacky way at the moment.
    // }
    // else{
    //     _tree->SetBranchAddress(start_str + "weight",&_data.weight);
    //     _tree->SetBranchAddress("genie_weight",&_data.genie_weight); //The weight is stored here in a hacky way at the moment.
    //     _tree->SetBranchAddress("nue_w",&_data.nue_w);
    //     _tree->SetBranchAddress("numu_w",&_data.numu_w);
    //     _tree->SetBranchAddress("xsec",&_data.xsec);
    // }
    
}

template<typename T>
void Reader<T>::GetPOT(){
    double pot = 0;
    _global_tree->SetBranchAddress("pot", &pot);
    uint nentries = _global_tree->GetEntries();
    for(uint i = 0; i < nentries; i++){
        _global_tree->GetEntry(i);
        _POT += pot;
    } 
}

template<typename T>
bool Reader<T>::GetEntry(int i){
    bool retVal;
    if(i == -1){
        retVal = _tree->GetEntry(++_entry) != 0;
    }
    else{
        retVal = _tree->GetEntry(i) != 0;
    }

    if(retVal){
        UpdateData();
    }

    return retVal;
}

template<typename T>
const Data<T>& Reader<T>::GetData(){
    return _data;
}

template<typename T>
void Reader<T>::UpdateData(){
    _data.ev = _sr->mc.nu[0].E;
    _data.NuMomX = _sr->mc.nu[0].momentum.x;
    _data.NuMomY = _sr->mc.nu[0].momentum.y;
    _data.NuMomZ = _sr->mc.nu[0].momentum.z;
    _data.nuPDG = _sr->mc.nu[0].pdg;
    _data.weight = _sr->mc.nu[0].genweight;
}

template<typename T>
TTree* Reader<T>::GetGlobalTree(){
    return _global_tree;
}

template<typename T>
TTree* Reader<T>::GetTree(){
    return _tree;
}

template<typename T>
TTree* Reader<T>::GetGenieTree(){
    return _genie_tree;
}

template class Reader<float>;
template class Reader<double>;
