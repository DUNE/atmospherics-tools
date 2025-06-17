#include "Reader.h"

template<typename T>
Reader<T>::Reader(std::string fname, std::string subfolder) {
  this->Open(fname, subfolder);
  _data = Data<T>();
  this->SetupTree();
  this->GetPOT();
  _entry = -1;
  _nentries = _tree->GetEntries();
}

template<typename T>
Reader<T>::Reader(TTree *tree) : _tree(tree) {
  _data = Data<T>();
  this->SetupTree();
  _entry = -1;
  _nentries = _tree->GetEntries();
}

template<typename T>
Reader<T>::~Reader() {
  if(_file){
    _file->Close();
  }
  if(_sr){
    delete _sr;
  }
}

template<typename T>
void Reader<T>::Open(std::string fname, std::string subfolder) {
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
  
  //Loading the global tree
  
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
const T Reader<T>::ReturnKinematicParameter(int Par) {
  switch (Par) {
  case kNuTrueCosZ:
    return _data.TrueCZ;
  case kNuRecoCosZ:
    return _data.RecoCZ;
  }

  std::cerr << "Invalid kinematic parameter requested:" << Par << std::endl;
  throw;
}

template<typename T>
const T Reader<T>::GetEventWeight() {
  return _data.weight;
}

template<typename T>
void Reader<T>::UpdateData(){
  
  _data.ev = _sr->mc.nu[0].E;
  _data.NuMomX = _sr->mc.nu[0].momentum.x;
  _data.NuMomY = _sr->mc.nu[0].momentum.y;
  _data.NuMomZ = _sr->mc.nu[0].momentum.z;
  _data.nuPDG = _sr->mc.nu[0].pdg;
  _data.weight = _sr->mc.nu[0].genweight;
  _data.mode = _sr->mc.nu[0].mode;

  _data.TrueCZ = -_sr->mc.nu[0].momentum.y;
  
  if(_sr->common.ixn.pandora.size() != 1){
    _data.erec = _BAD_VALUE_;
    _data.Selection = Unsel;
    return;
  }
  
  _data.Selection = Sel::SelNC;
  if(_sr->common.ixn.pandora[0].nuhyp.cvn.numu > cvn_numu){
    _data.Selection = Sel::SelNuMu;
  }
  else if(_sr->common.ixn.pandora[0].nuhyp.cvn.nue > cvn_nue){
    _data.Selection = Sel::SelNuE;
  }
  
  if (_data.Selection == Sel::SelNC) {
    _data.erec = _sr->common.ixn.pandora[0].Enu.calo;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.heshw.y;
  } else if (_data.Selection == Sel::SelNuMu) {
    _data.erec = _sr->common.ixn.pandora[0].Enu.lep_calo;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.lngtrk.y;
  } else if (_data.Selection == Sel::SelNuE) {
    _data.erec = _sr->common.ixn.pandora[0].Enu.e_calo;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.heshw.y;
  } else {
    std::cerr << "Invalid selection" << std::endl;
    throw;
  }


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

template<typename T>
int Reader<T>::GetNentries(){
  return _nentries;
};

template class Reader<float>;
template class Reader<double>;
