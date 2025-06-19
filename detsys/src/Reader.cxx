#include "Reader.h"

#include <filesystem>

template<typename T>
Reader<T>::Reader(std::string fname, std::string subfolder) {
  this->SetupTree();
  this->Open(fname, subfolder);
  _data = Data<T>();
  this->GetPOT();
  _entry = -1;
  _nentries = _Chain->GetEntries();
}

template<typename T>
Reader<T>::Reader(TChain *Chain) : _Chain(Chain) {
  _data = Data<T>();
  this->SetupTree();
  _entry = -1;
  _nentries = _Chain->GetEntries();
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
  
  std::string tree_name = "cafTree";
  std::string globtree_name = "meta";
  std::string genietree_name = "genieEvt";
  
  if(!subfolder.empty()){
    tree_name = subfolder + "/" + tree_name;
    globtree_name = subfolder + "/" + globtree_name;
    genietree_name = subfolder + "/" + genietree_name;
  }

  std::vector<std::string> ValidFilePaths;
  for (const auto & entry : std::filesystem::directory_iterator(fname)) {
    std::string filename = entry.path().string();

    TFile File = TFile(filename.c_str());
    if (File.IsZombie()) {
      std::cerr << "\t" << filename << " is not a valid root file" << std::endl;
      continue;
    }
    TTree* Tree = (TTree*)File.Get(tree_name.c_str());
    if (!Tree) {
      std::cerr << "\t" << filename << " does not contain" << tree_name << std::endl;
      continue;
    }

    ValidFilePaths.push_back(filename);
  }

  _Chain = new TChain(tree_name.c_str());
  for (auto filename : ValidFilePaths) {
    _Chain->Add(filename.c_str());
  }
  if (_Chain->GetEntries()==0) {
    std::cerr << "Not entries found - suggests invalid filename or corrupt file" << std::endl;
    abort();
  }
  _Chain->SetBranchAddress("rec", &_sr);

  _global_chain = new TChain(globtree_name.c_str());
  for (auto filename : ValidFilePaths) {
    _global_chain->Add(filename.c_str());
  }
  if (_global_chain->GetEntries()==0) {
    std::cerr << "Not entries found - suggests invalid filename or corrupt file" << std::endl;
    abort();
  }

}

template<typename T>
TFile* Reader<T>::GetFile(){
  return _file;
}

template<typename T>
void Reader<T>::SetupTree(){
  //_sr = new caf::StandardRecordProxy();
  _sr = new caf::StandardRecord();
}

template<typename T>
void Reader<T>::GetPOT(){
  double pot = 0;
  _global_chain->SetBranchAddress("pot", &pot);

  uint nentries = _global_chain->GetEntries();
  for(uint i = 0; i < nentries; i++){
    _global_chain->GetEntry(i);
    _POT += pot;
  } 
}

template<typename T>
bool Reader<T>::GetEntry(int i){
  bool retVal;
  if(i == -1){
    retVal = _Chain->GetEntry(++_entry) != 0;
  }
  else{
    retVal = _Chain->GetEntry(i) != 0;
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
  case kCVNNuMu:
    return _data.cvn_numu;
  case kCVNNuE:
    return _data.cvn_nue;
  }


  std::cerr << "Invalid kinematic parameter requested:" << Par << std::endl;
  throw;
}

template<typename T>
const T Reader<T>::GetEventWeight() {
  //return _data.weight / _POT;
  return 1;
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
    _data.RecoCZ = _BAD_VALUE_;
    _data.Selection = Unsel;
    _data.cvn_numu = _BAD_VALUE_;
    _data.cvn_nue = _BAD_VALUE_;
    return;
  }
  
  _data.Selection = Sel::SelNC;
  _data.cvn_numu = _sr->common.ixn.pandora[0].nuhyp.cvn.numu;
  _data.cvn_nue = _sr->common.ixn.pandora[0].nuhyp.cvn.nue;

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
TChain* Reader<T>::GetGlobalTree(){
  return _global_chain;
}

template<typename T>
TChain* Reader<T>::GetTree(){
  return _Chain;
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
