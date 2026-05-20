#include "Reader.h"

#include <filesystem>

template<typename T>
Reader<T>::Reader(std::string fname, std::string subfolder) {
  _fname = fname;
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

template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
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
  case kCVNNC:
    return _data.cvn_nc;
  case kSelection:
    return _data.Selection;
  case kNuETrue:
    return _data.ev;
  case kHadERec:
    return _data.had_erec;
  case kHadNuMuERec:
    return _data.had_numu_erec;
  case kHadNuEERec:
    return _data.had_nue_erec;
  case kCaloERec:
    return _data.calo_erec;
  case kLepERec:
    return _data.lep_erec;
  case kMuERec:
    return _data.mu_erec;
  case keERec:
    return _data.e_erec;
  case kNuERec:
    return _data.erec;
  case kNuMuERec:
    return _data.erec_numu;
  case kNuEERec:
    return _data.erec_nue;
  case kNuERes:
    return (_data.erec - _data.ev)/_data.ev;
  case kNuCosZRes:
    return (_data.RecoCZ - _data.TrueCZ)/_data.TrueCZ;
  case kNuCosZResAbs:
    return (_data.RecoCZ - _data.TrueCZ);
  case kNuVertX:
    return _data.vtx_x;
  case kNuVertY:
    return _data.vtx_y;
  case kNuVertZ:
    return _data.vtx_z;
  case kNuMomX:
    return _data.NuMomX;
  case kNuMomY:
    return _data.NuMomY;
  case kNuMomZ:
    return _data.NuMomZ;
  case kAnalysisBin:
    return _data.AnalysisBinIndex;
  case kTrueInt:
    return _data.trueInt;
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
  
  T ENu = _sr->mc.nu[0].E;
  T CosTh = _sr->mc.nu[0].momentum.y / sqrt(pow(_sr->mc.nu[0].momentum.x,2) + pow(_sr->mc.nu[0].momentum.y,2) + pow(_sr->mc.nu[0].momentum.z,2));
  T Phi = 60.; //Same as here: https://github.com/DUNE/atmospherics-tools/blob/41b9bb7afe248f9923aab6c3cc52900e76d2d3c1/weight_calculation/src/Calculator.cxx#L38
  
  Flavour nuE = (_sr->mc.nu[0].pdg > 0) ? Flavour::NuE : Flavour::NuEBar;
  Flavour nuMu = (_sr->mc.nu[0].pdg > 0) ? Flavour::NuMu : Flavour::NuMuBar;

  double nuE_flux = FlxMgr->GetFlux(nuE, ENu, CosTh, Phi);
  double nuMu_flux = FlxMgr->GetFlux(nuMu, ENu, CosTh, Phi);
  double ref_flux = FlxMgr->GetFlux(Flavour::Reference, ENu, CosTh, Phi);

  int NuSignSwitch = (_sr->mc.nu[0].pdg > 0) ? 1 : -1;
  int NuOscFlavIndex;
  switch(abs(_sr->mc.nu[0].pdg)) {
  case 12:
    NuOscFlavIndex = 1; break;
  case 14:
    NuOscFlavIndex = 2; break;
  case 16:
    NuOscFlavIndex = 3; break;
  default:
    throw;
  }
		      
  T Oscillation_FromNumu = OscillBase->ReturnOscillationProbability(2*NuSignSwitch,NuOscFlavIndex*NuSignSwitch,ENu,CosTh);
  T Oscillation_FromNuE  = OscillBase->ReturnOscillationProbability(1*NuSignSwitch,NuOscFlavIndex*NuSignSwitch,ENu,CosTh);

  _data.weight = (_sr->mc.nu[0].genweight/ref_flux/_POT*(400.0*3600.0*24.0*365.0/1.71958)) * (Oscillation_FromNuE*nuE_flux + Oscillation_FromNumu*nuMu_flux);
  //_data.weight = 1.0;

  _data.AnalysisBinIndex = -1;
  _data.ev = _sr->mc.nu[0].E;
  _data.NuMomX = _sr->mc.nu[0].momentum.x;
  _data.NuMomY = _sr->mc.nu[0].momentum.y;
  _data.NuMomZ = _sr->mc.nu[0].momentum.z;
  _data.isCC = _sr->mc.nu[0].iscc;
  _data.nuPDG = _sr->mc.nu[0].pdg;
  _data.mode = _sr->mc.nu[0].mode;

  _data.trueInt = Truth::NC;
  if (_data.isCC){
    if (std::abs(_data.nuPDG)==12) _data.trueInt = Truth::CCNuE;
    else if (std::abs(_data.nuPDG)==14) _data.trueInt = Truth::CCNuMu;
  }

  TVector3 TrueNuMomentumVector = (TVector3(_sr->mc.nu[0].momentum.X(),_sr->mc.nu[0].momentum.Y(),_sr->mc.nu[0].momentum.Z())).Unit();
  _data.TrueCZ = -TrueNuMomentumVector.Y();
  
  if(_sr->common.ixn.pandora.size() != 1){
    _data.had_erec = _BAD_VALUE_;
    _data.lep_erec = _BAD_VALUE_;
    _data.erec = _BAD_VALUE_;
    _data.had_numu_erec = _BAD_VALUE_;
    _data.mu_erec = _BAD_VALUE_;
    _data.erec_numu = _BAD_VALUE_;
    _data.had_nue_erec = _BAD_VALUE_;
    _data.e_erec = _BAD_VALUE_;
    _data.erec_nue = _BAD_VALUE_;
    _data.RecoCZ = _BAD_VALUE_;
    _data.Selection = Unsel;
    _data.cvn_numu = _BAD_VALUE_;
    _data.cvn_nue = _BAD_VALUE_;
    _data.trueInt = Truth::Other;
    //std::cout<<"unselected because of pandora size: "<<_sr->common.ixn.pandora.size()<<std::endl;
    return;
  }
 


  // TODO: fix this! Make this a mandatory argument in the yaml file
  if(_fname == "./data/Deprecated/Recombination_n1Sig"){ 
      _data.cvn_numu = _sr->common.ixn.pandora[0].nuhyp.cvn.numu;
  }
  else{
      _data.cvn_numu = _sr->common.ixn.pandora[0].nuhyp.cvn.nc;
  }

  //_data.cvn_nc = _sr->common.ixn.pandora[0].nuhyp.cvn.nc;
  //_data.cvn_numu = _sr->common.ixn.pandora[0].nuhyp.cvn.numu;
  _data.cvn_nue = _sr->common.ixn.pandora[0].nuhyp.cvn.nue;


  //std::vector<T> CVNScores = {_sr->common.ixn.pandora[0].nuhyp.cvn.nue, _sr->common.ixn.pandora[0].nuhyp.cvn.numu, _sr->common.ixn.pandora[0].nuhyp.cvn.nc};
  //_data.Selection = arg_max(CVNScores);
  

  _data.Selection = Sel::SelNC;
  if(_data.cvn_numu > cvn_numu){
    _data.Selection = Sel::SelNuMu;
  }
  else if(_data.cvn_nue > cvn_nue){
    _data.Selection = Sel::SelNuE;
  }

  if (_data.Selection == Sel::SelNC) {
    _data.had_erec = _sr->common.ixn.pandora[0].Enu.calo;
    _data.erec = _sr->common.ixn.pandora[0].Enu.calo;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.heshw.y;
  } else if (_data.Selection == Sel::SelNuMu) {
    _data.erec = _sr->common.ixn.pandora[0].Enu.lep_calo;
    _data.had_erec = _sr->common.ixn.pandora[0].Enu.mu_had;
    _data.lep_erec = _data.erec - _data.had_erec;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.lngtrk.y;
  } else if (_data.Selection == Sel::SelNuE) {
    _data.erec = _sr->common.ixn.pandora[0].Enu.e_calo;
    _data.had_erec = _sr->common.ixn.pandora[0].Enu.e_had;
    _data.lep_erec = _data.erec - _data.had_erec;
    _data.RecoCZ = -_sr->common.ixn.pandora[0].dir.heshw.y;
  } else {
    std::cerr << "Invalid selection" << std::endl;
    throw;
  }

  _data.had_numu_erec =  _sr->common.ixn.pandora[0].Enu.mu_had;
  _data.had_nue_erec =  _sr->common.ixn.pandora[0].Enu.e_had;
  _data.erec_numu =  _sr->common.ixn.pandora[0].Enu.lep_calo;
  _data.erec_nue =  _sr->common.ixn.pandora[0].Enu.e_calo;
  _data.mu_erec = _data.erec_numu - _data.had_numu_erec;
  _data.e_erec = _data.erec_nue - _data.had_nue_erec;
  _data.calo_erec = _sr->common.ixn.pandora[0].Enu.calo;

  if (std::isnan(_data.RecoCZ)) {
    return;
  }

  if (AnalysisBinning->CheckSelectionInAnalysisBinning(_data.Selection)) {
    std::vector<std::string> AnalysisBinningVars = AnalysisBinning->GetSelectionBinVars(_data.Selection);
    
    std::vector<T> EventDetails(AnalysisBinningVars.size());
    bool validEvent = true;

    for (size_t iBinVar=0;iBinVar<AnalysisBinningVars.size();iBinVar++) {
      EventDetails[iBinVar] = ReturnKinematicParameter(Kinematic_StringToInt(AnalysisBinningVars[iBinVar]));
      
      if (!std::isfinite(EventDetails[iBinVar])) {
        std::cerr << "Skipping event due to invalid kinematic parameter: "
                      << AnalysisBinningVars[iBinVar]
                      << " = " << EventDetails[iBinVar] << std::endl;
        validEvent = false;
        break;
      }
      if (EventDetails[iBinVar] < -1e8 || EventDetails[iBinVar] >1e8){
        std::cerr << "Skipping event due to invalid kinematic parameter: "
                      << AnalysisBinningVars[iBinVar]
                      << " = " << EventDetails[iBinVar] << std::endl;
        validEvent = false;
        break;
      } 
    }
    if (validEvent) {
      _data.AnalysisBinIndex = AnalysisBinning->GetBin(_data.Selection, EventDetails);
    } else {
      _data.AnalysisBinIndex = -1;  // mark as invalid
      _data.weight = 0;
      return;  // skip further processing
    }    
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
