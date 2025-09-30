#pragma once
#include <vector>


template <typename T>
struct Data
{
  int Selection;
  int AnalysisBinIndex;
  T ev;
  T erec;
  T erec_nue;
  int mode;
  T NuMomX;
  T NuMomY;
  T NuMomZ;
  T LepMomX;
  T LepMomY;
  T LepMomZ;
  T cvnnumu;
  T cvnnue;
  int isCC;
  int nuPDGunosc;
  int nuPDG;
  int run;
  int isFHC;
  T BeRPA_cvwgt = 1;
  T vtx_x;
  T vtx_y;
  T vtx_z; 
  int ipnu;
  T LepTheta; 
  T Q2;
  T RecoCZ;
  T TrueCZ;
  T cvn_numu;
  T cvn_nue;
  double weight; //Have to hardcode that here for now.
  double genie_weight; //Have to hardcode that here for now.
  T flux_nue;
  T flux_numu;
  T xsec;
  T nue_w;
  T numu_w;
  T BeRPA_A_cvwgt = 1;
  T osc_from_e_w = 0;
  T osc_from_mu_w = 0;
  T final_oscillated_w = 0;
  std::vector<std::vector<float>> *weightVec;
};


static inline Data<double> float2double(const Data<float>& data){
  Data<double> data_D;
  data_D.Selection = data.Selection;
  data_D.AnalysisBinIndex = data.AnalysisBinIndex;
  data_D.ev = data.ev;
  data_D.erec = data.erec;
  data_D.erec_nue = data.erec_nue;
  data_D.mode = data.mode;
  data_D.NuMomX = data.NuMomX;
  data_D.NuMomY = data.NuMomY;
  data_D.NuMomZ = data.NuMomZ;
  data_D.LepMomX = data.LepMomX;
  data_D.LepMomY = data.LepMomY;
  data_D.LepMomZ = data.LepMomZ;
  data_D.cvnnumu = data.cvnnumu;
  data_D.cvnnue = data.cvnnue;
  data_D.isCC = data.isCC;
  data_D.nuPDGunosc = data.nuPDGunosc;
  data_D.nuPDG = data.nuPDG;
  data_D.run = data.run;
  data_D.isFHC = data.isFHC;
  data_D.BeRPA_cvwgt = data.BeRPA_cvwgt;
  data_D.vtx_x = data.vtx_x;
  data_D.vtx_y = data.vtx_y;
  data_D.vtx_z = data.vtx_z; 
  data_D.ipnu = data.ipnu;
  data_D.LepTheta = data.LepTheta; 
  data_D.Q2 = data.Q2;
  data_D.weight = data.weight;
  data_D.cvn_numu = data.cvn_numu;
  data_D.cvn_nue = data.cvn_nue;
  data_D.flux_nue = data.flux_nue;
  data_D.flux_numu = data.flux_numu;
  data_D.xsec = data.xsec;
  data_D.nue_w = data.nue_w;
  data_D.numu_w = data.numu_w;
  data_D.BeRPA_A_cvwgt = data.BeRPA_A_cvwgt;
  data_D.weightVec = data.weightVec;
  data_D.genie_weight = data.genie_weight;
  data_D.RecoCZ = data.RecoCZ;
  data_D.TrueCZ = data.TrueCZ;

  return data_D;
}
