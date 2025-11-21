#pragma once

#include <iostream>

enum Sel{
  SelNuE = 0,
  SelNuMu = 1,
  SelNC = 2,
  Unsel = 3,
  nSelections
};

const double _BAD_VALUE_ = -999;

const double cvn_numu = 0.56;
const double cvn_nue = 0.55;

enum KinematicParameters {
  kNuRecoCosZ,
  kNuTrueCosZ,
  kCVNNuMu,
  kCVNNuE,
  kSelection,
  kNuERec,
  kNuETrue,
  kNuERes,
  kAnalysisBin,
  nKinPars
};

inline int Kinematic_StringToInt(std::string Str) {
  if (Str == "kNuRecoCosZ") {
    return kNuRecoCosZ;
  }
  if (Str == "kNuTrueCosZ") {
    return kNuTrueCosZ;
  }
  if (Str == "kCVNNuMu") {
    return kCVNNuMu;
  }
  if (Str == "kCVNNuE") {
    return kCVNNuE;
  }
  if (Str == "kSelection") {
    return kSelection;
  }
  if (Str == "kNuERec") {
    return kNuERec;
  }
  if (Str== "kNuETrue") {
    return kNuETrue;
  }
  if (Str == "kNuERes") {
    return kNuERes;
  }
  if (Str == "kAnalysisBin") {
    return kAnalysisBin;
  }

  std::cerr << "Did not find std::string -> int mapping for string:" << Str << std::endl;
  std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  throw;
}
