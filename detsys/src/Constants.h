#pragma once

#include <iostream>

enum Sel{
  SelNuE = 0,
  SelNuMu = 1,
  SelNC = 2,
  Unsel = 3,
  nSelections
};

enum Truth{
  Other = 0,
  CCNuE = 1,
  CCNuMu = 2,
  NC = 3
};

const double _BAD_VALUE_ = -999;

const double cvn_numu = 0.56;
const double cvn_nue = 0.55;

enum KinematicParameters {
  kNuRecoCosZ,
  kNuTrueCosZ,
  kCVNNuMu,
  kCVNNuE,
  kCVNNC,
  kSelection,
  kHadERec,
  kHadNuMuERec,
  kHadNuEERec,
  kCaloERec,
  kLepERec,
  kMuERec,
  keERec,
  kNuERec,
  kNuEERec,
  kNuMuERec,
  kNuETrue,
  kNuERes,
  kNuCosZRes,
  kNuCosZResAbs,
  kNuVertX,
  kNuVertY,
  kNuVertZ,
  kNuMomX,
  kNuMomY,
  kNuMomZ,
  kAnalysisBin,
  kTrueInt,
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
  if (Str == "kCVNNC") {
    return kCVNNC;
  }
  if (Str == "kSelection") {
    return kSelection;
  }
  if (Str == "kHadERec") {
    return kHadERec;
  }
  if (Str == "kHadNuMuERec") {
    return kHadNuMuERec;
  }
  if (Str == "kHadNuEERec") {
    return kHadNuEERec;
  }
  if (Str == "kCaloERec") {
    return kCaloERec;
  }
  if (Str == "kLepERec") {
    return kLepERec;
  }
  if (Str == "kMuERec") {
    return kMuERec;
  }
  if (Str == "keERec") {
    return keERec;
  }
  if (Str == "kNuERec") {
    return kNuERec;
  }
  if (Str == "kNuMuERec") {
    return kNuMuERec;
  }
  if (Str == "kNuEERec") {
    return kNuEERec;
  }
  if (Str== "kNuETrue") {
    return kNuETrue;
  }
  if (Str == "kNuERes") {
    return kNuERes;
  }
  if (Str == "kNuCosZRes") {
    return kNuCosZRes;
  }
  if (Str == "kNuCosZResAbs") {
    return kNuCosZRes;
  }
  if (Str == "kNuVertX") {
    return kNuVertX;
  }
  if (Str == "kNuVertY") {
    return kNuVertY;
  }
  if (Str == "kNuVertZ") {
    return kNuVertZ;
  }
  if (Str == "kNuMomX") {
    return kNuMomX;
  }
  if (Str == "kNuMomY") {
    return kNuMomY;
  }
  if (Str == "kNuMomZ") {
    return kNuMomZ;
  }
  if (Str == "kAnalysisBin") {
    return kAnalysisBin;
  }
  if (Str == "kTrueInt") {
    return kTrueInt;
  }

  std::cerr << "Did not find std::string -> int mapping for string:" << Str << std::endl;
  std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  throw;
}
