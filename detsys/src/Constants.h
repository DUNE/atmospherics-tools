#pragma once

#include <iostream>

enum Sel{
  SelNuE,
  SelNuMu,
  SelNC,
  Unsel
};

const double _BAD_VALUE_ = -999;

const double cvn_numu = 0.56;
const double cvn_nue = 0.55;

enum KinematicParameters {
  kNuRecoCosZ,
  kNuTrueCosZ,
  kCVNNuMu,
  kCVNNuE,
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

  std::cerr << "Did not find std::string -> int mapping for string:" << Str << std::endl;
  std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  throw;
}
