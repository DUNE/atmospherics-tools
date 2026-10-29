#pragma once

using FLOAT_T = double;

enum Flavours{
  NuE = 0,
  NuM = 1,
  ANuE = 2,
  ANuM = 3,
  nFlavs = 4
};

enum FlavRatios{
  Ratio_Total, // (NuMu+ANuMu)/(NuE+ANuE)
  Ratio_NuE, // NuE/ANuE
  Ratio_NuM, // NuMu/ANuMu
  nFlavRatios
};

const std::vector< std::string > FlavourNames = {"NuE","NuM","ANuE","ANuM"};
const std::vector< std::string > RatioFlavourNames = {"(NuM+ANuM)/(NuE+ANuE)","NuM/ANuM","NuE/ANuE"};
const std::vector< std::string > RatioFlavourShortNames = {"TotalRatio","NuMRatio","NuERatio"};
