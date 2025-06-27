#pragma once
#include "yaml-cpp/yaml.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <vector>
#include <string>

using FLOAT_T = double;

inline std::vector<std::string> SplitLine(const std::string& str) {
  std::vector<std::string> tokens;
  std::istringstream iss(str);
  std::string token;
  
  while (iss >> token) {
    tokens.push_back(token);
  }
  
  return tokens;
}

enum Flavours{
  NuE,
  NuM,
  ANuE,
  ANuM,
  nFlavs
};

class FluxReader {

private:

  std::string ModelName;
  std::vector<FLOAT_T> EnergyBinEdges;
  std::vector<FLOAT_T> CosineZBinEdges;
  std::vector<FLOAT_T> PhiBinEdges;

  // [Flavour][Point][Energy,CosZ,{Phi},Flux]
  std::vector< std::vector< std::vector< FLOAT_T > > > FluxPoints;

protected:
  const int NuE = 0;
  const int NuM = 1;
  const int ANuE = 2;
  const int ANuM = 3;
  const int nFlavours = 4;

  int MeasDimension;
  
  FluxReader(YAML::Node Config);
  void InitialiseFlux();
  virtual std::vector< std::vector< std::vector<FLOAT_T> > > ReadFlux() = 0;

  const std::vector< std::string > FlavourNames = {"NuE","NuM","ANuE","ANuM"};
  std::vector<TH1*> FluxHists = std::vector<TH1*>(nFlavours);
  YAML::Node Config;
  
public:
  void Plot2DFlux(std::string OutputName, std::string DrawOpts);

  std::string GetModelName() {return ModelName;}
  std::string GetFlavourName(int Flav) {return FlavourNames.at(Flav);}
  int GetNFlavours() {return nFlavours;}
};
