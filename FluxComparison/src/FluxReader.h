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

enum FlavRatios{
  Ratio_Total, // (NuMu+ANuMu)/(NuE+ANuE)
  Ratio_NuE, // NuE/ANuE
  Ratio_NuM, // NuMu/ANuMu
  nFlavRatios
};

class FluxReader {

private:

  std::string ModelName;
  std::vector<FLOAT_T> EnergyBinEdges;
  std::vector<FLOAT_T> CosineZBinEdges;
  std::vector<FLOAT_T> PhiBinEdges;

  // [Flavour][Point][Energy,CosZ,{Phi},Flux]
  std::vector< std::vector< std::vector< FLOAT_T > > > FluxPoints;

  FLOAT_T EnergyAxisMin;
  FLOAT_T EnergyAxisMax;
  int LineColor;
  int LineStyle;
  std::string FluxCaption;
  bool Smooth = false;
  
protected:
  const int NuE = 0;
  const int NuM = 1;
  const int ANuE = 2;
  const int ANuM = 3;
  const int nFlavours = 4;
  const std::vector< std::string > FlavourNames = {"NuE","NuM","ANuE","ANuM"};
  const std::vector< std::string > RatioFlavourNames = {"(NuM+ANuM)/(NuE+ANuE)","NuM/ANuM","NuE/ANuE"};
  
  int MeasDimension;
  
  FluxReader(YAML::Node Config);
  void InitialiseFlux();
  virtual std::vector< std::vector< std::vector<FLOAT_T> > > ReadFlux() = 0;
  void Build2DPlots();
  void BuildFlavourRatioPlots();
  
  std::vector<TH1*> FluxHists = std::vector<TH1*>(nFlavours);
  std::vector<TH1*> EnergyCosineZHists = std::vector<TH1*>(nFlavours);
  std::vector<TH1*> FlavourRatioHists = std::vector<TH1*>(nFlavRatios);
  YAML::Node Config;
  
public:
  void Plot2DFlux(std::string OutputName, std::string DrawOpts);

  std::string GetModelName() {return ModelName;}
  std::string GetFlavourName(int Flav) {return FlavourNames.at(Flav);}
  std::string GetFlavourRatioName(int Flav) {return RatioFlavourNames.at(Flav);}
  int GetNFlavours() {return nFlavours;}

  std::vector<TH1*> ReturnEnergyCosineZHists() {return EnergyCosineZHists;}
  std::vector<TH1*> ReturnFlavourRatioHists() {return FlavourRatioHists;}
};
