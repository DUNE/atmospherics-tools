#pragma once
#include "yaml-cpp/yaml.h"

#include "Utils.h"
#include "Constants.h"

#include "TH1.h"

#include <vector>

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
  int MeasDimension;
  
  FluxReader(YAML::Node Config);
  void InitialiseFlux();
  virtual std::vector< std::vector< std::vector<FLOAT_T> > > ReadFlux() = 0;
  void Build2DPlots();
  void BuildFlavourRatioPlots();
  
  std::vector<TH1*> FluxHists = std::vector<TH1*>(nFlavs);
  std::vector<TH1*> EnergyCosineZHists = std::vector<TH1*>(nFlavs);
  std::vector<TH1*> FlavourRatioHists = std::vector<TH1*>(nFlavRatios);
  YAML::Node Config;
  
public:

  void Plot2DFlux(std::string OutputName, std::string DrawOpts);
  std::string GetModelName() {return ModelName;}
  std::vector<TH1*> ReturnEnergyCosineZHists() {return EnergyCosineZHists;}
  std::vector<TH1*> ReturnFlavourRatioHists() {return FlavourRatioHists;}
};
