#pragma once

#include "yaml-cpp/yaml.h"

#include "FluxReader.h"

class BartolFluxReader : public FluxReader {
private:
  std::unordered_map<std::string, std::string> FilePathMap;
  
protected:
  std::vector< std::vector< std::vector<FLOAT_T> > > ReadFlux() override;
  
public:
  BartolFluxReader(YAML::Node Config);
};
