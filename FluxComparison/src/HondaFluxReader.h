#pragma once

#include "yaml-cpp/yaml.h"

#include "FluxReader.h"

class HondaFluxReader : public FluxReader {
private:
  std::string FilePath;
  
protected:
  std::vector< std::vector< std::vector<FLOAT_T> > > ReadFlux() override;
  
public:
  HondaFluxReader(YAML::Node Config);
};
