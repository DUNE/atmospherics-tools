#include <iostream>

#include "BartolFluxReader.h"
#include "HondaFluxReader.h"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: ./FluxComp Config.yaml" << std::endl;
    throw;
  }

  std::string ConfigName = argv[1];
  std::cout << "Config: " << ConfigName << std::endl;
  YAML::Node Config = YAML::LoadFile(ConfigName);

  if (!Config["FluxModels"]) {
    std::cerr << "No FluxModels to compare - Fix your YAML config" << std::endl;
    throw;
  }
  std::vector< std::string > ModelsFound;

  std::vector<FluxReader*> Fluxes;
  for (auto const &Model : Config["FluxModels"]) {

    std::string ModelName = Model["ModelName"].as<std::string>();
    if (std::find(ModelsFound.begin(), ModelsFound.end(), ModelName) != ModelsFound.end()) {
      std::cerr << "Already found model: " << ModelName << " in Config!" << std::endl;
      throw;
    }
    
    if (ModelName == "BARTOL" || ModelName == "bartol") {
      FluxReader* BARTOL = new BartolFluxReader(Model);
      Fluxes.push_back(BARTOL);
    } else if (ModelName == "HONDA" || ModelName == "honda") {
      FluxReader* HONDA = new HondaFluxReader(Model);
      Fluxes.push_back(HONDA);
    } else {
      std::cerr << "Found ModelName: " << ModelName << " which has not been associated with a FluxReader" << std::endl;
      throw;
    }

    ModelsFound.push_back(ModelName);
  }

  for (auto Flux: Fluxes) {
    Flux->Plot2DFlux(Flux->GetModelName()+".pdf","COLZ");
  }
  
  return 0;
}
