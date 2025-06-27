#include "HondaFluxReader.h"

#include <iostream>
#include <fstream>

HondaFluxReader::HondaFluxReader(YAML::Node Config) : FluxReader(Config) {
  if (!Config["FilePath"]) {
    std::cerr << "Did not find " << FilePath << " node in Config" << std::endl;
    throw;
  }                                                                                                                                                                                                                 
  FilePath = Config["FilePath"].as<std::string>();
  std::cout << "\t Found file path: " << FilePath << std::endl;

  InitialiseFlux();
}

std::vector< std::vector< std::vector<FLOAT_T> > > HondaFluxReader::ReadFlux() {
  std::vector< std::vector< std::vector<FLOAT_T> > > ReturnVec = std::vector< std::vector< std::vector<FLOAT_T> > >(GetNFlavours());

  FLOAT_T Energy = -999;
  FLOAT_T CosineZ = -999;
  FLOAT_T Phi = -999;
  FLOAT_T Flux_NuE = -999;
  FLOAT_T Flux_NuM = -999;
  FLOAT_T Flux_ANuE = -999;
  FLOAT_T Flux_ANuM = -999;
  
  std::ifstream input(FilePath);
  for(std::string line; getline( input, line );) {
    std::vector<std::string> Values_Str = SplitLine(line);
    if (Values_Str.size() == 0) continue;
    
    if (Values_Str.size() != 5 && Values_Str[0] == "average") {      

      if (Values_Str.size() == 13) {
	Values_Str[12].erase(std::remove(Values_Str[12].begin(), Values_Str[12].end(), ']'), Values_Str[12].end());
	
	FLOAT_T CosZ_LBinEdge = std::stod(Values_Str[5]);
        FLOAT_T CosZ_HBinEdge = std::stod(Values_Str[7]);
	FLOAT_T Phi_LBinEdge = std::stod(Values_Str[10]);
	FLOAT_T Phi_HBinEdge = std::stod(Values_Str[12]);

	CosineZ = (CosZ_HBinEdge+CosZ_LBinEdge)/2.;
	Phi = (Phi_HBinEdge+Phi_LBinEdge)/2.;

      }

      if (Values_Str.size() == 12) {
	Values_Str[4].erase(std::remove(Values_Str[4].begin(), Values_Str[4].end(), '='), Values_Str[4].end());
	Values_Str[11].erase(std::remove(Values_Str[11].begin(), Values_Str[11].end(), ']'), Values_Str[11].end());

	FLOAT_T CosZ_LBinEdge = std::stod(Values_Str[4]);
	FLOAT_T CosZ_HBinEdge = std::stod(Values_Str[6]);
	FLOAT_T Phi_LBinEdge = std::stod(Values_Str[9]);
	FLOAT_T Phi_HBinEdge = std::stod(Values_Str[11]);

	CosineZ = (CosZ_HBinEdge+CosZ_LBinEdge)/2.;
	Phi = (Phi_HBinEdge+Phi_LBinEdge)/2.;
      }
      
    }

    if (Values_Str.size() == 5) {
      Energy = std::stod(Values_Str[0]);
      FLOAT_T Energy3 = Energy*Energy*Energy;
      Flux_NuM = std::stod(Values_Str[1]) * Energy3;
      Flux_ANuM = std::stod(Values_Str[2]) * Energy3;
      Flux_NuE = std::stod(Values_Str[3]) * Energy3;
      Flux_ANuE = std::stod(Values_Str[4]) * Energy3;
      
      std::vector<FLOAT_T> FluxPoint_NuE = {Energy, CosineZ, Phi, Flux_NuE};
      std::vector<FLOAT_T> FluxPoint_NuM = {Energy, CosineZ, Phi, Flux_NuM};
      std::vector<FLOAT_T> FluxPoint_ANuE = {Energy, CosineZ, Phi, Flux_ANuE};
      std::vector<FLOAT_T> FluxPoint_ANuM = {Energy, CosineZ, Phi, Flux_ANuM};
      
      ReturnVec[NuE].push_back(FluxPoint_NuE);
      ReturnVec[NuM].push_back(FluxPoint_NuM);
      ReturnVec[ANuE].push_back(FluxPoint_ANuE);
      ReturnVec[ANuM].push_back(FluxPoint_ANuM);
    }
    
  }

  return ReturnVec;
}
