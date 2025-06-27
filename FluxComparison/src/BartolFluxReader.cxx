#include "BartolFluxReader.h"

#include <iostream>
#include <fstream>

BartolFluxReader::BartolFluxReader(YAML::Node Config) : FluxReader(Config) {
  std::cout << "Finding file path names:" << std::endl;
  for (int iFlav=0;iFlav<GetNFlavours();iFlav++) {
    std::string FlavourName = GetFlavourName(iFlav);
    if (!Config[FlavourName]) {
      std::cerr << "Did not find " << FlavourName << " node in Config" << std::endl;
      throw;
    }                                                                                                                                                                                                                 
    FilePathMap[FlavourName] = Config[FlavourName].as<std::string>();
    std::cout << "\t Found file path for flavour (" << FlavourName << " : " << FilePathMap[FlavourName] << std::endl;
  }                                                                                                                                                                                                                   
  std::cout << std::endl;

  InitialiseFlux();
}

std::vector< std::vector< std::vector<FLOAT_T> > > BartolFluxReader::ReadFlux() {
  std::vector< std::vector< std::vector<FLOAT_T> > > ReturnVec;

  for (int iFlav=0;iFlav<GetNFlavours();iFlav++) {
    std::vector< std::vector<FLOAT_T> > FlavourVec;
      
    std::string FlavourName = GetFlavourName(iFlav);
    std::string FilePath = FilePathMap[FlavourName];
    
    std::cout << "\tReading BARTOL flux (" << FlavourName << ") : " << FilePath << std::endl;
    
    std::ifstream input(FilePath);
    for(std::string line; getline( input, line );) {
      std::vector<std::string> Values_Str = SplitLine(line);
      if (Values_Str[0] == "#") continue;
      
      if (Values_Str.size() != 5) {
	std::cerr << "Expected format of BARTOL flux is 5 entries per line: Energy CosZ Flux MCStatUnc nUnweightMC" << std::endl;
	throw;
      }
      
      std::vector<FLOAT_T> Values(Values_Str.size());
      for (size_t i=0;i<Values.size();i++) {
	Values[i] = std::stod(Values_Str[i]);
      }
      
      FLOAT_T Energy_BinCentre = Values[0];
      FLOAT_T CosZ_BinCentre = Values[1];
      FLOAT_T FluxValue = Values[2];

      //Convert dn/dlogE [m^-2 x str^-1 x s^-1] to dn/dE [m^-2 x str^-1 x s^-1 * GeV^2]
      FluxValue *= Energy_BinCentre*Energy_BinCentre;
      /*
      int XBin = FluxHists[iFlav]->GetXaxis()->FindBin(Energy_BinCentre);
      FLOAT_T BinWidth = FluxHists[iFlav]->GetXaxis()->GetBinWidth(XBin);
      FLOAT_T Log_BinWidth = log(FluxHists[iFlav]->GetXaxis()->GetBinLowEdge(XBin)+FluxHists[iFlav]->GetXaxis()->GetBinWidth(XBin)) - log(FluxHists[iFlav]->GetXaxis()->GetBinLowEdge(XBin));
      FluxValue *= ((Log_BinWidth/BinWidth) * BinWidth * BinWidth);
      */
      
      std::vector<FLOAT_T> Point = {Energy_BinCentre,CosZ_BinCentre,FluxValue};
      FlavourVec.push_back(Point);
    }

    ReturnVec.push_back(FlavourVec);
  }

  return ReturnVec;
}
