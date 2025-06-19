#include "ObservableManager.h"
#include <iostream>

#include "Constants.h"

template<typename T>
Observable<T>::Observable(YAML::Node ObservableConfig) {
  Name = ObservableConfig["Name"].as<std::string>();
  
  nDimensions = 0;
  for (auto AxisNode: ObservableConfig["Axes"]) {
    nDimensions += 1;

    Axis<T> A = Axis<T>();

    A.Variable = AxisNode["Variable"].as<std::string>();
    A.Variable_Int = Kinematic_StringToInt(A.Variable);
    A.Label = AxisNode["Label"].as<std::string>();

    if (AxisNode["Binning"]) {
      A.Binning = AxisNode["Binning"].as<std::vector<T>>();
    } else if (AxisNode["BinDefinition"]) {
      std::vector<T> BinDef = AxisNode["BinDefinition"].as<std::vector<T>>();
      if (BinDef.size()!=3) {
	std::cerr << "BinDefinition Node found. Expected 3 items [nBins,LowEdge,HighEdge]" << std::endl;
	throw;
      }

      int nBins = (int)BinDef[0];
      T LowBinEdge = BinDef[1];
      T HighBinEdge = BinDef[2];

      std::vector<T> Binning;
      for (int iBin=0;iBin<nBins;iBin++) {
	Binning.push_back(LowBinEdge+(T)iBin*(HighBinEdge-LowBinEdge)/(T)nBins);
      }
      Binning.push_back(HighBinEdge);

      A.Binning = Binning;
    }
    
    Axes.emplace_back(A);
  }

  std::cout << "Observable:" << std::endl;
  std::cout << "\tName:" << Name << std::endl;
  std::cout << "\tnDimensions:" << nDimensions << std::endl;

  for (int i=0;i<nDimensions;i++) {
    std::cout << "\t\tDimension:"<< i << std::endl;
    
    std::cout << "\t\tVariable:" << Axes[i].Variable << "(" << Axes[i].Variable_Int << ")" << std::endl;
    std::cout << "\t\tLabel:" << Axes[i].Label << std::endl;

    std::cout << "\t\tBinning: [";
    for (auto BinEdge: Axes[i].Binning) {std::cout << BinEdge << ", ";}
    std::cout << "]" << std::endl;    
  }

  if (nDimensions == 0) {
    std::cerr << "Invalid number of axes defined for observable:" << Name << std::endl;
    throw;
  } else if (nDimensions == 1) {
    std::string HistName = Name+"_0";
    std::string HistTitle = Axes[0].Label+";"+Axes[0].Label+";Number of Events";
    int nBins = Axes[0].Binning.size()-1;
    
    if (typeid(T) == typeid(float)) {
      HistTemplate = new TH1F(HistName.c_str(),HistTitle.c_str(),nBins,Axes[0].Binning.data());
    } else {
      HistTemplate = new TH1D(HistName.c_str(),HistTitle.c_str(),nBins,Axes[0].Binning.data());
    }
  } else if (nDimensions == 2) {
    std::string HistName = Name+"_0";
    std::string HistTitle = Axes[0].Label+";"+Axes[1].Label+";Number of Events";
    int nXBins = Axes[0].Binning.size()-1;
    int nYBins = Axes[1].Binning.size()-1;
    
    if (typeid(T) == typeid(float)) {
      HistTemplate = new TH2F(HistName.c_str(),HistTitle.c_str(),nXBins,Axes[0].Binning.data(),nYBins,Axes[1].Binning.data());
    } else {
      HistTemplate = new TH2D(HistName.c_str(),HistTitle.c_str(),nXBins,Axes[0].Binning.data(),nYBins,Axes[1].Binning.data());
    }

  } else {
    std::cerr << "Currently only have support for 1 and 2 dimension observables" << std::endl;
    throw;
  }

  HistTemplate->SetDirectory(0);
  std::cout << std::endl;
}

template<typename T>
ObservableManager<T>::ObservableManager(YAML::Node Config) {
  for (auto ObservableNode: Config["Observables"]) {
    Observables.emplace_back(Observable<T>(ObservableNode));
  }
}
