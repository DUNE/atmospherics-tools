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

    A.isLog = false;
    if (AxisNode["IsLog"]) {
      A.isLog = AxisNode["IsLog"].as<bool>();
    }

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
    if (AxisNode["BinLabels"]){
      A.BinLabels = AxisNode["BinLabels"].as<std::vector<std::string>>();
      if (A.BinLabels.size()!=A.Binning.size()-1){
        std::cerr << "Number of bin labels does not match number of bins for variable "<<A.Variable<<std::endl;
        throw;
      }
    }
    Axes.emplace_back(A);
  }

  if (ObservableConfig["Cuts"]) {
    
    for (auto CutNode: ObservableConfig["Cuts"]) {
      Cut<T> ObsCut = Cut<T>();
      
      ObsCut.Variable = CutNode["Variable"].as<std::string>();
      ObsCut.Variable_Int = Kinematic_StringToInt(ObsCut.Variable);
      ObsCut.LowerBound = CutNode["LowerBound"].as<T>();
      ObsCut.UpperBound = CutNode["UpperBound"].as<T>();
      
      Cuts.push_back(ObsCut);
    }

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

  if (Cuts.size() > 0) {
    std::cout << "\tCuts applied:" << std::endl;
    for (size_t iCut=0;iCut<Cuts.size();iCut++) {
      std::cout << "\t\tVariable: " << Cuts[iCut].Variable << "(" << Cuts[iCut].Variable_Int << ")" << std::endl;
      std::cout << "\t\t\tBound: [" << Cuts[iCut].LowerBound << ", " << Cuts[iCut].UpperBound << ")" << std::endl;
    }
  }

  if (nDimensions == 0) {
    std::cerr << "Invalid number of axes defined for observable:" << Name << std::endl;
    throw;
  } else if (nDimensions == 1) {
    std::string HistName = Name+"_0";
    //std::string HistTitle = Axes[0].Label;
    std::string HistTitle = Name;
    /*for (size_t iCut=0;iCut<Cuts.size();iCut++) {
      HistTitle += " ("+std::string(Form("%4.2f",Cuts[iCut].LowerBound))+"<"+Cuts[iCut].Variable+"<"+Form("%4.2f",Cuts[iCut].UpperBound)+")";
    }*/
    HistTitle += ";"+Axes[0].Label+";Number of Events";
    int nBins = Axes[0].Binning.size()-1;
    
    if (typeid(T) == typeid(float)) {
      HistTemplate = new TH1F(HistName.c_str(),HistTitle.c_str(),nBins,Axes[0].Binning.data());
    } else {
      HistTemplate = new TH1D(HistName.c_str(),HistTitle.c_str(),nBins,Axes[0].Binning.data());
    }
    for (size_t iBin = 0; iBin < Axes[0].BinLabels.size(); iBin++){
      HistTemplate->GetXaxis()->SetBinLabel(iBin+1, Axes[0].BinLabels[iBin].c_str());
    }
    if (Axes[0].BinLabels.size() > 0){ 
      HistTemplate->GetXaxis()->CenterLabels(true);
      HistTemplate->GetXaxis()->SetLabelSize(0.05);
    }
  } else if (nDimensions == 2) {
    std::string HistName = Name+"_0";

    std::string HistTitle = Name+";"+Axes[0].Label;
    for (size_t iCut=0;iCut<Cuts.size();iCut++) {
      HistTitle += " ("+std::string(Form("%4.2f",Cuts[iCut].LowerBound))+"<"+Cuts[iCut].Variable+"<"+Form("%4.2f",Cuts[iCut].UpperBound)+")";
    }
    HistTitle += ";"+Axes[1].Label+";Number of Events";


    int nXBins = Axes[0].Binning.size()-1;
    int nYBins = Axes[1].Binning.size()-1;
    
    if (typeid(T) == typeid(float)) {
      HistTemplate = new TH2F(HistName.c_str(),HistTitle.c_str(),nXBins,Axes[0].Binning.data(),nYBins,Axes[1].Binning.data());
    } else {
      HistTemplate = new TH2D(HistName.c_str(),HistTitle.c_str(),nXBins,Axes[0].Binning.data(),nYBins,Axes[1].Binning.data());
    }
    for (size_t iBin = 0; iBin < Axes[0].BinLabels.size(); iBin++){
      HistTemplate->GetXaxis()->SetBinLabel(iBin+1, Axes[0].BinLabels[iBin].c_str());
    }
    for (size_t iBin = 0; iBin < Axes[1].BinLabels.size(); iBin++){
      HistTemplate->GetYaxis()->SetBinLabel(iBin+1, Axes[1].BinLabels[iBin].c_str());
    }
    if (Axes[0].BinLabels.size() > 0){ 
      HistTemplate->GetXaxis()->CenterLabels(true);
      HistTemplate->GetXaxis()->SetLabelSize(0.05);
    }
    if (Axes[1].BinLabels.size() > 0){ 
      HistTemplate->GetYaxis()->CenterLabels(true);
      HistTemplate->GetYaxis()->SetLabelSize(0.05);
    }
  } else {
    std::cerr << "Currently only have support for 1 and 2 dimension observables" << std::endl;
    throw;
  }

  HistTemplate->SetDirectory(0);

  if (nDimensions == 1) {

  std::string DiffHistName =
    Name + "_Difference";

  std::string DiffHistTitle =
    Name + " Difference;"
    + Axes[0].Label
    + " Difference;Matched Events";

  int nBins = Axes[0].Binning.size() - 1;

  T maxAbs = std::max(
    std::abs(Axes[0].Binning.front()),
    std::abs(Axes[0].Binning.back())
  );

  if (typeid(T) == typeid(float)) {

    DifferenceHistTemplate = new TH1F(
      DiffHistName.c_str(),
      DiffHistTitle.c_str(),
      nBins,
      -maxAbs,
      maxAbs
    );

  } else {

    DifferenceHistTemplate = new TH1D(
      DiffHistName.c_str(),
      DiffHistTitle.c_str(),
      nBins,
      -maxAbs,
      maxAbs
    );
  }

  DifferenceHistTemplate->SetDirectory(0);
  }


  std::cout << std::endl;
}

template<typename T>
ObservableManager<T>::ObservableManager(YAML::Node Config) {
  for (auto ObservableNode: Config["Observables"]) {
    Observables.emplace_back(Observable<T>(ObservableNode));
  }
}
