#include "AnalysisBinningManager.h"
#include <iostream>

template<typename T>
AnalysisBinningManager<T>::AnalysisBinningManager(std::string FilePath_) : AnalysisBinningManager(YAML::LoadFile(FilePath_)) {
}

template<typename T>
AnalysisBinningManager<T>::AnalysisBinningManager(YAML::Node Config_) {
  Config = Config_;

  if (!Config["AnalysisBinning"]) {
    std::cerr << "Config provided to AnalysisBinningManager does not contain the AnalysisBinning node" << std::endl;
    throw;
  }

  std::cout << "\nAnalysis Binning used for calculation of covariance/correlation matrices -" << std::endl;
  for (const auto& Selection : Config["AnalysisBinning"]) {

    int SelectionIndex = -1;
    std::string SelectionName = "";
    std::vector<std::string> BinVars;
    std::vector< std::vector<T> > BinEdges;
    
    for (const auto& Att : Selection) {
      std::string Key = Att.first.as<std::string>();

      if (Key == "Selection") {
	SelectionName = Att.second.as<std::string>();
      } else if (Key == "kSelection") {
        SelectionIndex = Att.second.as<int>();
      } else {
	BinVars.push_back(Key);
	std::vector<T> Vec = Att.second.as< std::vector<T> >();
	BinEdges.push_back(Vec);
      }
    }
    if (BinVars.size() != BinEdges.size()) {
      std::cerr << "Size of BinVars and BinEdges is different! Problem" << std::endl;
      throw;
    }

    SelectionBinning<T> SelecBinning;
    SelecBinning.SelectionIndex = SelectionIndex;
    SelecBinning.SelectionName = SelectionName;
    SelecBinning.BinVars = BinVars;
    SelecBinning.BinEdges = BinEdges;
   
    AnalysisSelectionBinning.push_back(SelecBinning);
  }

  for (size_t iSelecBinning=0;iSelecBinning<AnalysisSelectionBinning.size();iSelecBinning++) {
    std::cout << "\tSelectionIndex:" << AnalysisSelectionBinning[iSelecBinning].SelectionIndex << std::endl;
    std::cout << "\tSelectionName:" << AnalysisSelectionBinning[iSelecBinning].SelectionName << std::endl;
    for (size_t iVar=0;iVar<AnalysisSelectionBinning[iSelecBinning].BinVars.size();iVar++) {
      std::cout << "\t" << AnalysisSelectionBinning[iSelecBinning].BinVars[iVar] << ":" << std::endl;
      std::cout << "\t\t";
      for (size_t iBinEdge=0;iBinEdge<AnalysisSelectionBinning[iSelecBinning].BinEdges[iVar].size();iBinEdge++) {
	std::cout << AnalysisSelectionBinning[iSelecBinning].BinEdges[iVar][iBinEdge] << ", ";
      }
      std::cout << std::endl;
    }
  }
  std::cout << "\n" << std::endl;

}

template<typename T>
bool AnalysisBinningManager<T>::CheckSelectionInAnalysisBinning(int SelectionIndex) {
  for (size_t iSelecBinning=0;iSelecBinning<AnalysisSelectionBinning.size();iSelecBinning++) {
    if (AnalysisSelectionBinning[iSelecBinning].SelectionIndex == SelectionIndex) {
      return true;
    }
  }
  return false;
}

template<typename T>
std::vector<std::string> AnalysisBinningManager<T>::GetSelectionBinVars(int SelectionIndex) {
  for (size_t iSelecBinning=0;iSelecBinning<AnalysisSelectionBinning.size();iSelecBinning++) {
    if (AnalysisSelectionBinning[iSelecBinning].SelectionIndex == SelectionIndex) {
      return AnalysisSelectionBinning[iSelecBinning].BinVars;
    }
  }

  std::cerr << "Did not find Selection:" << SelectionIndex << " defined in the AnalysisBinning!" << std::endl;
  throw;
}

template<typename T>
int AnalysisBinningManager<T>::GetNBinsFromSelection(int SelectionIndex) {
  if (SelectionIndex < 0 || SelectionIndex>=static_cast<int>(AnalysisSelectionBinning.size())) {
    std::cerr << "Invalid SelectionIndex:" << SelectionIndex << std::endl;
    throw;
  }

  int nBins_Selection = 1;
  for (size_t iBinVar=0;iBinVar<AnalysisSelectionBinning[SelectionIndex].BinVars.size();iBinVar++) {
    nBins_Selection *= static_cast<int>((AnalysisSelectionBinning[SelectionIndex].BinEdges[iBinVar].size()-1));
  }

  return nBins_Selection;
}

template<typename T>
int AnalysisBinningManager<T>::GetNBins() {
  int nBins = 0;
  for (size_t iSelecBinning=0;iSelecBinning<AnalysisSelectionBinning.size();iSelecBinning++) {
    nBins += GetNBinsFromSelection(iSelecBinning);
  }

  return nBins;
}

template<typename T>
int AnalysisBinningManager<T>::GetBin(int SelectionIndex, std::vector<T> EventDetails) {
  int AnalysisSelectionIndex = -1;

  //First check which index selection we are working with
  for (size_t iSelecBinning=0;iSelecBinning<AnalysisSelectionBinning.size();iSelecBinning++) {
    if (AnalysisSelectionBinning[iSelecBinning].SelectionIndex == SelectionIndex) {
      AnalysisSelectionIndex = static_cast<int>(iSelecBinning);
    }
  }
  if (AnalysisSelectionIndex == -1) {
    //std::cerr << "Did not find SelectionIndex: " << SelectionIndex << " in the AnalysisBinning object" << std::endl;
    //throw;
    return -1;
  }

  //Second check we have the correct number of EventDetails for the given selection
  if (EventDetails.size() != AnalysisSelectionBinning[AnalysisSelectionIndex].BinVars.size()) {
    std::cerr << "Number of Event details provided for event does not equal the dimension of the selection the event was identified with" << std::endl;
    std::cerr << "EventDetails.size():" << EventDetails.size() << std::endl;
    std::cerr << "AnalysisSelectionBinning[AnalysisSelectionIndex].BinVars.size():" << AnalysisSelectionBinning[AnalysisSelectionIndex].BinVars.size() << std::endl;
    throw;
  }

  //Third now calculate the bin indexs in each BinVar
  std::vector<int> BinIndices(EventDetails.size(),-1);
  for (size_t iBinVar=0;iBinVar<AnalysisSelectionBinning[AnalysisSelectionIndex].BinVars.size();iBinVar++) {
    T EventAtt = EventDetails[iBinVar];
    
    for (size_t iBinEdge=0;iBinEdge<AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[iBinVar].size()-1;iBinEdge++) {
      T LowBinEdge = AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[iBinVar][iBinEdge];
      T UppBinEdge = AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[iBinVar][iBinEdge+1];

      if ((EventAtt >= LowBinEdge) && (EventAtt < UppBinEdge)) {
	BinIndices[iBinVar] = iBinEdge;
      }
    }
  }

  //Now check we got a bin for all dimensions
  for (size_t iBinVar=0;iBinVar<BinIndices.size();iBinVar++) {
    if (BinIndices[iBinVar] == -1) {
      std::cerr << "Did not find BinIndex for event in dimension:" << AnalysisSelectionBinning[AnalysisSelectionIndex].BinVars[iBinVar] << std::endl;
      std::cerr << "EventDetails[iBinVar]:" << EventDetails[iBinVar] << std::endl;
      std::cerr << "Analysis binning:\n\t";
      for(size_t iBinEdge=0;iBinEdge<AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[iBinVar].size();iBinEdge++) {
	std::cerr << AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[iBinVar][iBinEdge] << ", ";
      }
      std::cerr << std::endl;
      throw;
    }
  }

  //At this point, we have all the bin indexs we need - so first count all bins in selections previous to the one identified
  int GlobalBinIndex = 0;
  for (int iSelecBinning=0;iSelecBinning<AnalysisSelectionIndex;iSelecBinning++) {
    GlobalBinIndex += GetNBinsFromSelection(iSelecBinning);
  }

  //Next define 1st BinVar major axis, and 2nd BinVar minor axis
  int LocalBinIndex = 0;
  for (size_t iBinVar=0;iBinVar<BinIndices.size();iBinVar++) {
    int MajorIndex = BinIndices[iBinVar];
    int nBinsInMinorIndices = 1;
    for (size_t jBinVar=(iBinVar+1);jBinVar<BinIndices.size();jBinVar++) {
      nBinsInMinorIndices *= (AnalysisSelectionBinning[AnalysisSelectionIndex].BinEdges[jBinVar].size()-1);
    }
    LocalBinIndex += MajorIndex*nBinsInMinorIndices;
  }
  GlobalBinIndex += LocalBinIndex;

  if (GlobalBinIndex < 0 || GlobalBinIndex >= GetNBins()) {
    std::cerr << "Invalid Global Bin Index found!" << std::endl;
    throw;
  }

  return GlobalBinIndex;
}
