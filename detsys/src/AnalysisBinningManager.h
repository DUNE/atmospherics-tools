#pragma once
#include "yaml-cpp/yaml.h"
#include <TAxis.h>
#include <string>
#include <vector>


template <typename T>
struct SelectionBinning {
  int SelectionIndex;
  std::string SelectionName;
  std::vector<std::string> BinVars;
  std::vector<std::vector<T>> BinEdges;
};

template <typename T>
class AnalysisBinningManager {
 private:
  YAML::Node Config;
  std::vector<SelectionBinning<T>> AnalysisSelectionBinning;
  TAxis AnalysisBinningAxis;
  std::vector<std::string> ParamNames;

  void BuildAnalysisBinningMetadata();

 public:
  AnalysisBinningManager(std::string FilePath_);
  AnalysisBinningManager(YAML::Node Config_);

  int GetNBinsFromSelection(int SelectionIndex);
  int GetNBins();

  bool CheckSelectionInAnalysisBinning(int SelectionIndex);
  std::vector<std::string> GetSelectionBinVars(int Selection);
  int GetBin(int SelectionIndex, std::vector<T> EventDetails);

  const TAxis& GetAnalysisBinningAxis() const {return AnalysisBinningAxis;}

  const std::vector<std::string>& GetParamNames() const {return ParamNames;}

};

template struct SelectionBinning<float>;
template struct SelectionBinning<double>;
template class AnalysisBinningManager<float>;
template class AnalysisBinningManager<double>;
