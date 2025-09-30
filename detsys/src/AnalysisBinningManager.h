#pragma once
#include "yaml-cpp/yaml.h"

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
 public:
  AnalysisBinningManager(std::string FilePath_);
  AnalysisBinningManager(YAML::Node Config_);

  int GetNBinsFromSelection(int SelectionIndex);
  int GetNBins();

  bool CheckSelectionInAnalysisBinning(int SelectionIndex);
  std::vector<std::string> GetSelectionBinVars(int Selection);
  int GetBin(int SelectionIndex, std::vector<T> EventDetails);
};

template struct SelectionBinning<float>;
template struct SelectionBinning<double>;
template class AnalysisBinningManager<float>;
template class AnalysisBinningManager<double>;
