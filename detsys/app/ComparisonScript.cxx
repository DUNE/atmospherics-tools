#include "Reader.h"

#include "yaml-cpp/yaml.h"

#include "SampleManager.h"
#include "ObservableManager.h"
#include "AnalysisBinningManager.h"

using FLOAT_T = float;

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: ./ComparisonScript Config.yaml" << std::endl;
    throw; 
  }
  
  std::string ConfigFileName = argv[1];
  YAML::Node Config = YAML::LoadFile(ConfigFileName);

  AnalysisBinningManager<FLOAT_T> AnalysisBinning = AnalysisBinningManager<FLOAT_T>(Config);

  ObservableManager<FLOAT_T> Observables = ObservableManager<FLOAT_T>(Config);

  SampleManager<FLOAT_T> Samples = SampleManager<FLOAT_T>(Config);
  Samples.SetAnalysisBinning(&AnalysisBinning);
  Samples.SetObservables(&Observables);
  Samples.ReadData();
  std::cout << std::endl;

  if (Config["General"]["ScaleNormalisationTo"]) {
    Samples.ScaleToNormalisation(Config["General"]["ScaleNormalisationTo"].as<std::string>());
  }

  if (Config["DrawOptions"]["OneDimension"]) {
    YAML::Node OneDimensionDrawingConfig = Config["DrawOptions"]["OneDimension"];
    Samples.Plot1D(OneDimensionDrawingConfig);
  } else {
    std::cout << "Did not find the [DrawOptions][OneDimension] node in the provided config - not drawing 1D distributions" << std::endl;
  }

  if (Config["DrawOptions"]["AnalysisBinning"]) {
    YAML::Node AnalysisBinningDrawingConfig = Config["DrawOptions"]["AnalysisBinning"];
    Samples.PlotAnalysisBinning(AnalysisBinningDrawingConfig);
  } else {
    std::cout << "Did not find the [DrawOptions][OneDimension] node in the provided config - not drawing AnalysisBinning distributions" << std::endl;
  }

  if (Config["DrawOptions"]["TwoDimension"]) {
    YAML::Node TwoDimensionDrawingConfig = Config["DrawOptions"]["TwoDimension"];
    Samples.Plot2D(TwoDimensionDrawingConfig);
  } else {
    std::cout << "Did not find the [DrawOptions][TwoDimension] node in the provided config - not drawing 2D distributions" << std::endl;
  }

  return 0;
}
