#include "Reader.h"

#include <filesystem>
#include "yaml-cpp/yaml.h"

#include "SampleManager.h"
#include "ObservableManager.h"
#include "AnalysisBinningManager.h"
#include "FluxManager.h"
#include "Oscillator/OscillatorFactory.h"

using FLOAT_T = double;

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: ./ComparisonScript Config.yaml" << std::endl;
    throw; 
  }
  
  std::string ConfigFileName = argv[1];
  YAML::Node Config = YAML::LoadFile(ConfigFileName);

  //============================================================================================================================================================
  if (!Config["NuOscillator"]) {
    std::cerr << "Did not define the NuOscillator YAML node in the config" << std::endl;
    throw;
  }
  OscillatorFactory OscillFact = OscillatorFactory();
  OscillatorBase* OscillBase = OscillFact.CreateOscillator(Config["NuOscillator"]);
  if (!Config["NuOscillator"]["OscillationParameters"]) {
    std::cerr << "NuOscillator::General::OscillationParameters node is not defined in the Config!" << std::endl;
    throw;
  }
  std::vector<FLOAT_T> OscillationParameters;
  std::cout << "\nOscillation Parameters used: -" << std::endl;
  for (auto Pair : Config["NuOscillator"]["OscillationParameters"]) {
    std::cout << "\t" << std::setw(20) << Pair.first.as<std::string>() << " : " << Pair.second.as<FLOAT_T>() << std::endl;
    OscillationParameters.push_back(Pair.second.as<FLOAT_T>());
  }
  std::cout << "Setting up NuOscillator object..." << std::endl;
  OscillBase->Setup();
  OscillBase->CalculateProbabilities(OscillationParameters);
  std::cout << "\n" << std::endl;

  //============================================================================================================================================================

  std::filesystem::path fluxdir(Config["Fluxes"]["FluxDir"].as<std::string>());
  std::filesystem::path nue_file(Config["Fluxes"]["nue"].as<std::string>());
  std::filesystem::path nuebar_file(Config["Fluxes"]["nuebar"].as<std::string>());
  std::filesystem::path numu_file(Config["Fluxes"]["numu"].as<std::string>());
  std::filesystem::path numubar_file(Config["Fluxes"]["numubar"].as<std::string>());
  std::filesystem::path ref_file(Config["Fluxes"]["refflux"].as<std::string>());

  nue_file = fluxdir / nue_file;
  nuebar_file = fluxdir / nuebar_file;
  numu_file = fluxdir / numu_file;
  numubar_file = fluxdir / numubar_file;

  std::map<Flavour, std::string> fluxes = {
    {Flavour::NuE, nue_file},
    {Flavour::NuMu, numu_file},
    {Flavour::NuEBar, nuebar_file},
    {Flavour::NuMuBar, numubar_file},
    {Flavour::Reference, ref_file}
  };

  FluxManager manager(fluxes);

  //============================================================================================================================================================
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

  //============================================================================================================================================================
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
