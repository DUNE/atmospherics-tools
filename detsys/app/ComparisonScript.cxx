#include "Reader.h"

#include "yaml-cpp/yaml.h"

#include "SampleManager.h"
#include "ObservableManager.h"

using FLOAT_T = float;

int main(int argc, char const *argv[]) {

  std::string ConfigFileName = "Config.yaml";
  YAML::Node Config = YAML::LoadFile(ConfigFileName);

  ObservableManager<FLOAT_T> Observables = ObservableManager<FLOAT_T>(Config);

  SampleManager<FLOAT_T> Samples = SampleManager<FLOAT_T>(Config);
  Samples.SetObservables(&Observables);
  Samples.ReadData();

  if (Config["General"]["ScaleNormalisationTo"]) {
    Samples.ScaleToNormalisation(Config["General"]["ScaleNormalisationTo"].as<std::string>());
  }

  YAML::Node DrawOptions = Config["DrawOptions"];
  Samples.Plot1D(DrawOptions);

  return 0;
}
