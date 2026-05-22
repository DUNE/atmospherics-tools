#include "Reader.h"
#include "SampleManager.h"
#include <yaml-cpp/yaml.h>
#include "Oscillator/OscillatorFactory.h"
#include "ObservableManager.h"
#include "AnalysisBinningManager.h"
#include "FluxManager.h"
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <functional>

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"

template<typename T>
struct EventID {

  int run;
  int subrun;
  int event;
  /*int nuPDG;
  int mode;

  int px;
  int py;
  int pz;

  int vx;
  int vy;
  int vz;*/

  bool operator<(const EventID<T>& other) const {

    return std::tie(
      run,
      subrun,
      event
      /*nuPDG,
      mode,

      px,
      py,
      pz,

      vx,
      vy,
      vz*/
    )
    <
    std::tie(
      other.run,
      other.subrun,
      other.event
      /*other.nuPDG,
      other.mode,

      other.px,
      other.py,
      other.pz,

      other.vx,
      other.vy,
      other.vz*/
    );
  }
};



struct QuantityConfig {

  std::string name;

  int bins;
  double min;
  double max;
};


template<typename T>
struct EventData {

  int run;
  int subrun;
  int event;
  int nuPDG;
  int mode;

  T NuMomX;
  T NuMomY;
  T NuMomZ;

  T vtx_x;
  T vtx_y;
  T vtx_z;

  T erec;
  T RecoCZ;
  T had_erec;
  T lep_erec;

  T ev;

  T cvn_numu;
  T cvn_nue;

  int Selection;
  int trueInt;
};


template<typename T>
std::map<EventID<T>, EventData<T>>
BuildEventMap(Reader<T>& reader)
{
  std::map<EventID<T>, EventData<T>> eventMap;

  int nEntries = reader.GetNentries();

  std::cout << "Reading "
            << nEntries
            << " entries"
            << std::endl;

  for (int i = 0; i < nEntries; ++i) {

    if (!reader.GetEntry(i)) {
      continue;
    }

    const auto& data = reader.GetData();

    EventID<T> id {

      data.run,
      data.subrun,
      data.event
      /*data.nuPDG,
      data.mode,

      static_cast<int>(1e6 * data.NuMomX),
      static_cast<int>(1e6 * data.NuMomY),
      static_cast<int>(1e6 * data.NuMomZ),

      static_cast<int>(1e3 * data.vtx_x),
      static_cast<int>(1e3 * data.vtx_y),
      static_cast<int>(1e3 * data.vtx_z)
    */
    };

    EventData<T> ev;

    ev.run      = data.run;
    ev.subrun   = data.subrun;
    ev.event    = data.event;
    ev.nuPDG    = data.nuPDG;
    ev.mode     = data.mode;

    ev.NuMomX   = data.NuMomX;
    ev.NuMomY   = data.NuMomY;
    ev.NuMomZ   = data.NuMomZ;

    ev.vtx_x    = data.vtx_x;
    ev.vtx_y    = data.vtx_y;
    ev.vtx_z    = data.vtx_z;

    ev.erec      = data.erec;
    ev.RecoCZ    = data.RecoCZ;
    ev.had_erec  = data.had_erec;
    ev.lep_erec  = data.lep_erec;

    ev.ev        = data.ev;

    ev.cvn_numu  = data.cvn_numu;
    ev.cvn_nue   = data.cvn_nue;

    ev.Selection = data.Selection;
    ev.trueInt   = data.trueInt;

    eventMap[id] = ev;
  }

  return eventMap;
}

template<typename T>
using QuantityAccessor = std::function<T(const EventData<T>&)>;


template<typename T>
std::map<std::string, QuantityAccessor<T>>
BuildAccessorMap()
{
  std::map<std::string, QuantityAccessor<T>> accessors;

  accessors["erec"] =
    [](const EventData<T>& d){ return d.erec; };

  accessors["RecoCZ"] =
    [](const EventData<T>& d){ return d.RecoCZ; };

  accessors["had_erec"] =
    [](const EventData<T>& d){ return d.had_erec; };

  accessors["lep_erec"] =
    [](const EventData<T>& d){ return d.lep_erec; };

  accessors["cvn_numu"] =
    [](const EventData<T>& d){ return d.cvn_numu; };

  accessors["cvn_nue"] =
    [](const EventData<T>& d){ return d.cvn_nue; };

  accessors["trueE"] =
    [](const EventData<T>& d){ return d.ev; };

  accessors["NuMomX"] =
    [](const EventData<T>& d){ return d.NuMomX; };

  accessors["NuMomY"] =
    [](const EventData<T>& d){ return d.NuMomY; };

  accessors["NuMomZ"] =
    [](const EventData<T>& d){ return d.NuMomZ; };

  accessors["vtx_x"] =
    [](const EventData<T>& d){ return d.vtx_x; };

  accessors["vtx_y"] =
    [](const EventData<T>& d){ return d.vtx_y; };

  accessors["vtx_z"] =
    [](const EventData<T>& d){ return d.vtx_z; };

  return accessors;
}

template<typename T>
void compare_samples_yaml_impl(std::string yamlFile)
{
  YAML::Node config = YAML::LoadFile(yamlFile);

  if (!config["NuOscillator"]) {
    std::cerr << "Did not define the NuOscillator YAML node in the config" << std::endl;
    throw;
  }
  OscillatorFactory OscillFact = OscillatorFactory();
  OscillatorBase* OscillBase = OscillFact.CreateOscillator(config["NuOscillator"]);
  if (!config["NuOscillator"]["OscillationParameters"]) {
    std::cerr << "NuOscillator::General::OscillationParameters node is not defined in the Config!" << std::endl;
    throw;
  }
  std::vector<FLOAT_T> OscillationParameters;
  std::cout << "\nOscillation Parameters used: -" << std::endl;
  for (auto Pair : config["NuOscillator"]["OscillationParameters"]) {
    std::cout << "\t" << std::setw(20) << Pair.first.as<std::string>() << " : " << Pair.second.as<FLOAT_T>() << std::endl;
    OscillationParameters.push_back(Pair.second.as<FLOAT_T>());
  }
  std::cout << "Setting up NuOscillator object..." << std::endl;
  OscillBase->Setup();
  OscillBase->CalculateProbabilities(OscillationParameters);
  std::cout << "\n" << std::endl;



  std::filesystem::path fluxdir(config["Fluxes"]["FluxDir"].as<std::string>());
  std::filesystem::path nue_file(config["Fluxes"]["nue"].as<std::string>());
  std::filesystem::path nuebar_file(config["Fluxes"]["nuebar"].as<std::string>());
  std::filesystem::path numu_file(config["Fluxes"]["numu"].as<std::string>());
  std::filesystem::path numubar_file(config["Fluxes"]["numubar"].as<std::string>());
  std::filesystem::path ref_file(config["Fluxes"]["refflux"].as<std::string>());

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

  FluxManager FlxMgr(fluxes);

  AnalysisBinningManager<T> AnalysisBinning = AnalysisBinningManager<T>(config);
  ObservableManager<T> Observables = ObservableManager<T>(config);



  SampleManager<T> Samples = SampleManager<T>(config);
  Samples.SetAnalysisBinning(&AnalysisBinning);
  Samples.SetObservables(&Observables);
  Samples.SetFluxManager(&FlxMgr);
  Samples.SetOscillator(OscillBase);
  Samples.ReadDataMapping();
  Samples.PlotMatchedDifferences();
/*
  std::string pathA =
    config["sampleA"]["path"].as<std::string>();

  std::string subfolderA =
    config["sampleA"]["subfolder"].as<std::string>();

  std::string pathB =
    config["sampleB"]["path"].as<std::string>();

  std::string subfolderB =
    config["sampleB"]["subfolder"].as<std::string>();

  std::string outputFile =
    config["output"]["file"].as<std::string>();

  std::vector<QuantityConfig> quantities;

  for (const auto& q : config["quantities"]) {

    QuantityConfig qq;

    qq.name = q["name"].as<std::string>();
    qq.bins = q["bins"].as<int>();
    qq.min  = q["min"].as<double>();
    qq.max  = q["max"].as<double>();

    quantities.push_back(qq);
  }

  Reader<T> readerA(pathA, subfolderA);
  Reader<T> readerB(pathB, subfolderB);

  auto mapA = BuildEventMap(readerA);
  auto mapB = BuildEventMap(readerB);

  std::set<EventID<T>> allKeys;

  for (const auto& kv : mapA) {
    allKeys.insert(kv.first);
  }

  for (const auto& kv : mapB) {
    allKeys.insert(kv.first);
  }

  auto accessors = BuildAccessorMap<T>();

  std::map<std::string, TH1D*> hDiff;
  std::map<std::string, TH2D*> hCorr;

  for (const auto& q : quantities) {

    std::string hname = "h_d_" + q.name;

    hDiff[q.name] = new TH1D(
      hname.c_str(),
      (q.name + " A - B;Difference;Events").c_str(),
      q.bins,
      q.min,
      q.max
    );

    std::string h2name = "h2_" + q.name;

    hCorr[q.name] = new TH2D(
      h2name.c_str(),
      (q.name + " A vs B;A;B").c_str(),
      q.bins,
      q.min,
      q.max,
      q.bins,
      q.min,
      q.max
    );
  }

  int commonEvents = 0;
  int missingInA = 0;
  int missingInB = 0;

  for (const auto& key : allKeys) {

    bool inA = mapA.count(key);
    bool inB = mapB.count(key);

    if (!inA) {
      missingInA++;
      continue;
    }

    if (!inB) {
      missingInB++;
      continue;
    }

    commonEvents++;

    const auto& A = mapA.at(key);
    const auto& B = mapB.at(key);

    for (const auto& q : quantities) {

      if (!accessors.count(q.name)) {

        std::cerr << "Unknown quantity: "
                  << q.name
                  << std::endl;

        continue;
      }

      auto accessor = accessors[q.name];

      double valA = accessor(A);
      double valB = accessor(B);

      if (!std::isfinite(valA) ||
          !std::isfinite(valB)) {
        continue;
      }

      double diff = valA - valB;

      hDiff[q.name]->Fill(diff);
      hCorr[q.name]->Fill(valA, valB);
    }
  }

  std::cout << "\n===================================="
            << std::endl;

  std::cout << "Common events : "
            << commonEvents
            << std::endl;

  std::cout << "Missing in A  : "
            << missingInA
            << std::endl;

  std::cout << "Missing in B  : "
            << missingInB
            << std::endl;

  std::cout << "====================================\n"
            << std::endl;

  TFile fout(outputFile.c_str(), "RECREATE");

  for (const auto& kv : hDiff) {
    kv.second->Write();
  }

  for (const auto& kv : hCorr) {
    kv.second->Write();
  }

  fout.Close();

  std::cout << "Wrote "
            << outputFile
            << std::endl;*/
}



int main(int argc, char* argv[])
{
  if (argc != 2) {

    std::cerr << "Usage: "
              << argv[0]
              << " config.yaml"
              << std::endl;

    return 1;
  }

  std::string yamlFile = argv[1];
/*
  YAML::Node Config = YAML::LoadFile(yamlFile);
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

  FluxManager FlxMgr(fluxes);
*/
  compare_samples_yaml_impl<double>(yamlFile);

  return 0;
}
