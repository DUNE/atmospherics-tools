#pragma once
#include "yaml-cpp/yaml.h"
#include "TH1.h"
#include "TCanvas.h"

#include <vector>

#include "Reader.h"
#include "ObservableManager.h"
#include "AnalysisBinningManager.h"
#include "FluxManager.h"
#include "Oscillator/OscillatorFactory.h"

template <typename T>
struct Measurement {
  TH1* Histogram;
  int nDimensions;
  std::vector< Cut<T> > Cuts;
  std::vector<int> AxisVariables;
};

template <typename T>
class Sample {
 private:
  std::string Name;
  int SampleColour;
  std::string FilePath;
  std::string TupleName;

  Reader<T>* SampleReader;
  std::vector<Measurement<T>> Measurements;
  FluxManager* FlxMgr;
  AnalysisBinningManager<T>* AnalysisBinning;
  TH1D* AnalysisBinningHistogram;
  OscillatorBase* OscillBase;

 public:
  Sample(YAML::Node Config);
  void Scale(T ScaleFactor);
  void SetObservables(ObservableManager<T>* Observable);
  void SetAnalysisBinning(AnalysisBinningManager<T>* AnalysisBinning_);
  void ReadData();
  int GetNEvents() {return SampleReader->GetNentries();}
  std::string GetName() {return Name;}
  TH1* GetMeasurement(int iMeas);
  TH1* GetAnalysisBinningHistogram() {return AnalysisBinningHistogram;}
  void SetFluxManager(FluxManager* FlxMgr_) {
    FlxMgr = FlxMgr_;
    SampleReader->SetFluxManager(FlxMgr);
  }
  void SetOscillator(OscillatorBase* OscillBase_) {
    OscillBase = OscillBase_;
    SampleReader->SetOscillator(OscillBase);
  }
};

template <typename T>
class SampleManager {
 private:
  std::vector<Sample<T>*> Samples;
  ObservableManager<T>* ObsManager;
  AnalysisBinningManager<T>* AnalysisBinning;
  FluxManager* FlxMgr;
  OscillatorBase* OscillBase;

  std::string OutputFileName_1D;
  std::string OutputFileName_2D;
  std::string DrawOptions_1D;
  std::string DrawOptions_2D;
  std::string SampleNameToRatioTo;
  int IndexToRatioTo;
  T LegendHeight;
  T FontSize;
  T RatioYAxisMax;
  T RatioYAxisMin;

  void Plot1DRatioHists(TCanvas* Canv, std::vector<TH1*> Hists);
  void Plot1DHists(TCanvas* Canv, std::vector<TH1*> Hists);
  T CalculateCovariance(T XNominalBinContent, std::vector<T> XVariedBinContents, T YNominalBinContent, std::vector<T> YVariedBinContents);

 public:
  SampleManager(YAML::Node Config);
  void ScaleToNormalisation(std::string SampleNameToNormTo);

  void SetObservables(ObservableManager<T>* Observable_) {
    ObsManager = Observable_;

    for (Sample<T>* Samp : Samples) {
      Samp->SetObservables(Observable_);
    }
  }

  void ReadData() {
    for (Sample<T>* Samp : Samples) {
      Samp->ReadData();
    }
  }

  void PlotAnalysisBinning(YAML::Node Config);
  void Plot1D(YAML::Node Config);
  void Plot2D(YAML::Node Config);

  void SetAnalysisBinning(AnalysisBinningManager<T>* AnalysisBinning_) {
    AnalysisBinning = AnalysisBinning_;

    for (Sample<T>* Samp : Samples) {
      Samp->SetAnalysisBinning(AnalysisBinning_);
    }
  }

  void SetFluxManager(FluxManager* FlxMgr_) {
    FlxMgr = FlxMgr_;
    for (Sample<T>* Samp : Samples) {
      Samp->SetFluxManager(FlxMgr);
    }
  }

  void SetOscillator(OscillatorBase* OscillBase_) {
    OscillBase = OscillBase_;
    for (Sample<T>* Samp : Samples) {
      Samp->SetOscillator(OscillBase);
    }
  }

};

template class Sample<float>;
template class Sample<double>;
template class SampleManager<float>;
template class SampleManager<double>;
