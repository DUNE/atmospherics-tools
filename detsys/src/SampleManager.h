#pragma once
#include "yaml-cpp/yaml.h"
#include "TH1.h"

#include <vector>

#include "Reader.h"
#include "ObservableManager.h"

template <typename T>
struct Measurement {
  TH1* Histogram;
  int nDimensions;
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

 public:
  Sample(YAML::Node Config);
  void Scale(T ScaleFactor);
  void SetObservables(ObservableManager<T>* Observable);
  void ReadData();
  int GetNEvents() {return SampleReader->GetNentries();}
  std::string GetName() {return Name;}
  TH1* GetMeasurement(int iMeas);
};

template <typename T>
class SampleManager {
 private:
  std::vector<Sample<T>*> Samples;
  ObservableManager<T>* ObsManager;
  
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

  void Plot1D(YAML::Node Config);
};

template class Sample<float>;
template class Sample<double>;
template class SampleManager<float>;
template class SampleManager<double>;
