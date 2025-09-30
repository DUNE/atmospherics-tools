#pragma once

#include "yaml-cpp/yaml.h"
#include "TH1.h"
#include "TH2.h"

template <typename T>
struct Cut {
  std::string Variable;
  int Variable_Int;
  T LowerBound;
  T UpperBound;
};

template <typename T>
struct Axis {
  std::string Variable;
  std::string Label;
  int Variable_Int;
  std::vector<T> Binning;
  bool isLog;
};

template <typename T>
class Observable {
 private:
  std::string Name;
  int nDimensions; 
  TH1* HistTemplate;
  std::vector< Cut<T> > Cuts;
  
  std::vector<Axis<T>> Axes;
 public:
  Observable(YAML::Node Config);

  TH1* ReturnTemplateHistogram() {return HistTemplate;}
  std::vector< Cut<T> > GetCuts() {return Cuts;}
  int GetVariable(int iAxis) {return (Axes.at(iAxis)).Variable_Int;}
  bool GetIsLog(int iAxis) {return (Axes.at(iAxis)).isLog;}
  int GetNDimensions() {return nDimensions;}
};

template <typename T>
class ObservableManager {
 private:
  std::vector<Observable<T>> Observables;
  
 public:
  ObservableManager(YAML::Node Config);

  Observable<T>* GetObservable(int i) {
    return &(Observables.at(i));
  }

  int GetNObservables() {return Observables.size();}
};

template class Observable<float>;
template class Observable<double>;
template class ObservableManager<float>;
template class ObservableManager<double>;
