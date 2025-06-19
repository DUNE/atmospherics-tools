#include "SampleManager.h"

#include "TCanvas.h"
#include "TLegend.h"

template<typename T>
Sample<T>::Sample(YAML::Node SampleConfig) {
  Name = SampleConfig["Name"].as<std::string>();
  SampleColour = SampleConfig["Colour"].as<int>();
  FilePath = SampleConfig["FilePath"].as<std::string>();
  TupleName = SampleConfig["SubDirectoryName"].as<std::string>();

  std::cout << "Initialising sample:" << Name << std::endl;
  SampleReader = new Reader<T>(FilePath, TupleName);

  std::cout << std::endl;
}

template<typename T>
void Sample<T>::ReadData() {
  std::cout << "Reading data from Sample:" << Name << std::endl;

  for (int i=0;i<SampleReader->GetNentries();i++) {
    SampleReader->GetEntry(i);
    
    for (size_t iMeas=0;iMeas<Measurements.size();iMeas++) {
      T EventWeight = SampleReader->GetEventWeight();

      if (Measurements[iMeas].nDimensions == 1) {
	Measurements[iMeas].Histogram->Fill(SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[0]),EventWeight);
      } else if (Measurements[iMeas].nDimensions == 2) {
	T XVar = SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[0]);
	T YVar = SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[1]);

	int nBin = Measurements[iMeas].Histogram->FindBin(XVar,YVar);
	Measurements[iMeas].Histogram->Fill(nBin,EventWeight);
      } else {
	std::cerr << "Invalid number of axes! >2" << std::endl;
	throw;
      }

    }
  }
}

template<typename T>
void Sample<T>::SetObservables(ObservableManager<T>* ObsManager) {

  for (int iObs=0;iObs<ObsManager->GetNObservables();iObs++) {
    Measurement<T> Meas = Measurement<T>();

    TH1* Hist = ObsManager->GetObservable(iObs)->ReturnTemplateHistogram();
    Meas.Histogram = static_cast<TH1*>(Hist->Clone());
    Meas.Histogram->Sumw2();

    std::string HistName = std::string(Meas.Histogram->GetName())+"_"+Name;
    std::replace(HistName.begin(),HistName.end(),' ','_');
    Meas.Histogram->SetName(HistName.c_str());

    Meas.nDimensions = ObsManager->GetObservable(iObs)->GetNDimensions();

    std::vector<int> AVariables;
    for (int iAx=0;iAx<ObsManager->GetObservable(iObs)->GetNDimensions();iAx++) {
      AVariables.push_back(ObsManager->GetObservable(iObs)->GetVariable(iAx));
    }
    Meas.AxisVariables = AVariables;

    Measurements.push_back(Meas);
  }

  std::cout << "Sample: " << Name << std::endl;
  std::cout << "\tNumber of Measurements:" << Measurements.size() << std::endl;
  std::cout << std::endl;

  for (size_t iMeas=0;iMeas<Measurements.size();iMeas++) {
    std::cout << "\tMeasurement:" << iMeas << std::endl;

    int nAxes =  Measurements[iMeas].nDimensions;
    int nBins = 0;
    if (nAxes == 0) {
      std::cerr << "Invalid number of axes! = 0" << std::endl;
      throw;
    } else if (nAxes == 1) {
      nBins = (Measurements[iMeas].Histogram)->GetNbinsX();
    } else if (nAxes == 2) {
      nBins = (Measurements[iMeas].Histogram)->GetNbinsX() * (Measurements[iMeas].Histogram)->GetNbinsY();
    } else {
      std::cerr << "Invalid number of axes! >2" << std::endl;
      throw;
    }

    std::cout << "\tNumber of Axes:" << nAxes << std::endl;
    std::cout << "\tNumber of Histogram Bins:" << nBins << std::endl;

    std::cout << "\tVariables [per axis]:";
    for (int iAxis=0;iAxis<nAxes;iAxis++) {
      std::cout << Measurements[iMeas].AxisVariables[iAxis] << ", ";
    }
    std::cout << std::endl;

    std::cout << std::endl;
  }

}

template<typename T>
TH1* Sample<T>::GetMeasurement(int iMeas) {
  if (iMeas >= 0 && iMeas < (int)Measurements.size()) {
    TH1* Hist = Measurements[iMeas].Histogram;
    Hist->SetLineColor(SampleColour);
    return Hist;
  }

  std::cerr << "Invalid index provided:" << iMeas << std::endl;
  std::cerr << "Measurements.size():" << Measurements.size() << std::endl;
  throw;
}

template<typename T>
void Sample<T>::Scale(T ScaleFactor) {
  for (int iMeas=0;iMeas<(int)Measurements.size();iMeas++) {
    Measurements[iMeas].Histogram->Scale(ScaleFactor);
  }
}

//==========================================================================================================================================================

template<typename T>
SampleManager<T>::SampleManager(YAML::Node Config) {
  for (auto SampleNode: Config["Samples"]) {
    Samples.emplace_back(new Sample<T>(SampleNode));
  }
}

template<typename T>
void SampleManager<T>::ScaleToNormalisation(std::string SampleNameToNormTo) {
  int IndexToNormTo = -1;

  for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
    if (Samples[iSamp]->GetName() == SampleNameToNormTo) {
      IndexToNormTo = iSamp;
      break;
    }
  }
  if (IndexToNormTo == -1) {
    std::cerr << "Did not find:" << SampleNameToNormTo << std::endl;
    std::cerr << "Available:" << std::endl;
    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      std::cerr << "\t" << iSamp << " " << Samples[iSamp]->GetName() << std::endl;
    }
  }
  
  for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
    T Factor = (T)Samples[IndexToNormTo]->GetNEvents()/(T)Samples[iSamp]->GetNEvents();
    Samples[iSamp]->Scale(Factor);
  }
}

template<typename T>
void SampleManager<T>::Plot1D(YAML::Node Config) {

  std::string OutputFileName = Config["OutputName"].as<std::string>();
  std::string DrawOptions = Config["DrawOpts"].as<std::string>();
  std::string SampleNameToNormTo = Config["RatioDenominatorSample"].as<std::string>();

  T LegendHeight = 0.2;
  if (Config["LegendHeight"]) {
    LegendHeight = Config["LegendHeight"].as<T>();
  }

  T FontSize = _BAD_VALUE_;
  if (Config["LegendFontSize"]) {
    FontSize = Config["LegendFontSize"].as<T>();
  }

  T RatioYAxisMax = _BAD_VALUE_;
  if (Config["RatioMaximum"]) {
    RatioYAxisMax = Config["RatioMaximum"].as<T>();
  }
  T RatioYAxisMin = _BAD_VALUE_;
  if (Config["RatioMinimum"]) {
    RatioYAxisMin = Config["RatioMinimum"].as<T>();
  }

  int nSamples = Samples.size();
  int nObservables = ObsManager->GetNObservables();

  TCanvas* Canv = new TCanvas;
  Canv->SetRightMargin(0.2);
  Canv->Print((OutputFileName+"[").c_str());

  std::vector<TLegend*> Legends(nSamples);
  std::vector<TH1*> Hists(nSamples);

  for (int iObs=0;iObs<nObservables;iObs++) {
    if (ObsManager->GetObservable(iObs)->GetNDimensions() != 1) continue;

    for (int iSamp=0;iSamp<nSamples;iSamp++) {
      Hists[iSamp] = Samples[iSamp]->GetMeasurement(iObs);

      TH1* Hist = Hists[iSamp];
      Hist->SetStats(false);
      
      if (iSamp==0) {
	Hist->Draw(DrawOptions.c_str());
      } else {
	Hist->Draw((DrawOptions+" SAME").c_str());
      }

      Legends[iSamp] = new TLegend(0.8,0.9-(1+iSamp)*LegendHeight,0.99,0.9-iSamp*LegendHeight);
      Legends[iSamp]->SetTextSize(FontSize);
      Legends[iSamp]->AddEntry(Hist,(Samples[iSamp]->GetName()).c_str(),"l");
      Legends[iSamp]->AddEntry((TObject*)0,Form("Mean: %4.5f",Hist->GetMean()),"");
      Legends[iSamp]->AddEntry((TObject*)0,Form("RMS: %4.5f",Hist->GetRMS()),"");
      Legends[iSamp]->Draw();
    }

    Canv->Print(OutputFileName.c_str());

    int IndexToNormTo = -1;
    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      if (Samples[iSamp]->GetName() == SampleNameToNormTo) {
	IndexToNormTo = iSamp;
	break;
      }
    }
    if (IndexToNormTo == -1) {
      std::cerr << "Did not find:" << SampleNameToNormTo << std::endl;
      std::cerr << "Available:" << std::endl;
      for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
	std::cerr << "\t" << iSamp << " " << Samples[iSamp]->GetName() << std::endl;
      }
    }
    
    for (int iSamp=0;iSamp<nSamples;iSamp++) {
      if (iSamp == IndexToNormTo) continue;
      Hists[iSamp]->Divide(Hists[IndexToNormTo]);
    }
    Hists[IndexToNormTo]->Divide(Hists[IndexToNormTo]);

    if (RatioYAxisMax != _BAD_VALUE_) {
      for (int iSamp=0;iSamp<nSamples;iSamp++) {
	Hists[iSamp]->SetMaximum(RatioYAxisMax);
	Hists[iSamp]->SetMinimum(RatioYAxisMin);
      }
    } else {
      T Max = -1e8;
      T Min = 1e8;

      for (int iSamp=0;iSamp<nSamples;iSamp++) {
	if (Hists[iSamp]->GetMaximum() > Max) {Max = Hists[iSamp]->GetMaximum();}
	if (Hists[iSamp]->GetMinimum() < Min) {Min = Hists[iSamp]->GetMinimum();}
      }

      for (int iSamp=0;iSamp<nSamples;iSamp++) {
	Hists[iSamp]->SetMaximum(Max+1.3*(Max-Min));
	Hists[iSamp]->SetMinimum(Min-1.3*(Max-Min));
      }
    }

    for (int iSamp=0;iSamp<nSamples;iSamp++) {
      Hists[iSamp]->GetYaxis()->SetTitle(std::string("Ratio to "+SampleNameToNormTo).c_str());

      if (iSamp==0) {
        Hists[iSamp]->Draw(DrawOptions.c_str());
      } else {
        Hists[iSamp]->Draw((DrawOptions+" SAME").c_str());
      }

      Legends[iSamp] = new TLegend(0.8,0.9-(1+iSamp)*LegendHeight,0.99,0.9-iSamp*LegendHeight);
      Legends[iSamp]->SetTextSize(FontSize);
      Legends[iSamp]->AddEntry(Hists[iSamp],(Samples[iSamp]->GetName()).c_str(),"l");
      Legends[iSamp]->AddEntry((TObject*)0,Form("Mean: %4.5f",Hists[iSamp]->GetMean()),"");
      Legends[iSamp]->AddEntry((TObject*)0,Form("RMS: %4.5f",Hists[iSamp]->GetRMS()),"");
      Legends[iSamp]->Draw();
    }

    Canv->Print(OutputFileName.c_str());    
  }

  Canv->Print((OutputFileName+"]").c_str());
}
