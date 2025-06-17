#include "SampleManager.h"

#include "TCanvas.h"

template<typename T>
Sample<T>::Sample(YAML::Node SampleConfig) {
  Name = SampleConfig["Name"].as<std::string>();
  SampleColour = SampleConfig["Colour"].as<int>();
  FilePath = SampleConfig["FilePath"].as<std::string>();
  TupleName = SampleConfig["TupleName"].as<std::string>();

  std::cout << "Initialising sample:" << Name << std::endl;

  SampleReader = new Reader<T>(FilePath, TupleName);
  std::cout << "\tNEntries:" << SampleReader->GetNentries() << std::endl;

  std::cout << std::endl;
}

template<typename T>
void Sample<T>::ReadData() {
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

  std::cout << std::endl;
}

template<typename T>
void Sample<T>::SetObservables(ObservableManager<T>* ObsManager) {

  for (int iObs=0;iObs<ObsManager->GetNObservables();iObs++) {
    Measurement<T> Meas = Measurement<T>();

    TH1* Hist = ObsManager->GetObservable(iObs)->ReturnTemplateHistogram();
    Meas.Histogram = static_cast<TH1*>(Hist->Clone());

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

//==========================================================================================================================================================

template<typename T>
SampleManager<T>::SampleManager(YAML::Node Config) {
  for (auto SampleNode: Config["Samples"]) {
    Samples.emplace_back(new Sample<T>(SampleNode));
  }
}

template<typename T>
void SampleManager<T>::Plot1D(std::string OutputFileName) {
  int nSamples = Samples.size();
  int nObservables = ObsManager->GetNObservables();

  TCanvas* Canv = new TCanvas;
  Canv->Print((OutputFileName+"[").c_str());

  for (int iObs=0;iObs<nObservables;iObs++) {
    if (ObsManager->GetObservable(iObs)->GetNDimensions() != 1) continue;

    for (int iSamp=0;iSamp<nSamples;iSamp++) {
      TH1* Hist = Samples[iSamp]->GetMeasurement(iObs);
      
      if (iSamp==0) {
	Hist->Draw();
      } else {
	Hist->Draw("SAME");
      }
    }

    Canv->Print(OutputFileName.c_str());
  }

  Canv->Print((OutputFileName+"]").c_str());
}
