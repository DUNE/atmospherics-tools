#include "SampleManager.h"

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
    
    T EventWeight = SampleReader->GetEventWeight();
    AnalysisBinningHistogram->Fill(SampleReader->ReturnKinematicParameter(kAnalysisBin),EventWeight);

    for (size_t iMeas=0;iMeas<Measurements.size();iMeas++) {
      bool PassesCut = true;

      for (size_t iCut=0;iCut<Measurements[iMeas].Cuts.size();iCut++) {
	int Variable = Measurements[iMeas].Cuts[iCut].Variable_Int;
	T LowerBound = Measurements[iMeas].Cuts[iCut].LowerBound;
	T UpperBound = Measurements[iMeas].Cuts[iCut].UpperBound;
	T VariableValue = SampleReader->ReturnKinematicParameter(Variable);

	if (VariableValue < LowerBound) {
	  PassesCut = false;
	  break;
	}
	if (VariableValue >= UpperBound) {
	  PassesCut = false;
          break;
	}
      }
      if (!PassesCut) {
	continue;
      }

      if (Measurements[iMeas].nDimensions == 1) {
	Measurements[iMeas].Histogram->Fill(SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[0]),EventWeight);
      } else if (Measurements[iMeas].nDimensions == 2) {
	T XVar = SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[0]);
	T YVar = SampleReader->ReturnKinematicParameter(Measurements[iMeas].AxisVariables[1]);

	TH2* Hist2D = dynamic_cast<TH2*>(Measurements[iMeas].Histogram);
	Hist2D->Fill(XVar,YVar,EventWeight);
      } else {
	std::cerr << "Invalid number of axes! >2" << std::endl;
	throw;
      }

    }
  }

  for (size_t iMeas=0;iMeas<Measurements.size();iMeas++) {
    int nAxes =  Measurements[iMeas].nDimensions;
    if (nAxes == 0) {
      std::cerr << "Invalid number of axes! = 0" << std::endl;
      throw;
    } else if (nAxes == 1) {
      Measurements[iMeas].Histogram->GetYaxis()->SetTitle((std::string(Measurements[iMeas].Histogram->GetYaxis()->GetTitle())+"/Bin Width").c_str());
    } else if (nAxes == 2) {
      Measurements[iMeas].Histogram->GetZaxis()->SetTitle((std::string(Measurements[iMeas].Histogram->GetZaxis()->GetTitle())+"/Bin Width").c_str());
    } else {
      std::cerr << "Invalid number of axes! >2" << std::endl;
      throw;
    }

    Measurements[iMeas].Histogram->Scale(1.0,"width");
  }
}

template<typename T>
void Sample<T>::SetAnalysisBinning(AnalysisBinningManager<T>* AnalysisBinning_) {
  SampleReader->SetAnalysisBinning(AnalysisBinning_);
  AnalysisBinning = AnalysisBinning_;

  AnalysisBinningHistogram = new TH1D((Name+"_AnalysisBinning").c_str(),"Analysis Binning;Bin Number;Events",AnalysisBinning->GetNBins(),0,AnalysisBinning->GetNBins());
  AnalysisBinningHistogram->SetLineColor(SampleColour);
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

    Meas.Cuts = ObsManager->GetObservable(iObs)->GetCuts();

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
  AnalysisBinningHistogram->Scale(ScaleFactor);
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
    T Factor = static_cast<T>(Samples[IndexToNormTo]->GetNEvents())/static_cast<T>(Samples[iSamp]->GetNEvents());
    Samples[iSamp]->Scale(Factor);
  }
}

template<typename T>
void SampleManager<T>::Plot1DRatioHists(TCanvas* Canv, std::vector<TH1*> Hists) {

  std::vector<TLegend*> Legends(Hists.size());

  for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
    if (static_cast<int>(iSamp) == IndexToRatioTo) continue;
    Hists[iSamp]->Divide(Hists[IndexToRatioTo]);
  }
  Hists[IndexToRatioTo]->Divide(Hists[IndexToRatioTo]);

  if (RatioYAxisMax != _BAD_VALUE_) {
    for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
      Hists[iSamp]->SetMaximum(RatioYAxisMax);
      Hists[iSamp]->SetMinimum(RatioYAxisMin);
    }
  } else {
    T Max = -1e8;
    T Min = 1e8;

    for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
      for (int xBin=1;xBin<=Hists[iSamp]->GetNbinsX();xBin++) {
	if ((Hists[iSamp]->GetBinContent(xBin)+Hists[iSamp]->GetBinError(xBin)) > Max) {Max = (Hists[iSamp]->GetBinContent(xBin)+Hists[iSamp]->GetBinError(xBin));}
	if ((Hists[iSamp]->GetBinContent(xBin)-Hists[iSamp]->GetBinError(xBin)) < Min) {Min = (Hists[iSamp]->GetBinContent(xBin)-Hists[iSamp]->GetBinError(xBin));}
      }
    }

    for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
      Hists[iSamp]->GetYaxis()->SetRangeUser(Min-0.1*(Max-Min),Max+0.1*(Max-Min));
    }
  }


  for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
    Hists[iSamp]->GetYaxis()->SetTitle(std::string("Ratio to "+SampleNameToRatioTo).c_str());

    if (iSamp==0) {
      Hists[iSamp]->Draw((DrawOptions_1D).c_str());
    } else {
      Hists[iSamp]->Draw((DrawOptions_1D+" SAME").c_str());
    }

    Legends[iSamp] = new TLegend(0.8,0.9-(1.0+static_cast<T>(iSamp))*LegendHeight,0.99,0.9-static_cast<T>(iSamp)*LegendHeight);
    Legends[iSamp]->SetTextSize(FontSize);
    Legends[iSamp]->AddEntry(Hists[iSamp],(Samples[iSamp]->GetName()).c_str(),"l");
    Legends[iSamp]->AddEntry((TObject*)0,Form("Mean: %4.5f",Hists[iSamp]->GetMean()),"");
    Legends[iSamp]->AddEntry((TObject*)0,Form("RMS: %4.5f",Hists[iSamp]->GetRMS()),"");
    Legends[iSamp]->Draw();
  }

  Canv->Print(OutputFileName_1D.c_str());
}

template<typename T>
void SampleManager<T>::Plot1DHists(TCanvas* Canv, std::vector<TH1*> Hists) {

  std::vector<TLegend*> Legends(Hists.size());

  T Max = -1e8;
  T Min = 1e8;
  
  for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
    for (int xBin=1;xBin<=Hists[iSamp]->GetNbinsX();xBin++) {
      if ((Hists[iSamp]->GetBinContent(xBin)+Hists[iSamp]->GetBinError(xBin)) > Max) {Max = (Hists[iSamp]->GetBinContent(xBin)+Hists[iSamp]->GetBinError(xBin));}
      if ((Hists[iSamp]->GetBinContent(xBin)-Hists[iSamp]->GetBinError(xBin)) < Min) {Min = (Hists[iSamp]->GetBinContent(xBin)-Hists[iSamp]->GetBinError(xBin));}
    }
  }
  
  for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {
    Hists[iSamp]->GetYaxis()->SetRangeUser(Min-0.1*(Max-Min),Max+0.1*(Max-Min));
  }

  for (size_t iSamp=0;iSamp<Hists.size();iSamp++) {

    if (iSamp==0) {
      Hists[iSamp]->Draw((DrawOptions_1D).c_str());
    } else {
      Hists[iSamp]->Draw((DrawOptions_1D+" SAME").c_str());
    }

    Legends[iSamp] = new TLegend(0.8,0.9-(1.0+static_cast<T>(iSamp))*LegendHeight,0.99,0.9-static_cast<T>(iSamp)*LegendHeight);
    Legends[iSamp]->SetTextSize(FontSize);
    Legends[iSamp]->AddEntry(Hists[iSamp],(Samples[iSamp]->GetName()).c_str(),"l");
    Legends[iSamp]->AddEntry((TObject*)0,Form("Entries: %4.5f",Hists[iSamp]->GetEntries()),"");
    Legends[iSamp]->AddEntry((TObject*)0,Form("Integral: %4.5f",Hists[iSamp]->Integral()),"");
    Legends[iSamp]->AddEntry((TObject*)0,Form("Mean: %4.5f",Hists[iSamp]->GetMean()),"");
    Legends[iSamp]->AddEntry((TObject*)0,Form("RMS: %4.5f",Hists[iSamp]->GetRMS()),"");
    Legends[iSamp]->Draw();
  }

  Canv->Print(OutputFileName_1D.c_str());
}

template<typename T>
T SampleManager<T>::CalculateCovariance(T XNominalBinContent, std::vector<T> XVariedBinContents, T YNominalBinContent, std::vector<T> YVariedBinContents) {
  T Covariance = 0.;
  
  size_t nVariedSamples = XVariedBinContents.size();
  if (nVariedSamples != YVariedBinContents.size()) {
    std::cerr << "Invalid shape of XVariedBinContents and YVariedBinContents" << std::endl;
    throw;
  }
  
  for (size_t iMeas=0;iMeas<nVariedSamples;iMeas++) {
    Covariance += (XVariedBinContents[iMeas]-XNominalBinContent)*(YVariedBinContents[iMeas]-YNominalBinContent);
  }
  Covariance /= nVariedSamples;

  return Covariance;
}

template<typename T>
void SampleManager<T>::PlotAnalysisBinning(YAML::Node Config) {
  OutputFileName_1D = Config["OutputName"].as<std::string>();
  DrawOptions_1D = Config["DrawOpts"].as<std::string>();
  std::string NominalSample = Config["NominalSample"].as<std::string>();
  SampleNameToRatioTo = NominalSample;

  LegendHeight = 0.2;
  if (Config["LegendHeight"]) {
    LegendHeight = Config["LegendHeight"].as<T>();
  }

  FontSize = _BAD_VALUE_;
  if (Config["LegendFontSize"]) {
    FontSize = Config["LegendFontSize"].as<T>();
  }

  IndexToRatioTo = -1;
  for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
    if (Samples[iSamp]->GetName() == SampleNameToRatioTo) {
      IndexToRatioTo = iSamp;
      break;
    }
  }
  if (IndexToRatioTo == -1) {
    std::cerr << "Did not find:" << SampleNameToRatioTo << std::endl;
    std::cerr << "Available:" << std::endl;
    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      std::cerr << "\t" << iSamp << " " << Samples[iSamp]->GetName() << std::endl;
    }
  }
  int NominalIndex = IndexToRatioTo;

  TCanvas* Canv = new TCanvas;
  Canv->SetRightMargin(0.2);
  Canv->SetLeftMargin(0.15);
  Canv->Print((OutputFileName_1D+"[").c_str());

  std::vector<TH1*> Hists(Samples.size());

  //===============================================================================
  //AnalysisBinning

  for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
    Hists[iSamp] = Samples[iSamp]->GetAnalysisBinningHistogram();
    Hists[iSamp]->SetStats(false);
  }

  Plot1DHists(Canv, Hists);
  Plot1DRatioHists(Canv, Hists);

  //===============================================================================
  //Calculate the covariance matrix
  
  int nBins = AnalysisBinning->GetNBins();
  TH2* CovarianceMatrix;
  TH2* CorrelationMatrix;
  TH1* CovarianceMatrixDiag;
  if (typeid(T) == typeid(float)) {
    CovarianceMatrix = new TH2F("AnalysisBinningCovMat","Covariance Matrix;Analysis Binning;Analysis Binning",nBins,0,nBins,nBins,0,nBins);
    CorrelationMatrix = new TH2F("AnalysisBinningCorrMat","Correlation Matrix;Analysis Binning;Analysis Binning",nBins,0,nBins,nBins,0,nBins);
    CovarianceMatrixDiag = new TH1F("AnalysisBinningCovMatDiag","Covariance Matrix;Analysis Binning;#sqrt{Covariance Diagonal}",nBins,0,nBins);
  } else {
    CovarianceMatrix = new TH2D("AnalysisBinningCovMat","Covariance Matrix;Analysis Binning;Analysis Binning",nBins,0,nBins,nBins,0,nBins);
    CorrelationMatrix = new TH2D("AnalysisBinningCorrMat","Correlation Matrix;Analysis Binning;Analysis Binning",nBins,0,nBins,nBins,0,nBins);
    CovarianceMatrixDiag = new TH1D("AnalysisBinningCovMatDiag","Covariance Matrix;Analysis Binning;#sqrt{Covariance Diagonal}",nBins,0,nBins);
  }

  for (int xBin=0;xBin<nBins;xBin++) {
    T XNominalBinContent = Samples[NominalIndex]->GetAnalysisBinningHistogram()->GetBinContent(xBin+1);

    std::vector<T> XVariedBinContents;
    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      if (static_cast<int>(iSamp) == NominalIndex) continue;
      XVariedBinContents.push_back(Samples[iSamp]->GetAnalysisBinningHistogram()->GetBinContent(xBin+1));
    }
    
    for (int yBin=0;yBin<nBins;yBin++) {
      T YNominalBinContent = Samples[NominalIndex]->GetAnalysisBinningHistogram()->GetBinContent(yBin+1);
      
      std::vector<T> YVariedBinContents;
      for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
	if (static_cast<int>(iSamp) == NominalIndex) continue;
	YVariedBinContents.push_back(Samples[iSamp]->GetAnalysisBinningHistogram()->GetBinContent(yBin+1));
      }
      
      T Covariance = CalculateCovariance(XNominalBinContent,XVariedBinContents,YNominalBinContent,YVariedBinContents);
      CovarianceMatrix->SetBinContent(xBin+1,yBin+1,Covariance);
    }
  }
  CovarianceMatrix->SetStats(false);
  CovarianceMatrix->Draw("COLZ");
  Canv->Print(OutputFileName_1D.c_str());

  CovarianceMatrixDiag->SetStats(false);
  for (int xBin=0;xBin<nBins;xBin++) {
    CovarianceMatrixDiag->SetBinContent(xBin+1,TMath::Sqrt(CovarianceMatrix->GetBinContent(xBin+1,xBin+1)));
  }
  CovarianceMatrixDiag->Draw();
  Canv->Print(OutputFileName_1D.c_str());

  CorrelationMatrix->SetStats(false);
  for (int xBin=0;xBin<nBins;xBin++) {
    for (int yBin=0;yBin<nBins;yBin++) {
      if (CovarianceMatrixDiag->GetBinContent(xBin+1) > 0 && CovarianceMatrixDiag->GetBinContent(yBin+1) > 0) {
	CorrelationMatrix->SetBinContent(xBin+1,yBin+1,CovarianceMatrix->GetBinContent(xBin+1,yBin+1)/(CovarianceMatrixDiag->GetBinContent(xBin+1)*CovarianceMatrixDiag->GetBinContent(yBin+1)));
      } else {
	CorrelationMatrix->SetBinContent(xBin+1,yBin+1,0.);
      }
    }
  }
  CorrelationMatrix->Draw("COLZ");
  Canv->Print(OutputFileName_1D.c_str());

  Canv->Print((OutputFileName_1D+"]").c_str());
}

template<typename T>
void SampleManager<T>::Plot1D(YAML::Node Config) {

  OutputFileName_1D = Config["OutputName"].as<std::string>();
  DrawOptions_1D = Config["DrawOpts"].as<std::string>();
  SampleNameToRatioTo = Config["RatioDenominatorSample"].as<std::string>();

  LegendHeight = 0.2;
  if (Config["LegendHeight"]) {
    LegendHeight = Config["LegendHeight"].as<T>();
  }

  FontSize = _BAD_VALUE_;
  if (Config["LegendFontSize"]) {
    FontSize = Config["LegendFontSize"].as<T>();
  }

  RatioYAxisMax = _BAD_VALUE_;
  if (Config["RatioMaximum"]) {
    RatioYAxisMax = Config["RatioMaximum"].as<T>();
  }
  RatioYAxisMin = _BAD_VALUE_;
  if (Config["RatioMinimum"]) {
    RatioYAxisMin = Config["RatioMinimum"].as<T>();
  }

  IndexToRatioTo = -1;
  for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
    if (Samples[iSamp]->GetName() == SampleNameToRatioTo) {
      IndexToRatioTo = iSamp;
      break;
    }
  }
  if (IndexToRatioTo == -1) {
    std::cerr << "Did not find:" << SampleNameToRatioTo << std::endl;
    std::cerr << "Available:" << std::endl;
    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      std::cerr << "\t" << iSamp << " " << Samples[iSamp]->GetName() << std::endl;
    }
  }

  TCanvas* Canv = new TCanvas;
  Canv->SetRightMargin(0.2);
  Canv->SetLeftMargin(0.15);
  Canv->Print((OutputFileName_1D+"[").c_str());

  std::vector<TH1*> Hists(Samples.size());

  //===============================================================================
  //Observations

  for (int iObs=0;iObs<ObsManager->GetNObservables();iObs++) {
    if (ObsManager->GetObservable(iObs)->GetNDimensions() != 1) continue;

    bool IsLog = ObsManager->GetObservable(iObs)->GetIsLog(0);
    Canv->SetLogx(IsLog);

    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      Hists[iSamp] = Samples[iSamp]->GetMeasurement(iObs);
      Hists[iSamp]->SetStats(false);
    }

    Plot1DHists(Canv, Hists);
    Plot1DRatioHists(Canv, Hists);
  }

  Canv->Print((OutputFileName_1D+"]").c_str());
}

template<typename T>
void SampleManager<T>::Plot2D(YAML::Node Config) {

  OutputFileName_2D = Config["OutputName"].as<std::string>();
  DrawOptions_2D = Config["DrawOpts"].as<std::string>();

  TCanvas* Canv = new TCanvas;
  Canv->SetRightMargin(0.2);
  Canv->SetLeftMargin(0.15);
  Canv->Print((OutputFileName_2D+"[").c_str());

  std::vector<TH1*> Hists(Samples.size());

  //===============================================================================
  //Observations

  for (int iObs=0;iObs<ObsManager->GetNObservables();iObs++) {
    if (ObsManager->GetObservable(iObs)->GetNDimensions() != 2) continue;

    bool IsLogX = ObsManager->GetObservable(iObs)->GetIsLog(0);
    Canv->SetLogx(IsLogX);

    bool IsLogY = ObsManager->GetObservable(iObs)->GetIsLog(1);
    Canv->SetLogy(IsLogY);

    for (size_t iSamp=0;iSamp<Samples.size();iSamp++) {
      Hists[iSamp] = Samples[iSamp]->GetMeasurement(iObs);
      Hists[iSamp]->SetTitle(Samples[iSamp]->GetName().c_str());
      Hists[iSamp]->SetStats(false);

      Hists[iSamp]->Draw(DrawOptions_2D.c_str());
      Canv->Print(OutputFileName_2D.c_str());
    }

  }

  Canv->Print((OutputFileName_2D+"]").c_str());
}
