#include <iostream>

#include "BartolFluxReader.h"
#include "HondaFluxReader.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TFile.h"

TH2D* InterpolateHistogram(TH2* ModelHistogram, TH2* TemplateHistogram,  std::string HistName="") {
  if (HistName == "") {
    HistName = ModelHistogram->GetName();
  }
  
  TH2D* InterpHistogram = (TH2D*)TemplateHistogram->Clone((HistName+"_Interp").c_str());
  for (int xBin=1;xBin<=InterpHistogram->GetNbinsX();xBin++) {
    for (int yBin=1;yBin<=InterpHistogram->GetNbinsY();yBin++) {
      double XVal = InterpHistogram->GetXaxis()->GetBinCenter(xBin);
      double YVal = InterpHistogram->GetYaxis()->GetBinCenter(yBin);

      InterpHistogram->SetBinContent(xBin,yBin,ModelHistogram->Interpolate(XVal,YVal));
    }
  }

  InterpHistogram->GetXaxis()->SetTitle(ModelHistogram->GetXaxis()->GetTitle());
  InterpHistogram->GetYaxis()->SetTitle(ModelHistogram->GetYaxis()->GetTitle());
  InterpHistogram->GetZaxis()->SetTitle(ModelHistogram->GetZaxis()->GetTitle());
  InterpHistogram->SetLineColor(ModelHistogram->GetLineColor());
  InterpHistogram->SetLineStyle(ModelHistogram->GetLineStyle());
  
  return InterpHistogram;
}

TGraph* TSpline3_to_TGraph(TSpline3* spline){
  // Convert TSpline3->TGraph
  int nEvals = 30*spline->GetNp();
  std::vector<double> x_points(nEvals+1);
  std::vector<double> y_points(nEvals+1);

  double SplineMin = 1e8;
  double SplineMax = -1e8;

  int nKnots = spline->GetNp();
  for (int iKnot=0;iKnot<nKnots;iKnot++) {
    double x = -999;
    double y = -999;
    spline->GetKnot(iKnot, x, y);

    if (x < SplineMin) {SplineMin = x;}
    if (x > SplineMax) {SplineMax = x;}
  }
  
  for(int iEval=0;iEval<=nEvals;iEval++){
    double x = SplineMin + ((double)iEval/(double)nEvals) * (SplineMax-SplineMin);
    double y = spline->Eval(x);

    x_points[iEval] = x;
    y_points[iEval] = y;
  }
  TGraph* return_graph = new TGraph(nEvals+1, x_points.data(), y_points.data());

  return return_graph;
}

int main(int argc, char const *argv[]) {
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetOptStat(false);
  
  //================================================================================================================================================================
  //Grab things from the config

  std::string TwoDRatioOutputName = "2DComparison.pdf";
  std::string RatioOutputName = "Comparison.pdf";
  std::string InterpOutputName = "InterpOutput.pdf";
  std::string FlavourRatioInterpOutputName = "FlavourRatioInterpOutput.pdf";
  std::string TwoDFlavourRatioComparisons = "2DFlavourRatioComparisons.pdf";
  std::string SplineOutputName = "Splines.pdf";
  std::string SplineFileOutputName = "AtmosphericFluxShapeSplines.root";
  
  if (argc != 2) {
    std::cout << "Usage: ./FluxComp Config.yaml" << std::endl;
    throw;
  }

  std::string ConfigName = argv[1];
  std::cout << "Config: " << ConfigName << std::endl;
  YAML::Node Config = YAML::LoadFile(ConfigName);

  if (!Config["FluxModels"]) {
    std::cerr << "No FluxModels to compare - Fix your YAML config" << std::endl;
    throw;
  }

  std::string NominalFluxModel = Config["General"]["Nominal"].as<std::string>();
  bool  FoundNominalModel = false;
  for (auto const &Model : Config["FluxModels"]) {
    std::string ModelName = Model["ModelName"].as<std::string>();
    std::cout << NominalFluxModel << " " << ModelName << std::endl;
    if (NominalFluxModel == ModelName) {
      FoundNominalModel = true;
      break;
    }
  }
  if (!FoundNominalModel) {
    std::cerr << "Did not find the nominal model in the defined FluxModels!: " << NominalFluxModel << std::endl;
    throw;
  }

  std::vector<double> ComparisonBinning_Energy = Config["General"]["ComparisonBinning_Energy"].as< std::vector<double> >();
  std::vector<double> ComparisonBinning_CosineZ = Config["General"]["ComparisonBinning_CosineZ"].as< std::vector<double> >();
  TH2D* InterpTemplateHistogram = new TH2D("HistogramTemplate_Interp","Template Histogram;Energy (GeV);Cosine Z",ComparisonBinning_Energy.size()-1,ComparisonBinning_Energy.data(),ComparisonBinning_CosineZ.size()-1,ComparisonBinning_CosineZ.data());
  
  TCanvas* Canv = new TCanvas("Canv","");
  Canv->SetLogx(true);
  Canv->SetRightMargin(0.2);

  //================================================================================================================================================================
  //Build the flux predictions
  
  std::vector<FluxReader*> Fluxes;
  std::vector< std::string > ModelsFound;
  for (auto const &Model : Config["FluxModels"]) {

    std::string ModelName = Model["ModelName"].as<std::string>();
    if (std::find(ModelsFound.begin(), ModelsFound.end(), ModelName) != ModelsFound.end()) {
      std::cerr << "Already found model: " << ModelName << " in Config!" << std::endl;
      throw;
    }
    
    if (ModelName == "BARTOLSolMin" || ModelName == "BARTOLSolMax") {
      FluxReader* BARTOL = new BartolFluxReader(Model);
      Fluxes.push_back(BARTOL);
    } else if (ModelName == "HONDASolMin" || ModelName == "BARTOL25SolMin" || ModelName == "HONDASolMax" || ModelName == "BARTOL25SolMax") {
      FluxReader* HONDA = new HondaFluxReader(Model);
      Fluxes.push_back(HONDA);
    } else {
      std::cerr << "Found ModelName: " << ModelName << " which has not been associated with a FluxReader" << std::endl;
      throw;
    }

    ModelsFound.push_back(ModelName);
  }
  for (auto Flux: Fluxes) {
    Flux->Plot2DFlux(Flux->GetModelName()+".pdf","COLZ");
  }

  int NominalFluxIndex = -1;
  for (size_t iFlux=0;iFlux<Fluxes.size();iFlux++) {
    if (Fluxes[iFlux]->GetModelName() == NominalFluxModel) {
      NominalFluxIndex = iFlux;
      break;
    }
  }
  if (NominalFluxIndex == -1) {
    std::cerr << "Did not find nominal model flux histograms!" << std::endl;
    throw;
  }

  //================================================================================================================================================================
  //Compare the single flavour predictions

  //================================================================================
  //Compare the flux predictions
  
  std::vector< std::vector<TH2*> > SingleFlavourFluxHists(nFlavs);

  Canv->Print((InterpOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    SingleFlavourFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      SingleFlavourFluxHists[iFlav][iModel] = InterpolateHistogram((TH2*)(Fluxes[iModel]->ReturnEnergyCosineZHists())[iFlav], InterpTemplateHistogram);
      SingleFlavourFluxHists[iFlav][iModel]->SetTitle((std::string(SingleFlavourFluxHists[iFlav][iModel]->GetTitle())+" Interpolated "+Fluxes[iModel]->GetModelName()).c_str());
      
      SingleFlavourFluxHists[iFlav][iModel]->Draw("COLZ");
      Canv->Print(InterpOutputName.c_str());
    }
  }
  Canv->Print((InterpOutputName+"]").c_str());

  //================================================================================
  //Make Ratio Plots

  std::vector< std::vector<TH2*> > RatioFluxHists(nFlavs);
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    RatioFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFluxHists[iFlav][iModel] = (TH2*)SingleFlavourFluxHists[iFlav][iModel]->Clone();
      RatioFluxHists[iFlav][iModel]->Divide(SingleFlavourFluxHists[iFlav][NominalFluxIndex]);
    }
  }

  Canv->Print((TwoDRatioOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFluxHists[iFlav][iModel]->SetTitle(("Ratio of Interpolated "+Fluxes[iModel]->GetModelName()+" to Interpolated "+NominalFluxModel+" : "+FlavourNames.at(iFlav)).c_str());
      RatioFluxHists[iFlav][iModel]->GetZaxis()->SetTitle("Ratio");
      
      RatioFluxHists[iFlav][iModel]->Draw("COLZ");
      Canv->Print((TwoDRatioOutputName).c_str());
    }
  }
  Canv->Print((TwoDRatioOutputName+"]").c_str());

  //================================================================================================================================================================
  //Compare the flavour ratio predictions

  //================================================================================
  //Compare the flux predictions
  
  std::vector< std::vector<TH2*> > FlavourRatioFluxHists(nFlavRatios);

  Canv->Print((FlavourRatioInterpOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    FlavourRatioFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      FlavourRatioFluxHists[iFlav][iModel] = InterpolateHistogram((TH2*)(Fluxes[iModel]->ReturnFlavourRatioHists())[iFlav], InterpTemplateHistogram);
      FlavourRatioFluxHists[iFlav][iModel]->Smooth();
      FlavourRatioFluxHists[iFlav][iModel]->SetTitle((std::string(FlavourRatioFluxHists[iFlav][iModel]->GetTitle())+" Interpolated "+Fluxes[iModel]->GetModelName()).c_str());

      FlavourRatioFluxHists[iFlav][iModel]->Draw("COLZ");
      Canv->Print(FlavourRatioInterpOutputName.c_str());
    }
  }
  Canv->Print((FlavourRatioInterpOutputName+"]").c_str());
  //================================================================================
  //Make Ratio Plots

  std::vector< std::vector<TH2*> > RatioFlavourRatioFluxHists(nFlavRatios);
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    RatioFlavourRatioFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFlavourRatioFluxHists[iFlav][iModel] = (TH2*)FlavourRatioFluxHists[iFlav][iModel]->Clone();
      RatioFlavourRatioFluxHists[iFlav][iModel]->Divide(FlavourRatioFluxHists[iFlav][NominalFluxIndex]);
      RatioFlavourRatioFluxHists[iFlav][iModel]->Smooth();
    }
  }

  Canv->Print((TwoDFlavourRatioComparisons+"[").c_str());
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFlavourRatioFluxHists[iFlav][iModel]->SetTitle(("Ratio of Interpolated "+Fluxes[iModel]->GetModelName()+" to Interpolated "+NominalFluxModel+" : "+RatioFlavourNames.at(iFlav)).c_str());
      RatioFlavourRatioFluxHists[iFlav][iModel]->GetZaxis()->SetTitle("Ratio");
      
      RatioFlavourRatioFluxHists[iFlav][iModel]->Draw("COLZ TEXT0");
      Canv->Print((TwoDFlavourRatioComparisons).c_str());
    }
  }
  Canv->Print((TwoDFlavourRatioComparisons+"]").c_str());
  
  //========================================================================================================================================================================================
  //Build the Splines

  std::vector<FLOAT_T> ModelSplineIndex(Fluxes.size());
  std::vector<FLOAT_T> SplineXVals(Fluxes.size());
  int Counter = 0;

  //First set SplineXVal==0 as NominalModel
  ModelSplineIndex[Counter] = NominalFluxIndex;
  SplineXVals[Counter] = Counter;
  Counter += 1;
  //Now count the others
  for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
    if ((int)iModel == NominalFluxIndex) continue;
    ModelSplineIndex[Counter] = iModel;
    SplineXVals[Counter] = Counter;
    Counter += 1;
  }

  std::cout << "Spline XValues:" << std::endl;
  for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
    std::cout << "\t" << iModel << " = " << Fluxes[ModelSplineIndex[iModel]]->GetModelName() << std::endl;
  }

  TFile* File = TFile::Open(SplineFileOutputName.c_str(),"RECREATE");
  File->cd();
  InterpTemplateHistogram->Write();

  TLegend* Leg = new TLegend(0.8,0.9-0.05*Fluxes.size(),0.99,0.9);
  Leg->SetTextSize(0.018);
  for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
    Leg->AddEntry((TObject*)0,(Form("%4.1f = ",SplineXVals[iModel])+Fluxes[ModelSplineIndex[iModel]]->GetModelName()).c_str(),"");
  }

  Canv->SetLogx(false);
  Canv->Print((SplineOutputName+"[").c_str());  
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    for (int xBin=1;xBin<=InterpTemplateHistogram->GetNbinsX();xBin++) {
      for (int yBin=1;yBin<=InterpTemplateHistogram->GetNbinsY();yBin++) {


	std::vector<FLOAT_T> SplineYVals(Fluxes.size());
	for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
	  SplineYVals[iModel] = RatioFlavourRatioFluxHists[iFlav][ModelSplineIndex[iModel]]->GetBinContent(xBin,yBin);
	}

	TSpline3* Spline = new TSpline3(Form("Syst_%i_XBin_%i_YBin_%i",iFlav,xBin,yBin),SplineXVals.data(),SplineYVals.data(),Fluxes.size());
	TGraph* Graph = TSpline3_to_TGraph(Spline);
	Graph->SetTitle((std::string(Spline->GetTitle())+";Dial Value;Weight").c_str());

	Canv->Clear();
	Graph->Draw();
	Leg->Draw("SAME");
	Canv->Print((SplineOutputName).c_str());

	Spline->Write(Form("Syst_%s_XBin_%i_YBin_%i",RatioFlavourShortNames[iFlav].c_str(),xBin,yBin));
      }
    }
    
  }
  Canv->Print((SplineOutputName+"]").c_str());
  File->Write();

  return 0;
}
