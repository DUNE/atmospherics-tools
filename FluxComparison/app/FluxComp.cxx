#include <iostream>

#include "Utils.h"
#include "Constants.h"
#include "BartolFluxReader.h"
#include "HondaFluxReader.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TFile.h"

int main(int argc, char const *argv[]) {
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetOptStat(false);

  TCanvas* Canv = new TCanvas("Canv","");
  Canv->SetLogx(true);
  Canv->SetRightMargin(0.2);
  
  //================================================================================================================================================================
  //Grab things from the config

  if (argc != 2) {
    std::cout << "Usage: ./FluxComp Config.yaml" << std::endl;
    throw;
  }

  std::string ConfigName = argv[1];
  std::cout << "Config: " << ConfigName << std::endl;
  YAML::Node Config = YAML::LoadFile(ConfigName);

  //===============
  //File Names
  
  std::string SingleModelFlavourPlotPrefix = Config["General"]["SingleModelFlavourPlotPrefix"].as<std::string>();
  std::string SingleModelFlavourPlotSuffix = Config["General"]["SingleModelFlavourPlotSuffix"].as<std::string>();
  std::string SingleModelFlavourPlotDrawOptions = Config["General"]["SingleModelFlavourPlotDrawOptions"].as<std::string>();

  std::string InterpSingleModelFlavourPlotOutputName = Config["General"]["InterpSingleModelFlavourPlotOutputName"].as<std::string>();
  std::string InterpSingleModelFlavourPlotDrawOptions = Config["General"]["InterpSingleModelFlavourPlotDrawOptions"].as<std::string>();

  std::string RatioInterpModelFlavourPlotOutputName = Config["General"]["RatioInterpModelFlavourPlotOutputName"].as<std::string>();
  std::string RatioInterpModelFlavourPlotDrawOptions = Config["General"]["RatioInterpModelFlavourPlotDrawOptions"].as<std::string>();

  std::string InterpSingleModelFlavourCombPlotOutputName = Config["General"]["InterpSingleModelFlavourCombPlotOutputName"].as<std::string>();
  std::string InterpSingleModelFlavourCombPlotDrawOptions = Config["General"]["InterpSingleModelFlavourCombPlotDrawOptions"].as<std::string>();

  std::string RatioInterpModelFlavourCombPlotOutputName = Config["General"]["RatioInterpModelFlavourCombPlotOutputName"].as<std::string>();
  std::string RatioInterpModelFlavourCombPlotDrawOptions = Config["General"]["RatioInterpModelFlavourCombPlotDrawOptions"].as<std::string>();

  std::string SplinePlotOutputName = Config["General"]["SplinePlotOutputName"].as<std::string>();
  std::string SplineFileOutputName = Config["General"]["SplineFileOutputName"].as<std::string>();

  //===============
  //Flux Models
  
  if (!Config["FluxModels"]) {
    std::cerr << "No FluxModels to compare - Fix your YAML config" << std::endl;
    throw;
  }

  std::string NominalFluxModel = Config["General"]["NominalModel"].as<std::string>();
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

  //===============
  //Comparison Binning

  std::vector<double> ComparisonBinning_Energy = Config["General"]["ComparisonBinning_Energy"].as< std::vector<double> >();
  std::vector<double> ComparisonBinning_CosineZ = Config["General"]["ComparisonBinning_CosineZ"].as< std::vector<double> >();
  TH2D* InterpTemplateHistogram = new TH2D("HistogramTemplate_Interp","Template Histogram;Energy (GeV);Cosine Z",ComparisonBinning_Energy.size()-1,ComparisonBinning_Energy.data(),ComparisonBinning_CosineZ.size()-1,ComparisonBinning_CosineZ.data());
  
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
    Flux->Plot2DFlux((SingleModelFlavourPlotPrefix+Flux->GetModelName()+SingleModelFlavourPlotSuffix).c_str(),SingleModelFlavourPlotDrawOptions);
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

  Canv->Print((InterpSingleModelFlavourPlotOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    SingleFlavourFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      SingleFlavourFluxHists[iFlav][iModel] = InterpolateHistogram((TH2*)(Fluxes[iModel]->ReturnEnergyCosineZHists())[iFlav], InterpTemplateHistogram);
      SingleFlavourFluxHists[iFlav][iModel]->SetTitle((std::string("Interpolated ")+Fluxes[iModel]->GetModelName()+" : "+FlavourNames.at(iFlav)).c_str());
      
      SingleFlavourFluxHists[iFlav][iModel]->Draw(InterpSingleModelFlavourPlotDrawOptions.c_str());
      Canv->Print(InterpSingleModelFlavourPlotOutputName.c_str());
    }
  }
  Canv->Print((InterpSingleModelFlavourPlotOutputName+"]").c_str());

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

  Canv->Print((RatioInterpModelFlavourPlotOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFluxHists[iFlav][iModel]->SetTitle(("Ratio of Interpolated "+Fluxes[iModel]->GetModelName()+" to Interpolated "+NominalFluxModel+" : "+FlavourNames.at(iFlav)).c_str());
      RatioFluxHists[iFlav][iModel]->GetZaxis()->SetTitle("Ratio");
      
      RatioFluxHists[iFlav][iModel]->Draw(RatioInterpModelFlavourPlotDrawOptions.c_str());
      Canv->Print((RatioInterpModelFlavourPlotOutputName).c_str());
    }
  }
  Canv->Print((RatioInterpModelFlavourPlotOutputName+"]").c_str());

  //================================================================================================================================================================
  //Compare the flavour ratio predictions

  //================================================================================
  //Compare the flux predictions
  
  std::vector< std::vector<TH2*> > FlavourRatioFluxHists(nFlavRatios);

  Canv->Print((InterpSingleModelFlavourCombPlotOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    FlavourRatioFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      FlavourRatioFluxHists[iFlav][iModel] = InterpolateHistogram((TH2*)(Fluxes[iModel]->ReturnFlavourRatioHists())[iFlav], InterpTemplateHistogram);
      FlavourRatioFluxHists[iFlav][iModel]->Smooth();
      FlavourRatioFluxHists[iFlav][iModel]->SetTitle((std::string("Interpolated ")+Fluxes[iModel]->GetModelName()+" : "+RatioFlavourNames.at(iFlav)).c_str());

      FlavourRatioFluxHists[iFlav][iModel]->Draw(InterpSingleModelFlavourCombPlotDrawOptions.c_str());
      Canv->Print(InterpSingleModelFlavourCombPlotOutputName.c_str());
    }
  }
  Canv->Print((InterpSingleModelFlavourCombPlotOutputName+"]").c_str());
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

  Canv->Print((RatioInterpModelFlavourCombPlotOutputName+"[").c_str());
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      RatioFlavourRatioFluxHists[iFlav][iModel]->SetTitle(("Ratio of Interpolated "+Fluxes[iModel]->GetModelName()+" to Interpolated "+NominalFluxModel+" : "+RatioFlavourNames.at(iFlav)).c_str());
      RatioFlavourRatioFluxHists[iFlav][iModel]->GetZaxis()->SetTitle("Ratio");
      
      RatioFlavourRatioFluxHists[iFlav][iModel]->Draw(RatioInterpModelFlavourCombPlotDrawOptions.c_str());
      Canv->Print((RatioInterpModelFlavourCombPlotOutputName).c_str());
    }
  }
  Canv->Print((RatioInterpModelFlavourCombPlotOutputName+"]").c_str());
  
  //========================================================================================================================================================================================
  //Build the Splines

  size_t nSplineKnots = 3*Fluxes.size();
  
  std::vector<FLOAT_T> ModelSplineIndex(nSplineKnots);
  std::vector<FLOAT_T> SplineXVals(nSplineKnots);

  for (size_t iModel=0;iModel<nSplineKnots;iModel++) {
    ModelSplineIndex[iModel] = iModel % Fluxes.size();
    SplineXVals[iModel] = iModel;
  }

  std::cout << "Spline XValues:" << std::endl;
  for (size_t iModel=0;iModel<nSplineKnots;iModel++) {
    std::cout << "\t" << iModel << " = " << Fluxes[ModelSplineIndex[iModel]]->GetModelName() << std::endl;
  }

  TFile* File = TFile::Open(SplineFileOutputName.c_str(),"RECREATE");
  File->cd();
  InterpTemplateHistogram->Write();

  TLegend* Leg = new TLegend(0.8,0.9-0.05*(Fluxes.size()+1),0.99,0.9);
  Leg->SetTextSize(0.018);
  for (size_t iModel=0;iModel<nSplineKnots;iModel++) {
    if ((iModel < Fluxes.size()) || (iModel > 2*Fluxes.size())) continue;
    Leg->AddEntry((TObject*)0,(Form("%4.1f = ",SplineXVals[iModel])+Fluxes[ModelSplineIndex[iModel]]->GetModelName()).c_str(),"");
  }

  Canv->SetLogx(false);
  Canv->Print((SplinePlotOutputName+"[").c_str());  
  for (int iFlav=0;iFlav<nFlavRatios;iFlav++) {
    for (int xBin=1;xBin<=InterpTemplateHistogram->GetNbinsX();xBin++) {
      for (int yBin=1;yBin<=InterpTemplateHistogram->GetNbinsY();yBin++) {


	std::vector<FLOAT_T> SplineYVals(nSplineKnots);
	for (size_t iModel=0;iModel<nSplineKnots;iModel++) {
	  SplineYVals[iModel] = RatioFlavourRatioFluxHists[iFlav][ModelSplineIndex[iModel]]->GetBinContent(xBin,yBin);
	}

	TSpline3* Spline = new TSpline3(Form("Syst_%i_XBin_%i_YBin_%i",iFlav,xBin,yBin),SplineXVals.data(),SplineYVals.data(),nSplineKnots);
	TGraph* Graph = TSpline3_to_TGraph(Spline);
	Graph->SetTitle((std::string(Spline->GetTitle())+";Dial Value;Weight").c_str());
	Graph->GetXaxis()->SetLimits(Fluxes.size(),2*Fluxes.size());
	
	Canv->Clear();
	Graph->Draw();
	Leg->Draw("SAME");
	Canv->Print((SplinePlotOutputName).c_str());

	Spline->Write(Form("Syst_%s_XBin_%i_YBin_%i",RatioFlavourShortNames[iFlav].c_str(),xBin,yBin));
      }
    }
    
  }
  Canv->Print((SplinePlotOutputName+"]").c_str());
  File->Write();

  return 0;
}
