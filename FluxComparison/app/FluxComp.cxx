#include <iostream>

#include "BartolFluxReader.h"
#include "HondaFluxReader.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

TH2D* InterpolateHistogram(TH2* ModelHistogram, std::vector<double> EnergyBinning, std::vector<double> CosineZBinning, std::string HistName="") {
  if (HistName == "") {
    HistName = ModelHistogram->GetName();
  }
  
  TH2D* InterpHistogram = new TH2D((HistName+"_Interp").c_str(),ModelHistogram->GetTitle(),EnergyBinning.size()-1,EnergyBinning.data(),CosineZBinning.size()-1,CosineZBinning.data());
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

int main(int argc, char const *argv[]) {
  gStyle->SetOptStat(false);
  
  //================================================================================
  //Grab things from the config

  std::string RatioOutputName = "Comparison.pdf";
  std::string InterpOutputName = "InterpOutput.pdf";
  
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
  std::vector< std::string > ModelsFound;

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

  //================================================================================
  //Build the flux predictions
  
  std::vector<FluxReader*> Fluxes;
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

  //================================================================================
  //Compare the flux predictions
  
  std::vector<TH2*> NominalModelFluxHists(nFlavs);
  for (size_t iFlux=0;iFlux<Fluxes.size();iFlux++) {
    if (Fluxes[iFlux]->GetModelName() == NominalFluxModel) {
      for (int iFlav=0;iFlav<nFlavs;iFlav++) {
	TH1* Hist = (Fluxes[iFlux]->ReturnEnergyCosineZHists())[iFlav];
	NominalModelFluxHists[iFlav] = InterpolateHistogram((TH2*)Hist, ComparisonBinning_Energy, ComparisonBinning_CosineZ, std::string(Hist->GetName())+"_Nom");
      }
      break;
    }
  }
  if (NominalModelFluxHists.size() == 0) {
    std::cerr << "Did not find nominal model flux histograms!" << std::endl;
    throw;
  }

  TCanvas* Canv = new TCanvas;
  Canv->SetLogx(true);
  Canv->SetRightMargin(0.2);
  
  Canv->Print((InterpOutputName+"[").c_str());
  std::vector< std::vector<TH2*> > AlternativeFluxHists(nFlavs);
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    AlternativeFluxHists[iFlav].resize(Fluxes.size());
    for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
      AlternativeFluxHists[iFlav][iModel] = InterpolateHistogram((TH2*)(Fluxes[iModel]->ReturnEnergyCosineZHists())[iFlav], ComparisonBinning_Energy, ComparisonBinning_CosineZ);
      AlternativeFluxHists[iFlav][iModel]->SetTitle((std::string(AlternativeFluxHists[iFlav][iModel]->GetTitle())+" Interpolated "+Fluxes[iModel]->GetModelName()).c_str());
      
      AlternativeFluxHists[iFlav][iModel]->Draw("COLZ");
      Canv->Print(InterpOutputName.c_str());
    }
  }
  Canv->Print((InterpOutputName+"]").c_str());

  Canv->SetRightMargin(0.1);
  Canv->Print((RatioOutputName+"[").c_str());
  TLegend* Leg = new TLegend(0.1,0.9,0.9,1.0);
  Leg->SetNColumns(2);
  for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
    Leg->AddEntry(AlternativeFluxHists[0][iModel],Fluxes[iModel]->GetModelName().c_str(),"l");
  }
  
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    for (size_t iDim=0;iDim<2;iDim++) {
      if (iDim==0) {
	Canv->SetLogx(true);
      } else {
	Canv->SetLogx(false);
      }

      for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
	TH2* AlternativeModelFluxHist = AlternativeFluxHists[iFlav][iModel];
	TH2* NominalModelFluxHist = NominalModelFluxHists[iFlav];
	
	TH1* NominalModelFluxHist_Proj;
	TH1* AlternativeModelFluxHist_Proj;

	if (iDim==0) {
	  NominalModelFluxHist_Proj = (TH1*)NominalModelFluxHist->ProjectionX();
	  AlternativeModelFluxHist_Proj = (TH1*)AlternativeModelFluxHist->ProjectionX();

	  AlternativeModelFluxHist_Proj->GetXaxis()->SetTitle(AlternativeModelFluxHist->GetXaxis()->GetTitle());
	} else {
	  NominalModelFluxHist_Proj = (TH1*)NominalModelFluxHist->ProjectionY();
          AlternativeModelFluxHist_Proj = (TH1*)AlternativeModelFluxHist->ProjectionY();

	  AlternativeModelFluxHist_Proj->GetXaxis()->SetTitle(AlternativeModelFluxHist->GetYaxis()->GetTitle());
	}
	AlternativeModelFluxHist_Proj->SetTitle("");
	AlternativeModelFluxHist_Proj->GetYaxis()->SetTitle((std::string("Ratio of ")+AlternativeModelFluxHist->GetZaxis()->GetTitle()).c_str());

	AlternativeModelFluxHist_Proj->Divide(NominalModelFluxHist_Proj);
	AlternativeModelFluxHist_Proj->GetYaxis()->SetRangeUser(0.5,1.5);
	AlternativeModelFluxHist_Proj->SetLineStyle(AlternativeFluxHists[0][iModel]->GetLineStyle());

	if (iModel==0) {
	  AlternativeModelFluxHist_Proj->Draw();
	} else {
	  AlternativeModelFluxHist_Proj->Draw("SAME");
	}
      }

      Leg->Draw("SAME");
      
      TLatex FlavTitle(.8,.8,(Fluxes[0]->GetFlavourName(iFlav)).c_str());  
      FlavTitle.SetNDC(kTRUE);
      FlavTitle.Draw("SAME");
      
      Canv->Print(RatioOutputName.c_str());
    }    
  }

  Canv->Print((RatioOutputName+"]").c_str());
  
  return 0;
}
