#include <iostream>

#include "BartolFluxReader.h"
#include "HondaFluxReader.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"

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

// fit polynomial to absolute envelope
// This function fits a polynomial to the absolute envelope of the ratios of alternative flux models to the nominal model.
// It takes the absolute deviation histogram, fits a polynomial to it, and saves the fit

//it is called inside the function below PlotEnvelope
void FitAbsoluteEnvelope(TH1* hAbsDev, const std::string& outName, const std::string& xLabel)
{
  TF1* fitFunc = new TF1("fitFunc", "pol2", hAbsDev->GetXaxis()->GetXmin(), hAbsDev->GetXaxis()->GetXmax());
  hAbsDev->Fit(fitFunc, "R");

  TCanvas* c = new TCanvas("c", "Fit", 800, 600);
  hAbsDev->SetLineColor(kMagenta + 2);
  hAbsDev->SetTitle(("Fit to Absolute Envelope vs " + xLabel).c_str());
  hAbsDev->GetXaxis()->SetTitle(xLabel.c_str());
  hAbsDev->GetYaxis()->SetTitle("Max |Model/Nominal - 1|");
  hAbsDev->Draw("hist");
  fitFunc->SetLineColor(kGreen + 2);
  fitFunc->Draw("same");

  c->SaveAs((outName + "_absdev_fit.pdf").c_str());

  // Optional: write to file
  TFile* fout = new TFile((outName + "_fit.root").c_str(), "RECREATE");
  fitFunc->Write("AbsDevPolyFit");
  fout->Close();

  delete c;
}

//This function plots the envelope of the ratios of alternative flux models to the nominal model for a given flavor
//It takes the nominal histogram and a vector of alternative histograms, and plots the maximum and minimum ratios
//as a function of energy and cosine theta. The envelope is plotted as two histograms, one for the maximum ratio and one for the minimum ratio.
//The function is called for each flavor and dimension (energy and cosine theta).
void PlotEnvelope(int iFlav,
                  const std::string& flavLabel,
                  const TH2* hNominal,
                  const std::vector<TH2*>& altHists,
                  TCanvas* canvas,
                  const std::string& outPrefix)
{
  for (int iDim = 0; iDim < 2; ++iDim) {
    // Nominal projection
    TH1* hNomProj = (iDim == 0) ? hNominal->ProjectionX() : hNominal->ProjectionY();
    hNomProj->SetDirectory(nullptr);

    const int nBins = hNomProj->GetNbinsX();

    // Prepare ratio envelopes
    TH1* hMaxRatio = (TH1*)hNomProj->Clone("hMaxRatio");
    TH1* hMinRatio = (TH1*)hNomProj->Clone("hMinRatio");
    TH1* hAbsDev   = (TH1*)hNomProj->Clone("hAbsDev");

    hMaxRatio->Reset();
    hMinRatio->Reset();
    hAbsDev->Reset();

    for (int ibin = 1; ibin <= nBins; ++ibin) {
      const double nomVal = hNomProj->GetBinContent(ibin);
      if (nomVal <= 0) {
        hMaxRatio->SetBinContent(ibin, 1.0);
        hMinRatio->SetBinContent(ibin, 1.0);
        hAbsDev->SetBinContent(ibin, 0.0);
        continue;
      }

      double maxRatio = -1e9;
      double minRatio = +1e9;
      double maxAbsDev = 0.0;

      for (const auto& hAlt2D : altHists) {
        TH1* hAltProj = (iDim == 0) ? hAlt2D->ProjectionX() : hAlt2D->ProjectionY();
        hAltProj->SetDirectory(nullptr);

        const double altVal = hAltProj->GetBinContent(ibin);
        const double ratio = altVal / nomVal;
        const double absDev = std::abs(ratio - 1.0);

        if (ratio > maxRatio) maxRatio = ratio;
        if (ratio < minRatio) minRatio = ratio;
        if (absDev > maxAbsDev) maxAbsDev = absDev;

        delete hAltProj;
      }

      hMaxRatio->SetBinContent(ibin, maxRatio);
      hMinRatio->SetBinContent(ibin, minRatio);
      hAbsDev->SetBinContent(ibin, 1.0 + maxAbsDev); // draw on same scale
    }

    // Draw envelope
    canvas->Clear();
    hMaxRatio->SetLineColor(kRed);
    hMinRatio->SetLineColor(kBlue);
    hAbsDev->SetLineColor(kMagenta);
    hMaxRatio->SetLineWidth(2);
    hMinRatio->SetLineWidth(2);
    hAbsDev->SetLineWidth(2);

    hMaxRatio->SetTitle(Form("Envelope for %s vs %s", flavLabel.c_str(), iDim == 0 ? "E" : "cos(#theta)"));
    hMaxRatio->GetYaxis()->SetTitle("Model / Nominal");
    hMaxRatio->GetXaxis()->SetTitle(iDim == 0 ? "Energy (GeV)" : "cos(#theta)");
    hMaxRatio->GetYaxis()->SetRangeUser(0.5, 1.5);

    hMaxRatio->Draw("hist");
    hMinRatio->Draw("hist same");
    hAbsDev->Draw("hist same");

    TLegend* leg = new TLegend(0.15, 0.75, 0.5, 0.9);
    leg->AddEntry(hMaxRatio, "Max(Model/Nominal)", "l");
    leg->AddEntry(hMinRatio, "Min(Model/Nominal)", "l");
    leg->AddEntry(hAbsDev,   "Max|Model/Nominal - 1|", "l");
    leg->Draw();

    canvas->SaveAs((outPrefix + Form("_%s_env_%s.pdf", flavLabel.c_str(), iDim == 0 ? "E" : "cos")).c_str());


    // Fit polynomial to absolute envelope, calling the function FitAbsoluteEnvelope
    // This function fits a polynomial to the absolute envelope of the ratios of alternative flux models to the nominal model.
    // It takes the absolute deviation histogram, fits a polynomial to it, and saves the fit
    std::string xLabel = (iDim == 0) ? "Energy (GeV)" : "cos(#theta)";
    std::string outName = outPrefix + Form("_%s_env_%s", flavLabel.c_str(), iDim == 0 ? "E" : "cos");
    FitAbsoluteEnvelope(hAbsDev, outName, xLabel);



    // Clean up
    delete hMaxRatio;
    delete hMinRatio;
    delete hAbsDev;
    delete hNomProj;
  }
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

  // Reads which model to use as the reference/nominal
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
  
  // Reads the common interpolation binning for comparison
  std::vector<double> ComparisonBinning_Energy = Config["General"]["ComparisonBinning_Energy"].as< std::vector<double> >();
  std::vector<double> ComparisonBinning_CosineZ = Config["General"]["ComparisonBinning_CosineZ"].as< std::vector<double> >();

  //================================================================================
  //Build the flux predictions
  
  std::vector<FluxReader*> Fluxes; //iterate through flux models in FluxModels
  for (auto const &Model : Config["FluxModels"]) {

    std::string ModelName = Model["ModelName"].as<std::string>();
    if (std::find(ModelsFound.begin(), ModelsFound.end(), ModelName) != ModelsFound.end()) {
      std::cerr << "Already found model: " << ModelName << " in Config!" << std::endl;
      throw;
    }
    
    //Bartol flux reader and honda flux reader are different. Honda reads single table file while bartol reads separate file per flavor
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

  // loop over all flavors and flux models; interpolate 2D flux hist to common binning
  // store interpolated hist in nested vector AlternativeFluxHists[iFlav][iModel]
  // draw each interpolated hist and save to file
  // then draw ratio plots for each flavor and each dimension (x and y) of the flux histograms
  // ratio plots are drawn for each model compared to the nominal model
  // save ratio plots to file
  // finally, print the canvas with all the ratio plots
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
	Canv->SetLogx(true); // project along energy axis -> flux vs energy
      } else {
	Canv->SetLogx(false); // project along cosinetheta axis -> flux vs cosinetheta
      }

      for (size_t iModel=0;iModel<Fluxes.size();iModel++) {
	TH2* AlternativeModelFluxHist = AlternativeFluxHists[iFlav][iModel];
	TH2* NominalModelFluxHist = NominalModelFluxHists[iFlav];
	
	TH1* NominalModelFluxHist_Proj;
	TH1* AlternativeModelFluxHist_Proj;

  // previously initializing the histograms
  // but now we will project them directly from the 2D histograms
	if (iDim==0) { // project along energy axis -> flux vs energy
	  NominalModelFluxHist_Proj = (TH1*)NominalModelFluxHist->ProjectionX();
	  AlternativeModelFluxHist_Proj = (TH1*)AlternativeModelFluxHist->ProjectionX();

	  AlternativeModelFluxHist_Proj->GetXaxis()->SetTitle(AlternativeModelFluxHist->GetXaxis()->GetTitle());
	} else { // project along cosinetheta axis -> flux vs cosinetheta
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

// Plot the envelope of the ratios for each flavor
for (int iFlav = 0; iFlav < nFlavs; ++iFlav) {
PlotEnvelope(iFlav,
              Fluxes[0]->GetFlavourName(iFlav),
              NominalModelFluxHists[iFlav],
              AlternativeFluxHists[iFlav],
              Canv,
              "EnvelopePlot");
}
  
  return 0;
}
