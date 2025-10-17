#include "FluxReader.h"

#include <iostream>

#include "TCanvas.h"

FluxReader::FluxReader(YAML::Node Config_) {
  std::cout << std::endl;
  Config = Config_;
  ModelName = Config["ModelName"].as<std::string>();

  std::cout << "\nInitialising Flux:" << ModelName << std::endl;
  
  EnergyBinEdges = Config["EnergyBinEdges"].as< std::vector<FLOAT_T> >();
  std::cout << "Energy Bin Edges:" << std::endl;
  std::cout << "[";
  for (size_t iEn=0;iEn<EnergyBinEdges.size();iEn++) {
    std::cout << EnergyBinEdges[iEn] << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << std::endl;

  CosineZBinEdges = Config["CosineZBinEdges"].as< std::vector<FLOAT_T> >();
  std::cout << "CosineZ Bin Edges:" << std::endl;
  std::cout << "[";
  for (size_t iCZ=0;iCZ<CosineZBinEdges.size();iCZ++) {
    std::cout << CosineZBinEdges[iCZ] << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << std::endl;

  if (Config["PhiBinEdges"]) {
    MeasDimension = 3;

    PhiBinEdges = Config["PhiBinEdges"].as< std::vector<FLOAT_T> >();
    std::cout << "Phi Bin Edges:" << std::endl;
    std::cout << "[";
    for (size_t iPhi=0;iPhi<PhiBinEdges.size();iPhi++) {
      std::cout << PhiBinEdges[iPhi] << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;

  } else {
    MeasDimension = 2;
  }

  if (Config["Smooth"]) {
    Smooth = Config["Smooth"].as<bool>();
  }
  if (Smooth) {std::cout << "Warning!: Smoothing 2D distribution\n" << std::endl;}

  EnergyAxisMin = Config["EnergyAxisMin"].as<FLOAT_T>();
  EnergyAxisMax = Config["EnergyAxisMax"].as<FLOAT_T>();
  LineColor = Config["LineColor"].as<int>();
  LineStyle = Config["LineStyle"].as<int>();
  
  FluxCaption = Config["FluxCaption"].as<std::string>();
  std::string HistAxisCaptions = Config["HistAxisCaptions"].as<std::string>();
  
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    std::string HistName = ModelName+"_"+FlavourNames[iFlav];    
    std::string HistTitle = FlavourNames[iFlav]+";"+HistAxisCaptions;

    if (MeasDimension == 2) {
      FluxHists[iFlav] = new TH2D(HistName.c_str(),HistTitle.c_str(),EnergyBinEdges.size()-1,EnergyBinEdges.data(),CosineZBinEdges.size()-1,CosineZBinEdges.data());
      FluxHists[iFlav]->GetZaxis()->SetTitleOffset(1.3);
    } else if (MeasDimension == 3) {
      FluxHists[iFlav] = new TH3D(HistName.c_str(),HistTitle.c_str(),EnergyBinEdges.size()-1,EnergyBinEdges.data(),CosineZBinEdges.size()-1,CosineZBinEdges.data(),PhiBinEdges.size()-1,PhiBinEdges.data());
    } else {
      std::cerr << "Unsupported number of dimensions!" << std::endl;
    }
    FluxHists[iFlav]->SetStats(false);
    
  }

}

void FluxReader::InitialiseFlux() {
  std::cout << "Reading flux from tables for model:" << ModelName << std::endl;
  FluxPoints = ReadFlux();

  if (FluxPoints.size()==0) {
    std::cerr << "Found no flavour indices in vector returned from ReadFlux()" << std::endl;
    throw;
  }
  
  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    if (FluxPoints[iFlav].size()==0) {
      std::cerr << "No Flux points found for Model:" << ModelName << std::endl;
      throw;
    }
    
    for (size_t iP=0;iP<FluxPoints[iFlav].size();iP++) {
      if (static_cast<int>(FluxPoints[iFlav][iP].size()) != MeasDimension+1) {
	std::cerr << "Found a point with a different number of dimensions!" << std::endl;
	std::cerr << "FluxPoints[iFlav][iP].size():" << FluxPoints[iFlav][iP].size() << std::endl;
	std::cerr << "MeasDimension+1:" << MeasDimension+1 << std::endl;
	throw;
      }
    }

    if (MeasDimension == 2) {
      for (size_t iP=0;iP<FluxPoints[iFlav].size();iP++) {
	FLOAT_T Energy = FluxPoints[iFlav][iP][0];
	FLOAT_T CosineZ = FluxPoints[iFlav][iP][1];
	FLOAT_T Flux = FluxPoints[iFlav][iP][2];

	int XBin = FluxHists[iFlav]->GetXaxis()->FindBin(Energy);
	int YBin = FluxHists[iFlav]->GetYaxis()->FindBin(CosineZ);

	if (FluxHists[iFlav]->GetBinContent(XBin,YBin) != 0.) {
	  std::cerr << "Found bin that has already been set!" << std::endl;
	  throw;
	}

	FluxHists[iFlav]->SetBinContent(XBin,YBin,Flux);
      }
    } else if (MeasDimension == 3) {
      for (size_t iP=0;iP<FluxPoints[iFlav].size();iP++) {
        FLOAT_T Energy = FluxPoints[iFlav][iP][0];
        FLOAT_T CosineZ = FluxPoints[iFlav][iP][1];
	FLOAT_T Phi = FluxPoints[iFlav][iP][2];
	FLOAT_T Flux = FluxPoints[iFlav][iP][3];
	
        int XBin = FluxHists[iFlav]->GetXaxis()->FindBin(Energy);
        int YBin = FluxHists[iFlav]->GetYaxis()->FindBin(CosineZ);
	int ZBin = FluxHists[iFlav]->GetZaxis()->FindBin(Phi);

        if (FluxHists[iFlav]->GetBinContent(XBin,YBin,ZBin) != 0.) {
          std::cerr << "Found bin that has already been set!" << std::endl;
          throw;
        }

        FluxHists[iFlav]->SetBinContent(XBin,YBin,ZBin,Flux);
      }
    } else {
      std::cerr << "Unsupported number of dimensions!: " << MeasDimension << std::endl;
      throw;
    }
  }

  Build2DPlots();
  BuildFlavourRatioPlots();
  std::cout << std::endl;
}

void FluxReader::Build2DPlots() {
  if (MeasDimension == 2) {
    for (int iFlav=0;iFlav<nFlavs;iFlav++) {
      EnergyCosineZHists[iFlav] = (TH1*)(FluxHists[iFlav]->Clone());
      if (Smooth) {
	EnergyCosineZHists[iFlav]->Smooth();
      }
    }
  } else if (MeasDimension == 3) {
    for (int iFlav=0;iFlav<nFlavs;iFlav++) {
      EnergyCosineZHists[iFlav] = ((TH3*)FluxHists[iFlav])->Project3D("yx");

      FLOAT_T Max = -1;
      for (int xBin=1;xBin<=FluxHists[iFlav]->GetNbinsX();xBin++) {
        for (int yBin=1;yBin<=FluxHists[iFlav]->GetNbinsY();yBin++) {
          FLOAT_T BinContent = 0.;
          for (int zBin=1;zBin<=FluxHists[iFlav]->GetNbinsZ();zBin++) {
            BinContent += FluxHists[iFlav]->GetBinContent(xBin,yBin,zBin);
          }
          BinContent /= FluxHists[iFlav]->GetNbinsZ();
          if (BinContent > Max) {Max = BinContent;}
          EnergyCosineZHists[iFlav]->SetBinContent(xBin,yBin,BinContent);
        }
      }
      if (Smooth) {
        EnergyCosineZHists[iFlav]->Smooth();
      }
      EnergyCosineZHists[iFlav]->GetZaxis()->SetRangeUser(0,Max);
    }
  }
}

void FluxReader::BuildFlavourRatioPlots() {
  for (int iFlavRatio=0;iFlavRatio<nFlavRatios;iFlavRatio++) {
    FlavourRatioHists[iFlavRatio] = static_cast<TH1*>(EnergyCosineZHists[0]->Clone((ModelName+"_"+RatioFlavourNames[iFlavRatio]).c_str()));
    FlavourRatioHists[iFlavRatio]->SetTitle(RatioFlavourNames[iFlavRatio].c_str());
    FlavourRatioHists[iFlavRatio]->Reset();
  }

  TH2* Hist;
  for (int xBin=1;xBin<=((TH2*)FlavourRatioHists[0])->GetNbinsX();xBin++) {
    for (int yBin=1;yBin<=((TH2*)FlavourRatioHists[0])->GetNbinsY();yBin++) {
      FLOAT_T NuEBinContent = EnergyCosineZHists[NuE]->GetBinContent(xBin,yBin);
      FLOAT_T NuMBinContent = EnergyCosineZHists[NuM]->GetBinContent(xBin,yBin);
      FLOAT_T ANuEBinContent = EnergyCosineZHists[ANuE]->GetBinContent(xBin,yBin);
      FLOAT_T ANuMBinContent = EnergyCosineZHists[ANuM]->GetBinContent(xBin,yBin);
      
      FLOAT_T Ratio_T = (NuMBinContent+ANuMBinContent)/(NuEBinContent+ANuEBinContent);
      FLOAT_T Ratio_M = (NuMBinContent/ANuMBinContent);
      FLOAT_T Ratio_E = (NuEBinContent/ANuEBinContent);
      
      Hist = (TH2*)FlavourRatioHists[Ratio_Total];
      Hist->SetBinContent(xBin,yBin,Ratio_T);
      
      Hist = (TH2*)FlavourRatioHists[Ratio_NuM];
      Hist->SetBinContent(xBin,yBin,Ratio_M);

      Hist = (TH2*)FlavourRatioHists[Ratio_NuE];
      Hist->SetBinContent(xBin,yBin,Ratio_E);
    }
  }

  for (int iFlavRatio=0;iFlavRatio<nFlavRatios;iFlavRatio++) {
    Hist = (TH2*)FlavourRatioHists[iFlavRatio];
    
    FLOAT_T Min = 1e8;
    FLOAT_T Max = -1e8;
    for (int xBin=1;xBin<=Hist->GetNbinsX();xBin++) {
      for (int yBin=1;yBin<=Hist->GetNbinsY();yBin++) {
	if (Hist->GetBinContent(xBin,yBin) > Max) {Max = Hist->GetBinContent(xBin,yBin);}
	if (Hist->GetBinContent(xBin,yBin) < Min) {Min = Hist->GetBinContent(xBin,yBin);}
      }
    }
    Hist->GetZaxis()->SetRangeUser(Min,Max);
  }
}

void FluxReader::Plot2DFlux(std::string OutputName, std::string DrawOpts) {
  TCanvas* Canv = new TCanvas;
  Canv->SetLogx();
  
  Canv->SetRightMargin(0.2);
  Canv->Print((OutputName+"[").c_str());

  for (int iFlav=0;iFlav<nFlavs;iFlav++) {
    EnergyCosineZHists[iFlav]->SetTitle(FluxHists[iFlav]->GetTitle());
    EnergyCosineZHists[iFlav]->SetStats(false);      
    EnergyCosineZHists[iFlav]->GetZaxis()->SetTitle(FluxCaption.c_str());
    EnergyCosineZHists[iFlav]->GetXaxis()->SetRangeUser(EnergyAxisMin,EnergyAxisMax);
    EnergyCosineZHists[iFlav]->Draw(DrawOpts.c_str());
    EnergyCosineZHists[iFlav]->SetLineColor(LineColor);
    EnergyCosineZHists[iFlav]->SetLineStyle(LineStyle);
    
    Canv->Print(OutputName.c_str());
  }

  for (int iFlavRatio=0;iFlavRatio<nFlavRatios;iFlavRatio++) {
    FlavourRatioHists[iFlavRatio]->SetStats(false);
    FlavourRatioHists[iFlavRatio]->GetZaxis()->SetTitle(FluxCaption.c_str());
    FlavourRatioHists[iFlavRatio]->Draw(DrawOpts.c_str());
    FlavourRatioHists[iFlavRatio]->SetLineColor(LineColor);
    FlavourRatioHists[iFlavRatio]->SetLineStyle(LineStyle);

    Canv->Print(OutputName.c_str());
  }
  
  Canv->Print((OutputName+"]").c_str());
  delete Canv;
}
