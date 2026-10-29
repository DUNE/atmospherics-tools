#pragma once

#include "TH2D.h"
#include "TSpline.h"
#include "TGraph.h"

#include <sstream>
#include <string>

inline std::vector<std::string> SplitLine(const std::string& str) {
  std::vector<std::string> tokens;
  std::istringstream iss(str);
  std::string token;

  while (iss >> token) {
    tokens.push_back(token);
  }

  return tokens;
}

inline TH2D* InterpolateHistogram(TH2* ModelHistogram, TH2* TemplateHistogram,  std::string HistName="") {
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

inline TGraph* TSpline3_to_TGraph(TSpline3* spline){
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
