#ifndef TOYMC_TOYANALYSIS_H
#define TOYMC_TOYANALYSIS_H

#include "TH2D.h"
#include "vector"
#include "TMinuit.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TFile.h"
using namespace std;

class ToyAnalysis {
 public:
  ToyAnalysis(TH2D* ToyVShape);
  virtual ~ToyAnalysis();

 private:
  TH2D* inputToyVShape;

  Double_t xAxisPoint;
  Double_t yAxisPoint;
  Double_t xAxisPointError;
  Double_t yAxisPointError;
  Double_t sigma;
  Double_t sigmaError;
  Int_t leftBin;
  Int_t rightBin;
  Double_t Max1, Max2;

  Double_t second;
  Double_t first;

  TH1D *binProjection;
  Double_t leftFitRange;
  Double_t rightFitRange;

  TGraphErrors *g;
  TGraphErrors *gnew;
  TGraphErrors *Sigma;
  TGraphErrors *U_res;
  TGraphErrors *resolgeom;
  TGraph *CHI_NDF;
  TGraph *NDF;
  TGraphErrors *derivative_gr;
  TGraphErrors *straw_resol;

  Double_t chi, ndf;
  vector<Double_t > derivative_vec;
  Double_t deriv_vertex, deriv_vertex_error;
  vector<Double_t > ndf_vect;
  vector<Double_t > straw_coord;
  Double_t deriv_max_l, deriv_max_r;
  Double_t straw_res, straw_res_error;
  Double_t check_1_res, check_1_res_error;
  Double_t check_2_res, check_2_res_error;
  Double_t check_3_res, check_3_res_error;

  TString fit;
  TString res = {"CONVERGED"};

  virtual void init();
  virtual void bins_filling(Int_t i);
  virtual void finding_fit_range(Double_t axisPoint);
  virtual void fitting_bin_hist();
  virtual void finding_range();
  virtual void main_algorithm(Int_t start, Int_t end);
  virtual void geometric_resol();
  virtual void straw_resolution();
  virtual void deriv_calc();
  virtual void writingHists();
};

#endif //TOYMC_TOYANALYSIS_H
