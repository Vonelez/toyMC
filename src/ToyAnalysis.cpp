
#include "../include/ToyAnalysis.h"
ToyAnalysis::ToyAnalysis(TH2D* ToyVShape) {
  inputToyVShape = (TH2D *)ToyVShape->Clone();
  init();
  main_algorithm(inputToyVShape->FindFirstBinAbove(), inputToyVShape->FindLastBinAbove());
  finding_range();
  geometric_resol();
  writingHists();
}

ToyAnalysis::~ToyAnalysis() = default;

void ToyAnalysis::init() {
  xAxisPoint = 0.;
  yAxisPoint = 0.;
  deriv_max_l = 0.;
  deriv_max_r = 0.;

  g = new TGraphErrors();
  g->GetXaxis()->SetTitle("U (mm)");
  g->GetYaxis()->SetTitle("T (ns)");
  g->SetTitle("T = f(U)");
  Sigma = new TGraphErrors();
  Sigma->GetXaxis()->SetTitle("U (mm)");
  Sigma->GetYaxis()->SetTitle("#sigma_{T} (ns)");
  Sigma->SetTitle("#sigma_{T} = f(U)");
  U_res = new TGraphErrors();
  U_res->GetXaxis()->SetTitle("U (mm)");
  U_res->GetYaxis()->SetTitle("#sigma_{U} (mm)");
  U_res->SetTitle("#sigma_{U} = f(#sigma_{T})");
  resolgeom = new TGraphErrors();
  resolgeom->GetXaxis()->SetTitle("U (mm)");
  resolgeom->GetYaxis()->SetTitle("#sigma_{U} (mm)");
  resolgeom->SetTitle("#sigma_{U} = f(#sigma_{T})");
  straw_resol = new TGraphErrors();
  straw_resol->GetXaxis()->SetTitle("U (mm)");
  straw_resol->GetYaxis()->SetTitle("#sigma_{U} (mm)");
  straw_resol->SetTitle("#sigma_{U} = f(#sigma_{T})");
  CHI_NDF = new TGraph();
  CHI_NDF->GetXaxis()->SetTitle("U (mm)");
  CHI_NDF->GetYaxis()->SetTitle("#chi^{2} / NDF");
  NDF = new TGraph();
  NDF->GetXaxis()->SetTitle("U (mm)");
  NDF->GetYaxis()->SetTitle("NDF");
  derivative_gr = new TGraphErrors();
  derivative_gr->GetXaxis()->SetTitle("U (mm)");
  derivative_gr->GetYaxis()->SetTitle("#frac{df(U)}{dU}");
  derivative_gr->SetTitle("#frac{df(U)}{dU} = f(U)");
}

void ToyAnalysis::bins_filling(Int_t i) {
  xAxisPoint = inputToyVShape->GetXaxis()->GetBinCenter(i);
  TString name = {"bin_XAxis_"};
  TString num;
  num.Form("%f", xAxisPoint);
  TString full = name + num;
  binProjection = new TH1D();
  binProjection = inputToyVShape->ProjectionY(full, i, i+1);
  binProjection->GetXaxis()->SetTitle("T (ns)");
  binProjection->GetXaxis()->SetTitleOffset(1.4);
  binProjection->GetYaxis()->SetTitle("N entries");
}

void ToyAnalysis::finding_fit_range(Double_t axisPoint) {
  Int_t maximumBin = binProjection->GetMaximumBin();
  leftBin = 0;
  rightBin = 150;
  Int_t i = 1;
  Int_t j = 1;
  if (axisPoint >= -1. && axisPoint <= 3) {
    do {
      rightBin = maximumBin + i;
      i++;
    } while (binProjection->GetBinContent(maximumBin) * 0.2 < binProjection->GetBinContent(maximumBin + i));
  } else {
    do {
      rightBin = maximumBin + i;
      i++;
    } while (binProjection->GetBinContent(maximumBin) * 0.1 < binProjection->GetBinContent(maximumBin + i));
  }

  do {
    leftBin = maximumBin - j;
    j++;
  } while (binProjection->GetBinContent(maximumBin) * 0.1 < binProjection->GetBinContent(maximumBin - j));

  leftFitRange = binProjection->GetBinCenter(leftBin - 1);
  rightFitRange = binProjection->GetBinCenter(rightBin + 1);
}

void ToyAnalysis::fitting_bin_hist() {
  binProjection->Fit("gaus", "QR", "SAME", leftFitRange, rightFitRange);
  TF1 *fitFunction = binProjection->GetFunction("gaus");
  yAxisPoint = fitFunction->GetParameter(1);
  sigma = fitFunction->GetParameter(2);
  yAxisPointError = fitFunction->GetParError(1);
  sigmaError = fitFunction->GetParError(2);
  chi = fitFunction->GetChisquare();
  ndf = fitFunction->GetNDF();
  fit = gMinuit->fCstatu;
}

void ToyAnalysis::finding_range() {
  Int_t N = g->GetN();
  Double_t *x, *y;
  x = g->GetX();
  y = g->GetY();
  Int_t max_bin = TMath::LocMax(N, y);
  Max1 = x[max_bin];
  Max2 = 0;

  if (Max1 < 0.) {
    for (int i = max_bin; i < N; i++) {

      if (abs(y[i] / y[i+1]) > 5.) {
        Max2 = x[i];
      }
      if (abs(x[i] - x[max_bin]) > 21.) break;
    }
  } else {
    for (int i = max_bin; i < N; --i) {
      if (abs(y[i] / y[i-1]) > 5.) {
        Max2 = x[i];
        break;
      }
      if (abs(x[i] - x[max_bin]) > 21.) break;
    }
  }

  if (Max1 > Max2) {
    Double_t c = Max1;
    Max1 = Max2;
    Max2 = c;
  }
}

void ToyAnalysis::main_algorithm(Int_t start, Int_t end) {
  gStyle->SetOptFit(1111);
  first = inputToyVShape->GetXaxis()->GetBinCenter(1);
  second = inputToyVShape->GetXaxis()->GetBinCenter(2);

  xAxisPointError = (second - first) / 2.0;
  Int_t j = 0;
  for (int i = 0; i < end - start; ++i) {
    bins_filling(start + i);
    if (binProjection->Integral() < 150) continue;
    finding_fit_range(xAxisPoint);
    fitting_bin_hist();
    if (strncmp(fit.Data(), res.Data(), 5) != 0) continue;

    g->SetPoint(j, xAxisPoint, yAxisPoint);
    g->SetPointError(j, xAxisPointError, yAxisPointError);
    Sigma->SetPoint(j, xAxisPoint, sigma);
    Sigma->SetPointError(j, xAxisPointError, sigmaError);
    if (j % 10 == 0) {
      TCanvas *c1 = new TCanvas("Test", "Test", 1400, 900);
      TString folder("img");
      TString path("/test/h_x_");
      TString name = folder + path;
      Long_t num = j;
      TString end = {".root"};
      TString full = name + num + end;
      c1->cd();
      binProjection->Draw();
      c1->SaveAs(full, "Q");
    }
    if (ndf == 0) {
      CHI_NDF->SetPoint(j, xAxisPoint, 50);
      NDF->SetPoint(j, xAxisPoint, ndf);
    } else {
      Double_t chindf = chi / ndf;
      CHI_NDF->SetPoint(j, xAxisPoint, chindf);
      CHI_NDF->SetPoint(j, xAxisPoint, ndf);
      ndf_vect.push_back(ndf);
      straw_coord.push_back(xAxisPoint);
    }
    j++;
  }
}

void ToyAnalysis::geometric_resol() {
  Double_t *x = Sigma->GetX();

  for (int l = 0; l < ndf_vect.size() - 1; ++l) {
    if (deriv_max_l == 0 && straw_coord.at(l) > -12.) {
      if (abs(ndf_vect.at(l) - ndf_vect.at(l+1)) > 15) deriv_max_l = straw_coord.at(l+1);
    }
    if (abs(ndf_vect.at(l) - ndf_vect.at(l+1)) > 15 && straw_coord.at(l) > 0) deriv_max_r = straw_coord.at(l);
  }
  if (deriv_max_l == 0) deriv_max_l = straw_coord.at(0);
  if (deriv_max_r == 0) deriv_max_r = straw_coord.back();

  cout << deriv_max_l << "\n" << deriv_max_r << endl;

  deriv_calc();

  for (int k = 1; k < derivative_vec.size(); ++k) {
    if (derivative_vec.at(k-1) < 0 && derivative_vec.at(k+1) > 0  && x[k-1] > 0.0) {
      deriv_vertex = x[k];
      deriv_vertex_error = abs(x[k] - x[k+1]) / 2;
      break;
    }
  }

  straw_resolution();

  ofstream myfile;
  TString folder("img");
  TString suf("OUTPUT.txt");
  TString name = suf;
  myfile.open(name);
  myfile << "R-L \t C \t \t L \t R \n";
  myfile << deriv_max_r - deriv_max_l << "\t" << deriv_vertex << "+/-" << deriv_vertex_error << "\t" << deriv_max_l << "\t" << deriv_max_r << endl;
  myfile << "\n\nRESOLUTION " << straw_res << " Error " << straw_res_error << endl;
  myfile << "\n\nRESOLUTION check 1 " << check_1_res << " Error " << check_1_res_error << endl;
  myfile << "\n\nRESOLUTION check 2 " << check_2_res << " Error " << check_2_res_error << endl;
  myfile << "\n\nRESOLUTION check 3 " << check_3_res << " Error " << check_3_res_error << endl;
  myfile.close();
}

void ToyAnalysis::straw_resolution() {
  straw_resol = (TGraphErrors *)resolgeom->Clone();
  straw_resol->GetXaxis()->SetRangeUser(deriv_max_l, deriv_max_r);

  Double_t *x = straw_resol->GetX();

  Double_t sum = 0;
  Double_t weight_sum = 0;
  Double_t sum_1 = 0;
  Double_t weight_sum_1 = 0;
  Double_t sum_2 = 0;
  Double_t weight_sum_2 = 0;
  Double_t sum_3 = 0;
  Double_t weight_sum_3 = 0;

  Double_t *y_new = straw_resol->GetY();

  Int_t l = 0;
  Int_t r = straw_resol->GetN();
  for (int m = 0; m < straw_resol->GetN(); ++m) {
    if (x[m] <= deriv_max_l) l++;
    if (x[m] >= deriv_max_r) r--;
  }

  for (int k = l + 2; k < r - 2; ++k) {
    Double_t error = straw_resol->GetErrorY(k);
    if (y_new[k] < 0.05) continue;
    sum += y_new[k] / (error*error);
    weight_sum += 1 / (error*error);

    if (x[k] >= -6. && x[k] <= -4.) {
      sum_1 += y_new[k] / (error*error);
      weight_sum_1 += 1 / (error*error);
    }
    if (x[k] >= -1. && x[k] <= 1.) {
      sum_2 += y_new[k] / (error*error);
      weight_sum_2 += 1 / (error*error);
    }
    if (x[k] >= 4. && x[k] <= 6.) {
      sum_3 += y_new[k] / (error*error);
      weight_sum_3 += 1 / (error*error);
    }
  }

  straw_res = sum / weight_sum;
  straw_res_error = sqrt(1 / weight_sum);

  check_1_res = sum_1 / weight_sum_1;
  check_1_res_error = sqrt(1 / weight_sum_1);
  check_2_res = sum_2 / weight_sum_2;
  check_2_res_error = sqrt(1 / weight_sum_2);
  check_3_res = sum_3 / weight_sum_3;
  check_3_res_error = sqrt(1 / weight_sum_3);

  cout << straw_res << "--------- FFFFF" << endl;
}

void ToyAnalysis::deriv_calc() {
  gnew = (TGraphErrors *)g->Clone();
  Double_t *x = Sigma->GetX();
  Double_t *y = Sigma->GetY();
  Double_t *y_2 = gnew->GetY();

  Int_t Num = Sigma->GetN();

  for (int i = 1; i < Num - 1; ++i) {
    Double_t DY = y_2[i + 1] - y_2[i - 1];
    Double_t DX = x[i + 1] - x[i];
    Double_t Y_1_err = gnew->GetErrorY(i+1);
    Double_t Y_0_err = Sigma->GetErrorY(i);
    Double_t Y_2_err = gnew->GetErrorY(i-1);
    Double_t DY_error = sqrt(Y_1_err*Y_1_err + Y_2_err*Y_2_err);
    Double_t error = (1./2.) * sqrt((1 / DX)*(1 / DX)*DY_error*DY_error);
    Double_t derivative = (1./2.) * (DY / DX);
    derivative_gr->SetPoint(i-1, x[i], derivative);
    derivative_gr->SetPointError(i-1, 0., error);
    derivative_vec.push_back(derivative);
    Double_t resol = y[i] / abs(derivative);
    Double_t resol_error = (y[i] / (derivative*derivative))*(y[i] / (derivative*derivative))*error*error + (1 / derivative)*(1 / derivative)*Y_0_err*Y_0_err;

    resolgeom->SetPoint(i-1, x[i], resol);
    resolgeom->SetPointError(i-1, 0., sqrt(resol_error));
  }
}

void ToyAnalysis::writingHists() {

  TFile myfile("output.root", "RECREATE");
  g->Write("Parabola");
  gnew->Write("Parabola2");
  Sigma->Write("Sigma(x)");
  U_res->Write("RESOL");
  resolgeom->Write("RESOLgeom");
  CHI_NDF->Write("chi/ndf");
  NDF->Write("ndf");
  derivative_gr->Write("Derivative");
  straw_resol->Write("Straw_resol");
}





