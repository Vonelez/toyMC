#include "../include/GeneratingHits.h"



GeneratingHits::GeneratingHits(Int_t numOfHits, Int_t binning) {
  init();
  randomizeHits(numOfHits);
  fillingToyVShape(binning, numOfHits, uAxisPoints, tAxisPoints);
}

GeneratingHits::~GeneratingHits() = default;

void GeneratingHits::init() {
  seed = 364;
  rand = new TRandom3(seed);

  parabolaFunction = new TF1("parabola", "[0]+[1]*x+[2]*x*x");
  parabolaFunction->SetParameters(-10.08, -19.71, 7.342);
}

void GeneratingHits::randomizeHits(Int_t numOfHits) {
  uAxisPoints = new Double_t[numOfHits];
  tAxisPoints = new Double_t[numOfHits];

  for (int i = 0; i < numOfHits; ++i) {
    uAxisPoints[i] = rand->Uniform(-10.0, 10.0);
    tAxisPoints[i] = parabolaFunction->Eval(uAxisPoints[i]);
  }

}

void GeneratingHits::fillingToyVShape(Int_t binning, Int_t arrayLength, Double_t x[], Double_t y[]) {
  Int_t bins = 22 * 1000 / binning;
  toyVShape = new TH2D("toyVShape", "toyVShape", bins, -11.0, 11.0, 1100, -100.0, 1000.0);
  toyVShape->GetXaxis()->SetTitle("U (mm)");
  toyVShape->GetYaxis()->SetTitle("T (ns)");

  for (int i = 0; i < arrayLength; ++i) {
    toyVShape->Fill(x[i], y[i]);
  }

  TCanvas *shape = new TCanvas("toyVShape", "Shape", 1400, 1000);
  shape->cd();
  toyVShape->Draw("COLZ");
  shape->SaveAs("img/toyVShape.pdf", "Q");
}

TH2D * GeneratingHits::GetToyVShape() {
  return toyVShape;
}

