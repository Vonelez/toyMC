#ifndef TOYMC_GENERATINGHITS_H
#define TOYMC_GENERATINGHITS_H

#include "TH2D.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"

class GeneratingHits {
 public:
  TH2D *toyVShape;
  GeneratingHits(Int_t numOfHits, Int_t binning);
  virtual TH2D * GetToyVShape();

 private:
  UInt_t seed;
  TRandom *rand;
  TF1 *parabolaFunction;
  Double_t *uAxisPoints;
  Double_t *tAxisPoints;

  virtual ~GeneratingHits();
  virtual void init();
  virtual void randomizeHits(Int_t numOfHits);
  virtual void fillingToyVShape(Int_t binning, Int_t arrayLength, Double_t x[], Double_t y[]);
};

#endif //TOYMC_GENERATINGHITS_H
