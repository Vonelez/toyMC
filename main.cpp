#include <iostream>
#include "TROOT.h"
#include "include/GeneratingHits.h"
#include "include/ToyAnalysis.h"

using namespace std;

int main() {
    Int_t numOfParticles = 1500000;
    Int_t binning = 100;
    GeneratingHits *hitsGenerator = new GeneratingHits(numOfParticles, binning);
    ToyAnalysis *ana = new ToyAnalysis(hitsGenerator->GetToyVShape());


    return 0;
}