#include <iostream>
#include "TROOT.h"
#include "include/GeneratingHits.h"
#include "include/ToyAnalysis.h"

using namespace std;

int main() {
    Int_t numOfParticles = 1500000;
    Int_t binning;
    cin >> binning;
    GeneratingHits *hitsGenerator = new GeneratingHits(numOfParticles, binning);
    new ToyAnalysis(hitsGenerator->GetToyVShape(), binning);
    return 0;
}