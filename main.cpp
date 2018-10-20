#include <iostream>
#include "TROOT.h"
#include "include/GeneratingHits.h"
int main() {
    Int_t numOfParticles = 1500000;
    Int_t binning = 100;
    GeneratingHits *hitsGenerator = new GeneratingHits(numOfParticles, binning);

    return 0;
}