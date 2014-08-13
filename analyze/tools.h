#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

namespace tools{

  float GetYield(TH1F *&hist);
  float GetYield(TH1F *&hist, float xmin, float xmax);
  double DeltaEta(const LorentzVector& p1, const LorentzVector& p2);
  void Fill1D(TH1F *&hist, double x, double w);
  void Fill2D(TH2F *&hist, double x, double y, double w);
  float GetFR(float pt, float eta, TH2F *&hist);
  float GetValue(float xvalue, float yvalue, TH2F *&hist);
  double DeltaR(const LorentzVector& p1, const LorentzVector& p2);
  double DeltaPhi(const LorentzVector& p1, const LorentzVector& p2);
  float DeltaR(float, float, float, float);
  float DeltaPhi(float, float);

}


#endif
