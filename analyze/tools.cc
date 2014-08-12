#include "tools.h"

using namespace std;

namespace tools{

  float GetYield(TH1F *&hist){
    return hist->Integral( hist->GetXaxis()->GetFirst(), hist->GetXaxis()->GetLast() );
  }

  float GetYield(TH1F *&hist, float xmin, float xmax){
    int bmin = hist->GetXaxis()->FindBin(xmin);
    int bmax = hist->GetXaxis()->FindBin(xmax);
    return hist->Integral(bmin, bmax);
  }

  double DeltaEta(const LorentzVector& p1, const LorentzVector& p2)
  {
      return fabs(p1.Eta()-p2.Eta());
  }

  void Fill1D(TH1F *&hist, double x, double w)
  {   
      x = min(hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->GetLast()) , x);
      x = max(hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->GetFirst()), x);
      hist->Fill(x, w);
  }

  void Fill2D(TH2F *&hist, double x, double y, double w)
  {   
      x = min(hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->GetLast()) , x); 
      x = max(hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->GetFirst()), x); 
      y = min(hist->GetYaxis()->GetBinCenter(hist->GetYaxis()->GetLast()) , y); 
      y = max(hist->GetYaxis()->GetBinCenter(hist->GetYaxis()->GetFirst()), y); 
      hist->Fill(x, y, w); 
  }

  float GetFR(float pt, float eta, TH2F *&hist){
    float max_pt = hist->GetYaxis()->GetXmax()-0.01;
    int pt_bin   = hist->GetYaxis()->FindBin(min(pt, max_pt));
    int eta_bin  = hist->GetXaxis()->FindBin(fabs(eta));
    return hist->GetBinContent(eta_bin, pt_bin);
  }

  float GetValue(float xvalue, float yvalue, TH2F *&hist){
    float xmax = hist->GetXaxis()->GetXmax()-0.01;
    int xbin   = hist->GetXaxis()->FindBin(min(xvalue, xmax));
    //float ymax = hist->GetYaxis()->GetYmax()-0.01;
    int ybin   = hist->GetYaxis()->FindBin(yvalue);
    return hist->GetBinContent(xbin, ybin);
  }

  double DeltaR(const LorentzVector& p1, const LorentzVector& p2)
  { 
      return ROOT::Math::VectorUtil::DeltaR(p1, p2); 
      //float deta = p1.Eta()-p2.Eta();
      //float dphi = DeltaPhi(p1.Phi(), p2.Phi());
      //return sqrt(dphi*dphi + deta*deta);
  }

  double DeltaPhi(const LorentzVector& p1, const LorentzVector& p2)
  {
      return ROOT::Math::VectorUtil::DeltaPhi(p1, p2); 
  }


}
