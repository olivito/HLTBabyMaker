// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"

// HLT
//#include "HLT.cc"
#include "HLT.h"
//#include "Include.C"
#include "tools.h"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  TH1F *h_ht   = new TH1F("h_ht", "PF H_{T}", 200,0,2000);
  TH1F *h_rate = new TH1F("h_rate", "HLT Rate (HTT200 L1 Seed, Calo Filter HT>200)", 200,0,2000);
  TH1F *h_ht_njets4   = new TH1F("h_ht_njets4", "PF H_{T}", 200,0,2000);
  TH1F *h_rate_njets4 = new TH1F("h_rate_njets4", "HLT Rate (HTT200 L1 Seed, Calo Filter HT>200, >= 4PFJets)", 200,0,2000);
  TH1F *h_nGenJets = new TH1F("h_nGenJets", "nGenJets pt > 40, |eta| < 3.0", 8, 0, 8);
  TH1F *h_nCaloJets = new TH1F("h_nCaloJets", "nCaloJets pt > 40, |eta| < 3.0", 8, 0, 8);
  TH1F *h_nPFJets = new TH1F("h_nPFJets", "nPFJets pt > 40, |eta| < 3.0", 8, 0, 8);
  TH1F *h_ht_met100 = new TH1F("h_ht_met100", "PF H_{T} MET > 100", 200,0,2000);
  TH1F *h_ht_met150 = new TH1F("h_ht_met150", "PF H_{T} MET > 150", 200,0,2000);
  TH1F *h_ht_met200 = new TH1F("h_ht_met200", "PF H_{T} MET > 200", 200,0,2000);
  TH1F *h_rate_met100 = new TH1F("h_rate_met100", "HLT Rate (HTT200 L1 Seed, Calo Filter HT>200, MET>100)", 200,0,2000);
  TH1F *h_rate_met150 = new TH1F("h_rate_met150", "HLT Rate (HTT200 L1 Seed, Calo Filter HT>200, MET>150)", 200,0,2000);
  TH1F *h_rate_met200 = new TH1F("h_rate_met200", "HLT Rate (HTT200 L1 Seed, Calo Filter HT>200, MET>200)", 200,0,2000);
  TH2F *h_rate_2d_met100 = new TH2F("h_rate_2d_met100", "h_rate_2d_met100", 24, 800, 2000, 12, 400, 1000);
  TH2F *h_rate_2d_met150 = new TH2F("h_rate_2d_met150", "h_rate_2d_met150", 24, 800, 2000, 12, 400, 1000);
  TH2F *h_rate_2d_met200 = new TH2F("h_rate_2d_met200", "h_rate_2d_met200", 24, 800, 2000, 12, 400, 1000);

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    hlt.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      hlt.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      HLT::progress( nEventsTotal, nEventsChain );

      // Analysis Code
      if(hlt.pf_ht() < 1.0) continue;

      float weight = hlt.scale1fb()*(0.5*TMath::Erf((1.35121e-02)*(hlt.pf_ht()-(3.02695e+02)))+0.5);
      tools::Fill1D(h_ht, hlt.pf_ht(), weight);

    
      //jets
      int njets = 0;
      for(unsigned int j=0; j< hlt.genjets_pt().size(); j++){
        if(hlt.genjets_pt().at(j) < 40.0) continue;
        if(fabs(hlt.genjets_eta().at(j)) > 3.0) continue;
        njets++;
      }
      tools::Fill1D(h_nGenJets, njets, hlt.scale1fb());

      njets = 0;
      for(unsigned int j=0; j< hlt.calojets_pt().size(); j++){
        if(hlt.calojets_pt().at(j) < 40.0) continue;
        if(fabs(hlt.calojets_eta().at(j)) > 3.0) continue;
        njets++;
      }
      tools::Fill1D(h_nCaloJets, njets, hlt.scale1fb());
      
      njets = 0;
      for(unsigned int j=0; j< hlt.pfjets_pt().size(); j++){
        if(hlt.pfjets_pt().at(j) < 40.0) continue;
        if(fabs(hlt.pfjets_eta().at(j)) > 3.0) continue;
        njets++;
      }
      tools::Fill1D(h_nPFJets, njets, hlt.scale1fb());


      if(njets >= 4) tools::Fill1D(h_ht_njets4, pf_ht(), weight);

      
      if(hlt.met_pt() < 100) continue;
      tools::Fill1D(h_ht_met100, hlt.pf_ht(), weight);

      if(hlt.met_pt() < 150) continue;
      tools::Fill1D(h_ht_met150, hlt.pf_ht(), weight);
      
      if(hlt.met_pt() < 200) continue;
      tools::Fill1D(h_ht_met200, hlt.pf_ht(), weight);
      

    }
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }



  //float lumi = 1.4e34 cm^-2 s^-1
  float lumi = 1.4e-5; //fb^-1 s^-1
  float rate = -999.9;
  float rate_njets4 = -999.9;
  float rate_met100 = -999.9;
  float rate_met150 = -999.9;
  float rate_met200 = -999.9;
  for(int i=0; i<h_ht->GetNbinsX(); i++){
    rate = tools::GetYield(h_ht, 10*i, 4000)*lumi;
    rate_njets4 = tools::GetYield(h_ht_njets4, 10*i, 4000)*lumi;
    rate_met100 = tools::GetYield(h_ht_met100, 10*i, 4000)*lumi;
    rate_met150 = tools::GetYield(h_ht_met150, 10*i, 4000)*lumi;
    rate_met200 = tools::GetYield(h_ht_met200, 10*i, 4000)*lumi;
    tools::Fill1D(h_rate, 10*i, rate);
    tools::Fill1D(h_rate_njets4, 10*i, rate_njets4);
    tools::Fill1D(h_rate_met100, 10*i, rate_met100);
    tools::Fill1D(h_rate_met150, 10*i, rate_met150);
    tools::Fill1D(h_rate_met200, 10*i, rate_met200);
  }



  for(int i=0; i<h_rate_2d_met100->GetNbinsX(); i++){
    for(int j=0; j<h_rate_2d_met100->GetNbinsY(); j++){
      float rate1 = h_rate->GetBinContent(h_rate->FindBin(50*i + 800));
      float rate2 = h_rate_met100->GetBinContent(h_rate_met100->FindBin(50*j + 400));
      tools::Fill2D(h_rate_2d_met100, 50*i + 800, 50*j + 400, rate1+rate2);
    }
  }

  for(int i=0; i<h_rate_2d_met150->GetNbinsX(); i++){
    for(int j=0; j<h_rate_2d_met150->GetNbinsY(); j++){
      float rate1 = h_rate->GetBinContent(h_rate->FindBin(50*i + 800));
      float rate2 = h_rate_met150->GetBinContent(h_rate_met150->FindBin(50*j + 400));
      tools::Fill2D(h_rate_2d_met150, 50*i + 800, 50*j + 400, rate1+rate2);
    }
  }

  for(int i=0; i<h_rate_2d_met200->GetNbinsX(); i++){
    for(int j=0; j<h_rate_2d_met200->GetNbinsY(); j++){
      float rate1 = h_rate->GetBinContent(h_rate->FindBin(50*i + 800));
      float rate2 = h_rate_met200->GetBinContent(h_rate_met200->FindBin(50*j + 400));
      tools::Fill2D(h_rate_2d_met200, 50*i + 800, 50*j + 400, rate1+rate2);
    }
  }



  TCanvas* c0 = new TCanvas;

  h_nGenJets->Draw();
  float normalizer = 1.0/tools::GetYield(h_nGenJets);
  h_nGenJets->Scale(normalizer);
  c0->Print("h_nGenJets.pdf");
  h_nCaloJets->Draw();
  normalizer = 1.0/tools::GetYield(h_nCaloJets);
  h_nCaloJets->Scale(normalizer);
  c0->Print("h_nCaloJets.pdf");
  h_nPFJets->Draw();
  normalizer = 1.0/tools::GetYield(h_nPFJets);
  h_nPFJets->Scale(normalizer);
  c0->Print("h_nPFJets.pdf");


  h_rate_2d_met100->GetXaxis()->SetTitle("High H_{T} Threshold[GeV]");
  h_rate_2d_met100->GetYaxis()->SetTitle("Low H_{T} Threshold[GeV]");
  h_rate_2d_met100->Draw("COLZ");
  c0->Print("h_rate_2d_met100.pdf");

  h_rate_2d_met150->GetXaxis()->SetTitle("High H_{T} Threshold[GeV]");
  h_rate_2d_met150->GetYaxis()->SetTitle("Low H_{T} Threshold[GeV]");
  h_rate_2d_met150->Draw("COLZ");
  c0->Print("h_rate_2d_met150.pdf");

  h_rate_2d_met200->GetXaxis()->SetTitle("High H_{T} Threshold[GeV]");
  h_rate_2d_met200->GetYaxis()->SetTitle("Low H_{T} Threshold[GeV]");
  h_rate_2d_met200->Draw("COLZ");
  c0->Print("h_rate_2d_met200.pdf");


  h_ht_met100->Draw();
  c0->Print("h_ht_met100.pdf");
  h_ht_met150->Draw();
  c0->Print("h_ht_met150.pdf");
  h_ht_met200->Draw();
  c0->Print("h_ht_met200.pdf");

  h_rate_met100->Draw();
  c0->Print("h_rate_met100.pdf");
  h_rate_met150->Draw();
  c0->Print("h_rate_met150.pdf");
  h_rate_met200->Draw();
  c0->Print("h_rate_met200.pdf");


  delete c0;

  gStyle->SetOptStat("");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameBorderMode(0);

  TCanvas* c1 = new TCanvas;

  h_ht->GetXaxis()->SetTitle("HLT PF H_{T} [GeV]");
  h_ht->GetYaxis()->SetTitle("Events/10 GeV");
  h_ht->Draw();
  c1->Print("h_ht.pdf");


  h_rate->GetXaxis()->SetTitle("HLT PF H_{T} [GeV]");
  h_rate->GetYaxis()->SetTitle("HLT Rate [Hz]");
  h_rate->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  //t->SetTextFont(42);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->DrawLatex(0.68, 0.80, "L_{inst} = 1.4E34 cm^{-2} s^{-1}");
  t->DrawLatex(0.68, 0.75, "#sqrt{s} = 13 TeV, PU40bx25");
  c1->SetLogy();
  c1->Print("h_rate.pdf");



  h_rate_njets4->GetXaxis()->SetTitle("HLT PF H_{T} [GeV]");
  h_rate_njets4->GetYaxis()->SetTitle("HLT Rate [Hz]");
  h_rate_njets4->Draw();

  t->SetNDC();
  //t->SetTextFont(42);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->DrawLatex(0.68, 0.80, "L_{inst} = 1.4E34 cm^{-2} s^{-1}");
  t->DrawLatex(0.68, 0.75, "#sqrt{s} = 13 TeV, PU40bx25");
  c1->SetLogy();
  c1->Print("h_rate_njets4.pdf");

  
  delete c1;
  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
