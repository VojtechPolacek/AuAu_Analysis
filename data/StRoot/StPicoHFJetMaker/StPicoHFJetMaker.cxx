#include "StPicoHFJetMaker.h"
#include "JetInfo.h"

#include "BemcNewCalib.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/defines.h"

#include "TRandom3.h"
#include "TVector2.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/config.h"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "MyJet.h"

using namespace std;

const char* kCentTag[4] = { // centrality classes
  "",            // index 0 unused
  "CENT_0_10",   // c3 = 1
  "MID_20_40",   // c3 = 2
  "PERI_60_80"   // c3 = 3
};

const std::vector<TString> StPicoHFJetMaker::fConfigNames = {"full", "charged"}; // configurations for jet finding and background estimation

std::vector<float> phi_vector;
std::vector<float> eta_vector;
std::vector<float> E_vector;
std::vector<float> pt_vector;
float c_px, c_py, c_pz, c_pt, c_E, c_phi, c_eta;
int c_charge;
int c_runid, c_eventid, c_ijet;
int eec_ijet, eec_eventid, eec_runid;
float eec_data, RL;
//float EEC_low, EEC_mid, EEC_high;

 

//MC constituents and EEC variables
std::vector<float> mc_phi_vector;
std::vector<float> mc_eta_vector;
std::vector<float> mc_E_vector;
std::vector<float> mc_pt_vector;
float c_mc_px, c_mc_py, c_mc_pz, c_mc_pt, c_mc_E, c_mc_phi, c_mc_eta;
int c_mc_charge;
int c_mc_runid, c_mc_eventid, c_mc_ijet;
int eec_mc_ijet, eec_mc_eventid, eec_mc_runid;
float eec_mc, RL_mc;

//Matched constituents and EEC variables
std::vector<float> matched_phi_vector;
std::vector<float> matched_eta_vector;
std::vector<float> matched_E_vector;
std::vector<float> matched_pt_vector;
float c_matched_px, c_matched_py, c_matched_pz, c_matched_pt, c_matched_E, c_matched_phi, c_matched_eta;
int c_matched_charge;
int c_matched_runid, c_matched_eventid, c_matched_ijet;
int eec_matched_ijet, eec_matched_eventid, eec_matched_runid;
float eec_matched, RL_matched;

//Unmatched constituents and EEC variables
std::vector<float> unmatched_phi_vector;
std::vector<float> unmatched_eta_vector;
std::vector<float> unmatched_E_vector;
std::vector<float> unmatched_pt_vector;
float c_unmatched_px, c_unmatched_py, c_unmatched_pz, c_unmatched_pt, c_unmatched_E, c_unmatched_phi, c_unmatched_eta;
int c_unmatched_charge;
int c_unmatched_runid, c_unmatched_eventid, c_unmatched_ijet;
int eec_unmatched_ijet, eec_unmatched_eventid, eec_unmatched_runid;
float eec_unmatched, RL_unmatched;

//---------------------Binning for EEC------------------------
//Limits for R=0,2
int N_bins_02 = 33;
double EEC_bounds_02[34] = {
  0.0, 0.00500, 0.00571, 0.00653, 0.00746, 0.00853, 0.00975, 0.01114, 0.01273, 0.01455,
  0.01663, 0.01901, 0.02172, 0.02482, 0.02837, 0.03242, 0.03706, 0.04235, 0.04840,
  0.05531, 0.06321, 0.07225, 0.08257, 0.09436, 0.10784, 0.12325, 0.14085, 0.16098,
  0.18397, 0.21000, 0.25000, 0.29748, 0.35373, 0.4};

//Limits for R=0,3
int N_bins_03 = 37;
double EEC_bounds_03[38] = {0.0, 0.005, 0.00571, 0.00653, 0.00746, 0.00853, 0.00975, 0.01114, 0.01273, 0.01455, 0.01663, 0.01901, 0.02172, 0.02482, 0.02837, 0.03242, 0.03706, 0.04235, 0.04840, 0.05531, 0.06321, 0.07225, 0.08257, 0.09436, 0.10784, 0.12325, 0.14085, 0.16098, 0.18397, 0.21025, 0.24029, 0.27462, 0.31385, 0.35868, 0.40992, 0.46849, 0.53541,0.6};
  
//Limits for R=0,4
int N_bins_04 = 37;
double EEC_bounds_04[38] = {
  0.0, 0.00500, 0.00571, 0.00653, 0.00746, 0.00853, 0.00975, 0.01114, 0.01273, 0.01455,
  0.01663, 0.01901, 0.02172, 0.02482, 0.02837, 0.03242, 0.03706, 0.04235, 0.04840,
  0.05531, 0.06321, 0.07225, 0.08257, 0.09436, 0.10784, 0.12325, 0.14085, 0.16098,
  0.18397, 0.21000, 0.25000, 0.29748, 0.35373, 0.42030, 0.49905, 0.59239, 0.70321, 0.8
};


const double CUT_AREA_02 = 0.07; // R = 0.2
const double CUT_AREA_03 = 0.20; // R = 0.3
const double CUT_AREA_04 = 0.40; // R = 0.4

const double CUT_NEUTRAL_FRACTION = 0.95; 



vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R); // match reco jets to MC jets based on eta-phi distance; return pairs of matched jets (deltaR < R)

ClassImp(StPicoHFJetMaker)

StPicoHFJetMaker::StPicoHFJetMaker(TString name, StPicoDstMaker *picoMaker,
                                   TString outputBaseFileName)
    : StPicoJetMaker(name, picoMaker, outputBaseFileName)
{
  mRefmultCorrUtil   = NULL;
  fMcSumPt           = 0.0f;
  fEventId = 0;
}


// _________________________________________________________
StPicoHFJetMaker::~StPicoHFJetMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFJetMaker::InitJets() {
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  TH1::SetDefaultSumw2();


  //---------------------QA histograms------------------------
  mOutList->SetName("QA_histograms"); // list to hold QA histograms; will be written to file

  TH1D* hcent9 = new TH1D("hcent9", "centrality9;bin;events", 10, -1, 9);
  mOutList->Add(hcent9); // centrality distribution with 9 bins (0-5%, 5-10%, ..., 70-80%) - adjust as needed for your dataset

  TAxis* ax = hcent9->GetXaxis();
  const char* lab9[10] = {
    "undef (-1)",
    "0: 0-5%",
    "1: 5-10%",
    "2: 10-20%",
    "3: 20-30%",
    "4: 30-40%",
    "5: 40-50%",
    "6: 50-60%",
    "7: 60-70%",
    "8: 70-80%"
  }; // centrality bin labels - adjust as needed for your dataset

  int b;
  for (b = 1; b <= 10; ++b) ax->SetBinLabel(b, lab9[b-1]);
  ax->CenterLabels(true);
  mOutList->Add(new TH1D("hTowE", "BEMC tower E;E (GeV);counts", 400, 0.0, 40.0));
  mOutList->Add(new TH2D("hTowEtaPhi", "BEMC towers;#eta;#phi",
                       40, -1.0, 1.0, 120, 0.0, 2.0*TMath::Pi()));

  mOutList->Add(new TH1D("hTrackPt", "Primary tracks;p_{T} (GeV/c);counts", 400, 0.0, 40.0));
  mOutList->Add(new TH2D("hTrackEtaPhi", "Primary tracks;#eta;#phi",
                       40, -1.0, 1.0, 120, 0.0, 2.0*TMath::Pi()));
  mOutList->Add(new TH2D("hRhoVsRefMultfull",
                       "#rho vs refMult;refMult;#rho_full (GeV/c per unit area)",
                       800, 0, 800,     // adjust refMult range/binning to your dataset
                       200, 0.0, 200.0  // adjust rho range/binning as needed
                       ));
                       
  mOutList->Add(new TH2D("hRhoVsRefMultcharged",
                       "#rho vs refMult;refMult;#rho_charged (GeV/c per unit area)",
                       800, 0, 800,     // adjust refMult range/binning to your dataset
                       200, 0.0, 200.0  // adjust rho range/binning as needed
                       ));
  

  //---------------------------Vectors of TTrees and histograms for different R and centrality classes---------------------------
  //JetTrees
  TDirectory* fileDir = gDirectory;
  fTreeRC.clear();
  fTreeRC.reserve(fR.size());

  MCJetTreeRC.clear();
  MCJetTreeRC.reserve(fR.size());
  //Constituent trees and EEC trees for data
  fConstituentTreeRC.clear();
  fConstituentTreeRC.reserve(fR.size());
  fEECTreeRC.clear();
  fEECTreeRC.reserve(fR.size());
  
  
  //EEC histograms for data
  //fHistEEC.clear();
  //fHistEEC.reserve(fR.size());
  fHistEEC_5_10.clear();
  fHistEEC_5_10.reserve(fR.size());
  fHistEEC_10_15.clear();
  fHistEEC_10_15.reserve(fR.size());
  fHistEEC_15_20.clear();
  fHistEEC_15_20.reserve(fR.size());
  fHistEEC_20_30.clear();
  fHistEEC_20_30.reserve(fR.size());
  fHistEEC_30_50.clear();
  fHistEEC_30_50.reserve(fR.size());
  

  //Trees for embedding: matched jets, unmatched jets, and MC jets
  fMatchedJetConstituentTreeRC.clear();
  fMatchedJetConstituentTreeRC.reserve(fR.size());
  fUnmatchedJetConstituentTreeRC.clear();
  fUnmatchedJetConstituentTreeRC.reserve(fR.size());
  fMCJetConstituentTreeRC.clear();
  fMCJetConstituentTreeRC.reserve(fR.size());
  fEECTreeunmatchedRC.clear();
  fEECTreeunmatchedRC.reserve(fR.size());
  fEECTreematchedRC.clear();
  fEECTreematchedRC.reserve(fR.size());
  fEECTree_MC_RC.clear();
  fEECTree_MC_RC.reserve(fR.size());

  //Histograms for embedding: matched jets, unmatched jets, and MC jets
  fHistEEC_MC_5_10.clear();
  fHistEEC_MC_5_10.reserve(fR.size());
  fHistEEC_MC_10_15.clear();
  fHistEEC_MC_10_15.reserve(fR.size());
  fHistEEC_MC_15_20.clear();
  fHistEEC_MC_15_20.reserve(fR.size());
  fHistEEC_MC_20_30.clear();
  fHistEEC_MC_20_30.reserve(fR.size());
  fHistEEC_MC_30_50.clear();
  fHistEEC_MC_30_50.reserve(fR.size());
  fHistEEC_MC_50_100.clear();
  fHistEEC_MC_50_100.reserve(fR.size());

  fHistEEC_matched_5_10.clear();
  fHistEEC_matched_5_10.reserve(fR.size());
  fHistEEC_matched_10_15.clear();
  fHistEEC_matched_10_15.reserve(fR.size());
  fHistEEC_matched_15_20.clear();
  fHistEEC_matched_15_20.reserve(fR.size());
  fHistEEC_matched_20_30.clear();
  fHistEEC_matched_20_30.reserve(fR.size());
  fHistEEC_matched_30_50.clear();
  fHistEEC_matched_30_50.reserve(fR.size());
  fHistEEC_matched_50_100.clear();
  fHistEEC_matched_50_100.reserve(fR.size());

  fHistEEC_unmatched_5_10.clear();
  fHistEEC_unmatched_5_10.reserve(fR.size());
  fHistEEC_unmatched_10_15.clear();
  fHistEEC_unmatched_10_15.reserve(fR.size());
  fHistEEC_unmatched_15_20.clear();
  fHistEEC_unmatched_15_20.reserve(fR.size());
  fHistEEC_unmatched_20_30.clear();
  fHistEEC_unmatched_20_30.reserve(fR.size());
  fHistEEC_unmatched_30_50.clear();
  fHistEEC_unmatched_30_50.reserve(fR.size());
  fHistEEC_unmatched_50_100.clear();
  fHistEEC_unmatched_50_100.reserve(fR.size());
  
  fHistEEC_unmatched_all.clear();
  fHistEEC_unmatched_all.reserve(fR.size());

  for (size_t iR = 0; iR < fR.size(); ++iR) { // loop over R values
    const TString rName = Form("R%.1f", fR[iR]); // directory name for this R
    TDirectory* rdir = fileDir->mkdir(rName); // create directory for this R
    if (!rdir) rdir = (TDirectory*)fileDir->Get(rName); // if it already exists (e.g. from a previous run), get it
    rdir->cd(); // change to this directory for creating subdirectories and objects

    int nBinsEEC = 0; // number of bins for EEC histograms for this R
    double* EEC_bounds = nullptr; // pointer to bin edges for EEC histograms for this R

    if(std::abs(fR[iR] - 0.2) < 1e-3) { //Choose EEC binning based on R value; adjust the tolerance as needed for floating-point comparisons
      nBinsEEC = N_bins_02;
      EEC_bounds = EEC_bounds_02;
    } else if(std::abs(fR[iR] - 0.3) < 1e-3) {
      nBinsEEC = N_bins_03;
      EEC_bounds = EEC_bounds_03;
    } else if(std::abs(fR[iR] - 0.4) < 1e-3) {
      nBinsEEC = N_bins_04;
      EEC_bounds = EEC_bounds_04;
    }
   
    
    // create vectors to hold TTrees and histograms for this R; we'll fill these and then add them to the main vectors after the loop over centrality classes
    std::vector<std::vector<TTree*>> treesConfig;
    treesConfig.reserve(2);


    std::vector<std::vector<TTree*>> MCJetTreeConfig;
    MCJetTreeConfig.reserve(2);

    std::vector<std::vector<TTree*>> ConstituentTreeConfig;
    ConstituentTreeConfig.reserve(2);
    std::vector<std::vector<TTree*>> EECTreeConfig;
    EECTreeConfig.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_5_10_Config;
    HistEEC_5_10_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_10_15_Config;
    HistEEC_10_15_Config.reserve(2);    
    std::vector<std::vector<TH1D*>> HistEEC_15_20_Config;
    HistEEC_15_20_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_20_30_Config;
    HistEEC_20_30_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_30_50_Config;
    HistEEC_30_50_Config.reserve(2);
    
    std::vector<std::vector<TTree*>> EECTreematchedConfig;
    EECTreematchedConfig.reserve(2);
    std::vector<std::vector<TTree*>> EECTreeunmatchedConfig;
    EECTreeunmatchedConfig.reserve(2);
    std::vector<std::vector<TTree*>> EECTree_MC_Config;
    EECTree_MC_Config.reserve(2);
    std::vector<std::vector<TTree*>> MatchedJetConstituentTreeConfig;
    MatchedJetConstituentTreeConfig.reserve(2);
    std::vector<std::vector<TTree*>> UnmatchedJetConstituentTreeConfig;
    UnmatchedJetConstituentTreeConfig.reserve(2);
    std::vector<std::vector<TTree*>> MCJetConstituentTreeConfig;
    MCJetConstituentTreeConfig.reserve(2);
    
    std::vector<std::vector<TH1D*>> HistEEC_MC_5_10_Config;
    HistEEC_MC_5_10_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_MC_10_15_Config;
    HistEEC_MC_10_15_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_MC_15_20_Config;
    HistEEC_MC_15_20_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_MC_20_30_Config;
    HistEEC_MC_20_30_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_MC_30_50_Config;
    HistEEC_MC_30_50_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_MC_50_100_Config;
    HistEEC_MC_50_100_Config.reserve(2);

    std::vector<std::vector<TH1D*>> HistEEC_matched_5_10_Config;
    HistEEC_matched_5_10_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_matched_10_15_Config;
    HistEEC_matched_10_15_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_matched_15_20_Config;
    HistEEC_matched_15_20_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_matched_20_30_Config;
    HistEEC_matched_20_30_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_matched_30_50_Config;
    HistEEC_matched_30_50_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_matched_50_100_Config;
    HistEEC_matched_50_100_Config.reserve(2);

    std::vector<std::vector<TH1D*>> HistEEC_unmatched_5_10_Config;
    HistEEC_unmatched_5_10_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_unmatched_10_15_Config;
    HistEEC_unmatched_10_15_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_unmatched_15_20_Config;
    HistEEC_unmatched_15_20_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_unmatched_20_30_Config;
    HistEEC_unmatched_20_30_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_unmatched_30_50_Config;
    HistEEC_unmatched_30_50_Config.reserve(2);
    std::vector<std::vector<TH1D*>> HistEEC_unmatched_50_100_Config;
    HistEEC_unmatched_50_100_Config.reserve(2);

    std::vector<std::vector<TH1D*>> HistEEC_unmatched_all_Config;
    HistEEC_unmatched_all_Config.reserve(2);

    for (size_t iConfig = 0; iConfig < fConfigNames.size(); ++iConfig){ // loop over configurations (full, charged)
      TDirectory* configDir = rdir->mkdir(fConfigNames[iConfig]); // create directory for this configuration
      if (!configDir) configDir = (TDirectory*)rdir->Get(fConfigNames[iConfig]); // if it already exists, get it
      configDir->cd();

      //Create vectors to hold TTrees and histograms for this configuration; we'll fill these and then add them to the main vectors after the loop over centrality classes
      std::vector<TTree*> treesC;  // 3 classes: 1..3 (we'll index 0..2)
      treesC.reserve(3);

      std::vector<TTree*> MCJetTreeC;  // 3 classes: 1..3 (we'll index 0..2)
      MCJetTreeC.reserve(3);

      std::vector<TTree*> ConstituentTreeC; // member of StPicoHFJetMaker
      ConstituentTreeC.reserve(3);
      std::vector<TTree*> EECTreeC; // member of StPicoHFJetMaker
      EECTreeC.reserve(3);
      //std::vector<TH1D*> Hist_EEC_C; // member of StPicoHFJetMaker
      //Hist_EEC_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_C_5_10; // member of StPicoHFJetMaker
      Hist_EEC_C_5_10.reserve(3);
      std::vector<TH1D*> Hist_EEC_C_10_15; // member of StPicoHFJetMaker
      Hist_EEC_C_10_15.reserve(3);
      std::vector<TH1D*> Hist_EEC_C_15_20; // member of StPicoHFJetMaker
      Hist_EEC_C_15_20.reserve(3);
      std::vector<TH1D*> Hist_EEC_C_20_30; // member of StPicoHFJetMaker
      Hist_EEC_C_20_30.reserve(3);
      std::vector<TH1D*> Hist_EEC_C_30_50; // member of StPicoHFJetMaker
      Hist_EEC_C_30_50.reserve(3);
      

      
      std::vector<TTree*> EECTreematchedC; // member of StPicoHFJetMaker
      EECTreematchedC.reserve(3);
      std::vector<TTree*> EECTreeunmatchedC; // member of StPicoHFJetMaker
      EECTreeunmatchedC.reserve(3);
      std::vector<TTree*> EECTree_MC_C; // member of StPicoHFJetMaker
      EECTree_MC_C.reserve(3);
      std::vector<TTree*> MatchedJetConstituentTreeC; // member of StPicoHFJetMaker
      MatchedJetConstituentTreeC.reserve(3);
      std::vector<TTree*> UnmatchedJetConstituentTreeC; // member of StPicoHFJetMaker
      UnmatchedJetConstituentTreeC.reserve(3);
      std::vector<TTree*> MCJetConstituentTreeC; // member of StPicoHFJetMaker
      MCJetConstituentTreeC.reserve(3);

      std::vector<TH1D*> Hist_EEC_MC_5_10_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_5_10_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_MC_10_15_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_10_15_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_MC_15_20_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_15_20_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_MC_20_30_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_20_30_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_MC_30_50_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_30_50_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_MC_50_100_C; // member of StPicoHFJetMaker
      Hist_EEC_MC_50_100_C.reserve(3);

      std::vector<TH1D*> Hist_EEC_matched_5_10_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_5_10_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_matched_10_15_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_10_15_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_matched_15_20_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_15_20_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_matched_20_30_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_20_30_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_matched_30_50_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_30_50_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_matched_50_100_C; // member of StPicoHFJetMaker
      Hist_EEC_matched_50_100_C.reserve(3);

      std::vector<TH1D*> Hist_EEC_unmatched_5_10_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_5_10_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_unmatched_10_15_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_10_15_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_unmatched_15_20_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_15_20_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_unmatched_20_30_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_20_30_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_unmatched_30_50_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_30_50_C.reserve(3);
      std::vector<TH1D*> Hist_EEC_unmatched_50_100_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_50_100_C.reserve(3);
      
      std::vector<TH1D*> Hist_EEC_unmatched_all_C; // member of StPicoHFJetMaker
      Hist_EEC_unmatched_all_C.reserve(3);

      for (int c3 = 1; c3 <= 3; ++c3) { // loop over centrality classes (1..3); adjust as needed for your dataset
        TDirectory* cdir = configDir->mkdir(kCentTag[c3]); // create directory for this centrality class
        if (!cdir) cdir = (TDirectory*)configDir->Get(kCentTag[c3]); // if it already exists, get it
        cdir->cd();
      
        //ATH1D* hEEC = new TH1D("hEEC", "EEC vs RL;RL;EEC", N_bins, EEC_bounds); // create EEC histogram for this centrality class; adjust name/title as needed

        TTree* jetTree = new TTree("JetTree", "JetTree");
        //jetTree->SetDirectory(0); // detach from file - won't be auto-saved
        jetTree->Branch("runId", &fRunNumber, "runId/I");
        jetTree->Branch("eventId", &fEventId,   "eventId/I");
        jetTree->Branch("centralityWeight", &fCentralityWeight, "centralityWeight/F"); 
        if (mIsEmbedding) {
          jetTree->Branch("xsecWeight", &fXsecWeight, "xsecWeight/F");
          jetTree->Branch("deltaR", &fDeltaR, "deltaR/F");
          jetTree->Branch("mc_pt", &fMcJet.pt, "mc_pt/F");
          jetTree->Branch("mc_eta", &fMcJet.eta, "mc_eta/F");
          jetTree->Branch("mc_phi", &fMcJet.phi, "mc_phi/F");
          jetTree->Branch("mc_area", &fMcJet.area, "mc_area/F");
          jetTree->Branch("mc_pt_lead", &fMcJet.pt_lead, "mc_pt_lead/F");
          jetTree->Branch("mc_n_constituents", &fMcJet.n_constituents, "mc_n_constituents/I");
          jetTree->Branch("mc_neutral_fraction", &fMcJet.neutral_fraction, "mc_neutral_fraction/F");
          jetTree->Branch("mc_sum_pt", &fMcSumPt, "mc_sum_pt/F");
        
        }
        jetTree->Branch("reco_pt", &fRecoJet.pt, "reco_pt/F");
        jetTree->Branch("reco_pt_corr", &fRecoJet.pt_corr, "reco_pt_corr/F");
        jetTree->Branch("reco_eta", &fRecoJet.eta, "reco_eta/F");
        jetTree->Branch("reco_phi", &fRecoJet.phi, "reco_phi/F");
        jetTree->Branch("reco_area", &fRecoJet.area, "reco_area/F");
        jetTree->Branch("reco_rho", &fRecoJet.rho, "reco_rho/F");
        jetTree->Branch("reco_pt_lead", &fRecoJet.pt_lead, "reco_pt_lead/F");
        jetTree->Branch("reco_n_constituents", &fRecoJet.n_constituents, "reco_n_constituents/I");
        jetTree->Branch("reco_n_constituents_real", &fRecoJet.n_constituents_real, "reco_n_constituents_real/I");
        jetTree->Branch("reco_neutral_fraction", &fRecoJet.neutral_fraction, "reco_neutral_fraction/F");
        jetTree->Branch("reco_trigger_match", &fRecoJet.trigger_match, "reco_trigger_match/O");

        treesC.push_back(jetTree);
        if(!mIsEmbedding) {    
          TTree* constituentTree = new TTree("ConstituentTree", "Jet Constituents");
          //constituentTree->SetDirectory(0); // detach from file - won't be auto-saved
          constituentTree->Branch("runid", &c_runid, "runid/I");
          constituentTree->Branch("eventid", &c_eventid, "eventid/I");
          constituentTree->Branch("ijet", &c_ijet, "ijet/I");
          constituentTree->Branch("px", &c_px, "px/F");
          constituentTree->Branch("py", &c_py, "py/F");
          constituentTree->Branch("pz", &c_pz, "pz/F");
          constituentTree->Branch("pt", &c_pt, "pt/F");
          constituentTree->Branch("E", &c_E, "E/F");
          constituentTree->Branch("eta", &c_eta, "eta/F");
          constituentTree->Branch("phi", &c_phi, "phi/F");
          constituentTree->Branch("charge", &c_charge, "charge/I");
          
          ConstituentTreeC.push_back(constituentTree);


          //Tree for EEC
          TTree* EECTree = new TTree("EECTree", "EEC");
          EECTree->SetDirectory(0); // detach from file - won't be auto-saved
          EECTree->Branch("ijet", &eec_ijet, "ijet/I");
          EECTree->Branch("runid", &eec_runid, "runid/I");
          EECTree->Branch("eventid", &eec_eventid, "eventid/I");
          EECTree->Branch("eec", &eec_data, "eec/F");
          EECTree->Branch("RL", &RL, "RL/F");

          EECTreeC.push_back(EECTree);

          //Histograms for EEC
          //TH1D* hEEC = new TH1D("hEEC", "EEC vs RL;RL;EEC", N_bins, EEC_bounds);
          //Hist_EEC_C.push_back(hEEC);
          TH1D* hEEC_5_10 = new TH1D("hEEC_5_10", "EEC vs RL for 5<=pT<10;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_C_5_10.push_back(hEEC_5_10);
          TH1D* hEEC_10_15 = new TH1D("hEEC_10_15", "EEC vs RL for 10<=pT<15;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_C_10_15.push_back(hEEC_10_15);
          TH1D* hEEC_15_20 = new TH1D("hEEC_15_20", "EEC vs RL for 15<=pT<20;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_C_15_20.push_back(hEEC_15_20);
          TH1D* hEEC_20_30 = new TH1D("hEEC_20_30", "EEC vs RL for 20<=pT<30;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_C_20_30.push_back(hEEC_20_30);
          TH1D* hEEC_30_50 = new TH1D("hEEC_30_50", "EEC vs RL for 30<=pT<50;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_C_30_50.push_back(hEEC_30_50);
        }

        if(mIsEmbedding){

          TTree* MCJetTree = new TTree("MCJetTree", "MC Jet Tree");
          //MCJetTree->SetDirectory(0); // detach from file - won't be auto
          MCJetTree->Branch("runId", &fRunNumber, "runId/I");
          MCJetTree->Branch("eventId", &fEventId,   "eventId/I");
          MCJetTree->Branch("centralityWeight", &fCentralityWeight, "centralityWeight/F"); 
          MCJetTree->Branch("xsecWeight", &fXsecWeight, "xsecWeight/F");
          MCJetTree->Branch("mc_pt", &fMcJet.pt, "mc_pt/F");
          MCJetTree->Branch("mc_eta", &fMcJet.eta, "mc_eta/F");
          MCJetTree->Branch("mc_phi", &fMcJet.phi, "mc_phi/F");
          MCJetTree->Branch("mc_area", &fMcJet.area, "mc_area/F");
          MCJetTree->Branch("mc_pt_lead", &fMcJet.pt_lead, "mc_pt_lead/F");
          MCJetTree->Branch("mc_n_constituents", &fMcJet.n_constituents, "mc_n_constituents/I");
          MCJetTree->Branch("mc_n_constituents_real", &fMcJet.n_constituents_real, "mc_n_constituents_real/I");

          MCJetTreeC.push_back(MCJetTree);

          TTree* MatchedJet_constituentTree = new TTree("MatchedJet_ConstituentTree", "Constituents of Matched Jets");
          //MatchedJet_constituentTree->SetDirectory(0); // detach from file - won't be auto-saved
          MatchedJet_constituentTree->Branch("runid", &c_matched_runid, "runid/I");
          MatchedJet_constituentTree->Branch("eventid", &c_matched_eventid, "eventid/I");
          MatchedJet_constituentTree->Branch("ijet", &c_matched_ijet, "ijet/I");
          MatchedJet_constituentTree->Branch("px", &c_matched_px, "px/F");
          MatchedJet_constituentTree->Branch("py", &c_matched_py, "py/F");
          MatchedJet_constituentTree->Branch("pz", &c_matched_pz, "pz/F");
          MatchedJet_constituentTree->Branch("pt", &c_matched_pt, "pt/F");
          MatchedJet_constituentTree->Branch("E", &c_matched_E, "E/F");
          MatchedJet_constituentTree->Branch("eta", &c_matched_eta, "eta/F");
          MatchedJet_constituentTree->Branch("phi", &c_matched_phi, "phi/F");
          MatchedJet_constituentTree->Branch("charge", &c_matched_charge, "charge/I");

          MatchedJetConstituentTreeC.push_back(MatchedJet_constituentTree);

          TTree* UnmatchedJet_constituentTree = new TTree("UnmatchedJet_ConstituentTree", "Constituents of Unmatched Jets");
          //UnmatchedJet_constituentTree->SetDirectory(0); // detach from file - won't be auto-saved
          UnmatchedJet_constituentTree->Branch("runid", &c_unmatched_runid, "runid/I");
          UnmatchedJet_constituentTree->Branch("eventid", &c_unmatched_eventid, "eventid/I");
          UnmatchedJet_constituentTree->Branch("ijet", &c_unmatched_ijet, "ijet/I");
          UnmatchedJet_constituentTree->Branch("px", &c_unmatched_px, "px/F");
          UnmatchedJet_constituentTree->Branch("py", &c_unmatched_py, "py/F");
          UnmatchedJet_constituentTree->Branch("pz", &c_unmatched_pz, "pz/F");
          UnmatchedJet_constituentTree->Branch("pt", &c_unmatched_pt, "pt/F");
          UnmatchedJet_constituentTree->Branch("E", &c_unmatched_E, "E/F");
          UnmatchedJet_constituentTree->Branch("eta", &c_unmatched_eta, "eta/F");
          UnmatchedJet_constituentTree->Branch("phi", &c_unmatched_phi, "phi/F");
          UnmatchedJet_constituentTree->Branch("charge", &c_unmatched_charge, "charge/I");

          UnmatchedJetConstituentTreeC.push_back(UnmatchedJet_constituentTree);

          TTree* MCJet_constituentTree = new TTree("MCJet_ConstituentTree", "Constituents of MC Jets");
          //MCJet_constituentTree->SetDirectory(0); // detach from file - won't be auto-saved
          MCJet_constituentTree->Branch("runid", &c_mc_runid, "runid/I");
          MCJet_constituentTree->Branch("eventid", &c_mc_eventid, "eventid/I");
          MCJet_constituentTree->Branch("ijet", &c_mc_ijet, "ijet/I");
          MCJet_constituentTree->Branch("px", &c_mc_px, "px/F");
          MCJet_constituentTree->Branch("py", &c_mc_py, "py/F");
          MCJet_constituentTree->Branch("pz", &c_mc_pz, "pz/F");
          MCJet_constituentTree->Branch("pt", &c_mc_pt, "pt/F");
          MCJet_constituentTree->Branch("E", &c_mc_E, "E/F");
          MCJet_constituentTree->Branch("eta", &c_mc_eta, "eta/F");
          MCJet_constituentTree->Branch("phi", &c_mc_phi, "phi/F");
          MCJet_constituentTree->Branch("charge", &c_mc_charge, "charge/I");

          MCJetConstituentTreeC.push_back(MCJet_constituentTree);

          TTree* EECTree_unmatchedReco = new TTree("EECTree_UnmatchedRecoJets", "EEC for Unmatched Reco Jets");
          EECTree_unmatchedReco->SetDirectory(0); // detach from file - won't be auto-saved
          EECTree_unmatchedReco->Branch("ijet", &eec_unmatched_ijet, "ijet/I");
          EECTree_unmatchedReco->Branch("runid", &eec_unmatched_runid, "runid/I");
          EECTree_unmatchedReco->Branch("eventid", &eec_unmatched_eventid, "eventid/I");
          EECTree_unmatchedReco->Branch("eec", &eec_unmatched, "eec/F");
          EECTree_unmatchedReco->Branch("RL", &RL_unmatched, "RL/F");

          EECTreeunmatchedC.push_back(EECTree_unmatchedReco);

          TTree* EECTree_matchedjets = new TTree("EECTree_MatchedJets", "EEC for Matched Jets");
          EECTree_matchedjets->SetDirectory(0); // detach from file - won't be auto-saved
          EECTree_matchedjets->Branch("ijet", &eec_matched_ijet, "ijet/I");
          EECTree_matchedjets->Branch("runid", &eec_matched_runid, "runid/I");
          EECTree_matchedjets->Branch("eventid", &eec_matched_eventid, "eventid/I");
          EECTree_matchedjets->Branch("eec", &eec_matched, "eec/F");
          EECTree_matchedjets->Branch("RL", &RL_matched, "RL/F");

          EECTreematchedC.push_back(EECTree_matchedjets);


          TTree* EECTree_MC = new TTree("EECTree_MCJets", "EEC for MC Jets");
          EECTree_MC->SetDirectory(0); // detach from file - won't be auto-saved
          EECTree_MC->Branch("ijet", &eec_mc_ijet, "ijet/I");
          EECTree_MC->Branch("runid", &eec_mc_runid, "runid/I");
          EECTree_MC->Branch("eventid", &eec_mc_eventid, "eventid/I");
          EECTree_MC->Branch("eec", &eec_mc, "eec/F");
          EECTree_MC->Branch("RL", &RL_mc, "RL/F");

          EECTree_MC_C.push_back(EECTree_MC);
          

          TH1D* hEEC_MC_5_10 = new TH1D("hEEC_MC_5_10", "EEC vs RL for MC Jets with 5<=pT<10;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_5_10_C.push_back(hEEC_MC_5_10);
          TH1D* hEEC_MC_10_15 = new TH1D("hEEC_MC_10_15", "EEC vs RL for MC Jets with 10<=pT<15;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_10_15_C.push_back(hEEC_MC_10_15);
          TH1D* hEEC_MC_15_20 = new TH1D("hEEC_MC_15_20", "EEC vs RL for MC Jets with 15<=pT<20;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_15_20_C.push_back(hEEC_MC_15_20);
          TH1D* hEEC_MC_20_30 = new TH1D("hEEC_MC_20_30", "EEC vs RL for MC Jets with 20<=pT<30;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_20_30_C.push_back(hEEC_MC_20_30);
          TH1D* hEEC_MC_30_50 = new TH1D("hEEC_MC_30_50", "EEC vs RL for MC Jets with 30<=pT<50;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_30_50_C.push_back(hEEC_MC_30_50);
          TH1D* hEEC_MC_50_100 = new TH1D("hEEC_MC_50_100", "EEC vs RL for MC Jets with 50<=pT<100;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_MC_50_100_C.push_back(hEEC_MC_50_100);
          

          TH1D* hEEC_matched_5_10 = new TH1D("hEEC_matched_5_10", "EEC vs RL for Matched Reco Jets with 5<=pT<10;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_5_10_C.push_back(hEEC_matched_5_10);
          TH1D* hEEC_matched_10_15 = new TH1D("hEEC_matched_10_15", "EEC vs RL for Matched Reco Jets with 10<=pT<15;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_10_15_C.push_back(hEEC_matched_10_15);
          TH1D* hEEC_matched_15_20 = new TH1D("hEEC_matched_15_20", "EEC vs RL for Matched Reco Jets with 15<=pT<20;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_15_20_C.push_back(hEEC_matched_15_20);
          TH1D* hEEC_matched_20_30 = new TH1D("hEEC_matched_20_30", "EEC vs RL for Matched Reco Jets with 20<=pT<30;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_20_30_C.push_back(hEEC_matched_20_30);
          TH1D* hEEC_matched_30_50 = new TH1D("hEEC_matched_30_50", "EEC vs RL for Matched Reco Jets with 30<=pT<50;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_30_50_C.push_back(hEEC_matched_30_50);
          TH1D* hEEC_matched_50_100 = new TH1D("hEEC_matched_50_100", "EEC vs RL for Matched Reco Jets with 50<=pT<100;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_matched_50_100_C.push_back(hEEC_matched_50_100);

          TH1D* hEEC_unmatched_5_10 = new TH1D("hEEC_unmatched_5_10", "EEC vs RL for Unmatched Reco Jets with 5<=pT<10;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_5_10_C.push_back(hEEC_unmatched_5_10);
          TH1D* hEEC_unmatched_10_15 = new TH1D("hEEC_unmatched_10_15", "EEC vs RL for Unmatched Reco Jets with 10<=pT<15;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_10_15_C.push_back(hEEC_unmatched_10_15);
          TH1D* hEEC_unmatched_15_20 = new TH1D("hEEC_unmatched_15_20", "EEC vs RL for Unmatched Reco Jets with 15<=pT<20;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_15_20_C.push_back(hEEC_unmatched_15_20);
          TH1D* hEEC_unmatched_20_30 = new TH1D("hEEC_unmatched_20_30", "EEC vs RL for Unmatched Reco Jets with 20<=pT<30;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_20_30_C.push_back(hEEC_unmatched_20_30);
          TH1D* hEEC_unmatched_30_50 = new TH1D("hEEC_unmatched_30_50", "EEC vs RL for Unmatched Reco Jets with 30<=pT<50;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_30_50_C.push_back(hEEC_unmatched_30_50);
          TH1D* hEEC_unmatched_50_100 = new TH1D("hEEC_unmatched_50_100", "EEC vs RL for Unmatched Reco Jets with 50<=pT<100;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_50_100_C.push_back(hEEC_unmatched_50_100);

          TH1D* hEEC_unmatched_all = new TH1D("hEEC_unmatched_all", "EEC vs RL for All Unmatched Reco Jets;RL;EEC", nBinsEEC, EEC_bounds);
          Hist_EEC_unmatched_all_C.push_back(hEEC_unmatched_all);
        }

        configDir->cd();
      }
    
      treesConfig.push_back(treesC);
      if (!mIsEmbedding) { //only create EEC trees and histograms for data, not embedding
        ConstituentTreeConfig.push_back(ConstituentTreeC);
        EECTreeConfig.push_back(EECTreeC);
        HistEEC_5_10_Config.push_back(Hist_EEC_C_5_10);
        HistEEC_10_15_Config.push_back(Hist_EEC_C_10_15);
        HistEEC_15_20_Config.push_back(Hist_EEC_C_15_20);
        HistEEC_20_30_Config.push_back(Hist_EEC_C_20_30);
        HistEEC_30_50_Config.push_back(Hist_EEC_C_30_50);
      }

      if (mIsEmbedding) {
        MCJetTreeConfig.push_back(MCJetTreeC);
        EECTreematchedConfig.push_back(EECTreematchedC);
        EECTreeunmatchedConfig.push_back(EECTreeunmatchedC);
        EECTree_MC_Config.push_back(EECTree_MC_C);
        MatchedJetConstituentTreeConfig.push_back(MatchedJetConstituentTreeC);
        UnmatchedJetConstituentTreeConfig.push_back(UnmatchedJetConstituentTreeC);
        MCJetConstituentTreeConfig.push_back(MCJetConstituentTreeC);

        HistEEC_MC_5_10_Config.push_back(Hist_EEC_MC_5_10_C);
        HistEEC_MC_10_15_Config.push_back(Hist_EEC_MC_10_15_C);
        HistEEC_MC_15_20_Config.push_back(Hist_EEC_MC_15_20_C);
        HistEEC_MC_20_30_Config.push_back(Hist_EEC_MC_20_30_C);
        HistEEC_MC_30_50_Config.push_back(Hist_EEC_MC_30_50_C);
        HistEEC_MC_50_100_Config.push_back(Hist_EEC_MC_50_100_C);

        HistEEC_matched_5_10_Config.push_back(Hist_EEC_matched_5_10_C);
        HistEEC_matched_10_15_Config.push_back(Hist_EEC_matched_10_15_C);
        HistEEC_matched_15_20_Config.push_back(Hist_EEC_matched_15_20_C);
        HistEEC_matched_20_30_Config.push_back(Hist_EEC_matched_20_30_C);
        HistEEC_matched_30_50_Config.push_back(Hist_EEC_matched_30_50_C);
        HistEEC_matched_50_100_Config.push_back(Hist_EEC_matched_50_100_C);

        HistEEC_unmatched_5_10_Config.push_back(Hist_EEC_unmatched_5_10_C);
        HistEEC_unmatched_10_15_Config.push_back(Hist_EEC_unmatched_10_15_C);
        HistEEC_unmatched_15_20_Config.push_back(Hist_EEC_unmatched_15_20_C);
        HistEEC_unmatched_20_30_Config.push_back(Hist_EEC_unmatched_20_30_C);
        HistEEC_unmatched_30_50_Config.push_back(Hist_EEC_unmatched_30_50_C);
        HistEEC_unmatched_50_100_Config.push_back(Hist_EEC_unmatched_50_100_C);

        HistEEC_unmatched_all_Config.push_back(Hist_EEC_unmatched_all_C);
      }

      rdir->cd();

    }


    fTreeRC.push_back(treesConfig);
    if(!mIsEmbedding) { //only create EEC trees and histograms for data, not embedding
      fConstituentTreeRC.push_back(ConstituentTreeConfig);
      fEECTreeRC.push_back(EECTreeConfig);
      //fHistEEC.push_back(Hist_EEC_C);
      fHistEEC_5_10.push_back(HistEEC_5_10_Config);
      fHistEEC_10_15.push_back(HistEEC_10_15_Config);
      fHistEEC_15_20.push_back(HistEEC_15_20_Config);
      fHistEEC_20_30.push_back(HistEEC_20_30_Config);
      fHistEEC_30_50.push_back(HistEEC_30_50_Config);
    }

    if(mIsEmbedding) {
      MCJetTreeRC.push_back(MCJetTreeConfig);

      fEECTreematchedRC.push_back(EECTreematchedConfig);
      fEECTreeunmatchedRC.push_back(EECTreeunmatchedConfig);
      fEECTree_MC_RC.push_back(EECTree_MC_Config);
      fMatchedJetConstituentTreeRC.push_back(MatchedJetConstituentTreeConfig);
      fUnmatchedJetConstituentTreeRC.push_back(UnmatchedJetConstituentTreeConfig);
      fMCJetConstituentTreeRC.push_back(MCJetConstituentTreeConfig);

      fHistEEC_MC_5_10.push_back(HistEEC_MC_5_10_Config);
      fHistEEC_MC_10_15.push_back(HistEEC_MC_10_15_Config);
      fHistEEC_MC_15_20.push_back(HistEEC_MC_15_20_Config);
      fHistEEC_MC_20_30.push_back(HistEEC_MC_20_30_Config);
      fHistEEC_MC_30_50.push_back(HistEEC_MC_30_50_Config);
      fHistEEC_MC_50_100.push_back(HistEEC_MC_50_100_Config);

      fHistEEC_matched_5_10.push_back(HistEEC_matched_5_10_Config);
      fHistEEC_matched_10_15.push_back(HistEEC_matched_10_15_Config);
      fHistEEC_matched_15_20.push_back(HistEEC_matched_15_20_Config);
      fHistEEC_matched_20_30.push_back(HistEEC_matched_20_30_Config);
      fHistEEC_matched_30_50.push_back(HistEEC_matched_30_50_Config);
      fHistEEC_matched_50_100.push_back(HistEEC_matched_50_100_Config);

      fHistEEC_unmatched_5_10.push_back(HistEEC_unmatched_5_10_Config);
      fHistEEC_unmatched_10_15.push_back(HistEEC_unmatched_10_15_Config);
      fHistEEC_unmatched_15_20.push_back(HistEEC_unmatched_15_20_Config);
      fHistEEC_unmatched_20_30.push_back(HistEEC_unmatched_20_30_Config);
      fHistEEC_unmatched_30_50.push_back(HistEEC_unmatched_30_50_Config);
      fHistEEC_unmatched_50_100.push_back(HistEEC_unmatched_50_100_Config);

      fHistEEC_unmatched_all.push_back(HistEEC_unmatched_all_Config);
    }

    fileDir->cd();
  }

  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt = "") { return; }

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  TDirectory* fileDir = gDirectory;
  if (!fileDir) return kStOK;

  const size_t nR = fR.size();

for (size_t iR = 0; iR < nR; ++iR) {
  TDirectory* rdir = dynamic_cast<TDirectory*>(fileDir->Get(Form("R%.1f", fR[iR])));
  if (!rdir) continue;

  for (size_t iConfig = 0; iConfig < fConfigNames.size(); ++iConfig) {
    TDirectory* configDir = dynamic_cast<TDirectory*>(rdir->Get(fConfigNames[iConfig]));
    if (!configDir) continue;

    for (int c3 = 1; c3 <= 3; ++c3) {
      const int ciTree = c3 - 1;  // 0..2

      TDirectory* cdir = dynamic_cast<TDirectory*>(configDir->Get(kCentTag[c3]));
      if (!cdir) continue;
      cdir->cd();

      // --- write the tree
      if (iR < fTreeRC.size() && ciTree >= 0 &&
        ciTree < (int)fTreeRC[iR][iConfig].size() &&
        fTreeRC[iR][iConfig][ciTree]) {
          fTreeRC[iR][iConfig][ciTree]->Write();
      }

      if (!mIsEmbedding) { //only write EEC trees and histograms for data, not embedding
        /*
        if (iR < fConstituentTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)fConstituentTreeRC[iR][iConfig].size() &&
          fConstituentTreeRC[iR][iConfig][ciTree]) {
          fConstituentTreeRC[iR][iConfig][ciTree]->Write();
        }

        if (iR < fEECTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)fEECTreeRC[iR][iConfig].size() &&
          fEECTreeRC[iR][iConfig][ciTree]) {
          fEECTreeRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }
        */
        // --- Normalize histograms for each R and centrality
        if(iR < fHistEEC_5_10.size() && ciTree >= 0 && ciTree < (int)fHistEEC_5_10[iR][iConfig].size()) {
          TH1D* h = fHistEEC_5_10[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_10_15.size() && ciTree >= 0 && ciTree < (int)fHistEEC_10_15[iR][iConfig].size()) {
          TH1D* h = fHistEEC_10_15[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_15_20.size() && ciTree >= 0 && ciTree < (int)fHistEEC_15_20[iR][iConfig].size()) {
          TH1D* h = fHistEEC_15_20[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_20_30.size() && ciTree >= 0 && ciTree < (int)fHistEEC_20_30[iR][iConfig].size()) {
          TH1D* h = fHistEEC_20_30[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_30_50.size() && ciTree >= 0 && ciTree < (int)fHistEEC_30_50[iR][iConfig].size()) {
          TH1D* h = fHistEEC_30_50[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        /*
        if (iR < fHistEEC.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC[iR].size() &&
          fHistEEC[iR][iConfig][ciTree]) {
          fHistEEC[iR][iConfig][ciTree]->Write();
        }
        */
    
        if (iR < fHistEEC_5_10.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_5_10[iR][iConfig].size() &&
          fHistEEC_5_10[iR][iConfig][ciTree]) {
          fHistEEC_5_10[iR][iConfig][ciTree]->Write();
        }
    
        if (iR < fHistEEC_10_15.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_10_15[iR][iConfig].size() &&
          fHistEEC_10_15[iR][iConfig][ciTree]) {
          fHistEEC_10_15[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_15_20.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_15_20[iR][iConfig].size() &&
          fHistEEC_15_20[iR][iConfig][ciTree]) {
          fHistEEC_15_20[iR][iConfig][ciTree]->Write();
        }
    
        if (iR < fHistEEC_20_30.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_20_30[iR][iConfig].size() &&
          fHistEEC_20_30[iR][iConfig][ciTree]) {
          fHistEEC_20_30[iR][iConfig][ciTree]->Write();
        }
    
        if (iR < fHistEEC_30_50.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_30_50[iR][iConfig].size() &&
          fHistEEC_30_50[iR][iConfig][ciTree]) {
          fHistEEC_30_50[iR][iConfig][ciTree]->Write();
        }
      }
      
      //Embedding trees and histograms are written for both data and embedding, but only if mIsEmbedding is true. This allows us to have the same output structure for both cases, while only filling the embedding-specific trees and histograms when we are actually running on embedding.
      if(mIsEmbedding){
        if(iR < MCJetTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)MCJetTreeRC[iR][iConfig].size() &&
          MCJetTreeRC[iR][iConfig][ciTree]) {
          MCJetTreeRC[iR][iConfig][ciTree]->Write();
        }

        if (iR < fEECTreematchedRC.size() && ciTree >= 0 &&
          ciTree < (int)fEECTreematchedRC[iR][iConfig].size() &&
          fEECTreematchedRC[iR][iConfig][ciTree]) {
          fEECTreematchedRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        if (iR < fEECTreeunmatchedRC.size() && ciTree >= 0 &&
          ciTree < (int)fEECTreeunmatchedRC[iR][iConfig].size() &&
          fEECTreeunmatchedRC[iR][iConfig][ciTree]) {
          fEECTreeunmatchedRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        if (iR < fEECTree_MC_RC.size() && ciTree >= 0 &&
          ciTree < (int)fEECTree_MC_RC[iR][iConfig].size() &&
          fEECTree_MC_RC[iR][iConfig][ciTree]) {
          fEECTree_MC_RC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        if (iR < fMatchedJetConstituentTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)fMatchedJetConstituentTreeRC[iR][iConfig].size() &&
          fMatchedJetConstituentTreeRC[iR][iConfig][ciTree]) {
          fMatchedJetConstituentTreeRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        if (iR < fUnmatchedJetConstituentTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)fUnmatchedJetConstituentTreeRC[iR][iConfig].size() &&
          fUnmatchedJetConstituentTreeRC[iR][iConfig][ciTree]) {
          fUnmatchedJetConstituentTreeRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        if (iR < fMCJetConstituentTreeRC.size() && ciTree >= 0 &&
          ciTree < (int)fMCJetConstituentTreeRC[iR][iConfig].size() &&
          fMCJetConstituentTreeRC[iR][iConfig][ciTree]) {
          fMCJetConstituentTreeRC[iR][iConfig][ciTree]->Write("", TObject::kOverwrite);
        }

        //Histogram normalization for embedding histograms
        if (iR < fHistEEC_MC_5_10.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_5_10[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_5_10[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_MC_10_15.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_10_15[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_10_15[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_MC_15_20.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_15_20[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_15_20[iR][iConfig][ciTree];
          double integral = h->Integral("width"); // includes bin width automatically
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_MC_20_30.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_20_30[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_20_30[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_MC_30_50.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_30_50[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_30_50[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_MC_50_100.size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_50_100[iR][iConfig].size()) {
          TH1D* h = fHistEEC_MC_50_100[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_5_10.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_5_10[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_5_10[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_10_15.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_10_15[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_10_15[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_15_20.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_15_20[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_15_20[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_20_30.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_20_30[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_20_30[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_30_50.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_30_50[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_30_50[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_matched_50_100.size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_50_100[iR][iConfig].size()) {
          TH1D* h = fHistEEC_matched_50_100[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_5_10.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_5_10[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_5_10[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_10_15.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_10_15[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_10_15[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_15_20.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_15_20[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_15_20[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_20_30.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_20_30[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_20_30[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_30_50.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_30_50[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_30_50[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_50_100.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_50_100[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_50_100[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }

        if (iR < fHistEEC_unmatched_all.size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_all[iR][iConfig].size()) {
          TH1D* h = fHistEEC_unmatched_all[iR][iConfig][ciTree];
          double integral = h->Integral("width");
          if (integral > 0) h->Scale(1.0 / integral);
        }


        //Histogram Writing for embedding histograms
        
        if (iR < fHistEEC_matched_5_10.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_5_10[iR][iConfig].size() &&
          fHistEEC_matched_5_10[iR][iConfig][ciTree]) {
          fHistEEC_matched_5_10[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_matched_10_15.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_10_15[iR][iConfig].size() &&
          fHistEEC_matched_10_15[iR][iConfig][ciTree]) {
          fHistEEC_matched_10_15[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_matched_15_20.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_15_20[iR][iConfig].size() &&
          fHistEEC_matched_15_20[iR][iConfig][ciTree]) {
          fHistEEC_matched_15_20[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_matched_20_30.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_20_30[iR][iConfig].size() &&
          fHistEEC_matched_20_30[iR][iConfig][ciTree]) {
          fHistEEC_matched_20_30[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_matched_30_50.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_30_50[iR][iConfig].size() &&
          fHistEEC_matched_30_50[iR][iConfig][ciTree]) {
          fHistEEC_matched_30_50[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_matched_50_100.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_matched_50_100[iR][iConfig].size() &&
          fHistEEC_matched_50_100[iR][iConfig][ciTree]) {
          fHistEEC_matched_50_100[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_5_10.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_5_10[iR][iConfig].size() &&
          fHistEEC_MC_5_10[iR][iConfig][ciTree]) {
          fHistEEC_MC_5_10[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_10_15.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_10_15[iR][iConfig].size() &&
          fHistEEC_MC_10_15[iR][iConfig][ciTree]) {
          fHistEEC_MC_10_15[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_15_20.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_15_20[iR][iConfig].size() &&
          fHistEEC_MC_15_20[iR][iConfig][ciTree]) {
          fHistEEC_MC_15_20[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_20_30.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_20_30[iR][iConfig].size() &&
          fHistEEC_MC_20_30[iR][iConfig][ciTree]) {
          fHistEEC_MC_20_30[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_30_50.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_30_50[iR][iConfig].size() &&
          fHistEEC_MC_30_50[iR][iConfig][ciTree]) {
          fHistEEC_MC_30_50[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_MC_50_100.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_MC_50_100[iR][iConfig].size() &&
          fHistEEC_MC_50_100[iR][iConfig][ciTree]) {
          fHistEEC_MC_50_100[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_5_10.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_5_10[iR][iConfig].size() &&
          fHistEEC_unmatched_5_10[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_5_10[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_10_15.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_10_15[iR][iConfig].size() &&
          fHistEEC_unmatched_10_15[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_10_15[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_15_20.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_15_20[iR][iConfig].size() &&
          fHistEEC_unmatched_15_20[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_15_20[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_20_30.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_20_30[iR][iConfig].size() &&
          fHistEEC_unmatched_20_30[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_20_30[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_30_50.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_30_50[iR][iConfig].size() &&
          fHistEEC_unmatched_30_50[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_30_50[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_50_100.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_50_100[iR][iConfig].size() &&
          fHistEEC_unmatched_50_100[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_50_100[iR][iConfig][ciTree]->Write();
        }

        if (iR < fHistEEC_unmatched_all.size() && ciTree >= 0 &&
          ciTree < (int)fHistEEC_unmatched_all[iR][iConfig].size() &&
          fHistEEC_unmatched_all[iR][iConfig][ciTree]) {
          fHistEEC_unmatched_all[iR][iConfig][ciTree]->Write();
        }
      }
    } // c3
 
  } // charged/full

}   // iR
  fileDir->cd();
  return kStOK;
}

// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {

    for (int i = 0; i < 4800; ++i) {
      Sump[i] = 0.0;
    }

  fMcSumPt = 0.0f; 
  Bool_t vetoReco = kFALSE; 
  
  TH1D *hcent9 = static_cast<TH1D *>(mOutList->FindObject("hcent9"));

  vector<fastjet::PseudoJet> jetTracks;
  vector<fastjet::PseudoJet> neutraljetTracks; // from bemc towers only
  vector<fastjet::PseudoJet> fullTracks;
  vector<fastjet::PseudoJet> MCjetTracks;

  fRunNumber = mPicoDst->event()->runId(); // runID
  fEventId = mPicoDst->event()->eventId(); // eventID
  int refMult = mPicoDst->event()->refMult(); // refMult for centrality determination
  double vz = mPrimVtx.z(); // vertex z for centrality determination
  mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(),
                             vz);
  fCentrality = mRefmultCorrUtil->centrality9(); // 0 = 0-5 %,..., 8 = 70-80
                                                 // %
  if (fCentrality == -1)
    return kStOK; // no fCentrality

  fCentralityWeight = mRefmultCorrUtil->weight();

  if (hcent9) hcent9->Fill(fCentrality, fCentralityWeight);

  TH1D* hTowE = (TH1D*)mOutList->FindObject("hTowE");
  TH2D* hTowEtaPhi = (TH2D*)mOutList->FindObject("hTowEtaPhi");
  TH1D* hTrackPt = (TH1D*)mOutList->FindObject("hTrackPt");
  TH2D* hTrackEtaPhi = (TH2D*)mOutList->FindObject("hTrackEtaPhi");

  // map to 3 classes: 1=0-10%, 2=20-40%, 3=60-80%, else 0 (=skip)
  int c3 = 0;
  if (fCentrality == 0 || fCentrality == 1)       c3 = 1; // 0-10%
  else if (fCentrality == 3 || fCentrality == 4)  c3 = 2; // 20-40%
  else if (fCentrality == 7 || fCentrality == 8)  c3 = 3; // 60-80%
  if (c3 == 0) {
  return kStOK;
  }

  // MC tracks
  int noMCtracks = mPicoDst->numberOfMcTracks();
  for (int i = 0; i < noMCtracks; i++) {
    StPicoMcTrack *mctrk = (StPicoMcTrack *)mPicoDst->mcTrack(i);
    if (mctrk->idVtxStart() > 1)
      continue; // only primary tracks
    int geantId = mctrk->geantId();
    double mcpt = mctrk->pt();
    double mceta = mctrk->eta();
    if ((geantId > 3 && geantId < 7) || fabs(mceta) > 1.0 || mcpt < 0.2)
      continue;
    fMcSumPt += (Float_t)mcpt;
    TVector3 mcmom = mctrk->p();
    double mcphi = mcmom.Phi();
    if (mcphi < 0.0)
      mcphi += 2.0 * TMath::Pi();
    if (mcphi > 2.0 * TMath::Pi())
      mcphi -= 2.0 * TMath::Pi();

    double mcpx, mcpy, mcpz;
    mcpx = mcmom.x();
    mcpy = mcmom.y();
    mcpz = mcmom.z();

    double mcE = mctrk->energy();

    fastjet::PseudoJet inputMcParticle(mcpx, mcpy, mcpz, mcE);
    if (mctrk->charge() == 0) {
      inputMcParticle.set_user_index(0);
    } else
      inputMcParticle.set_user_index(i);
    MCjetTracks.push_back(inputMcParticle);
  }

if (mIsEmbedding && fpThatmax > 0.0 && !MCjetTracks.empty()) {
  float Rcheck = fR.empty() ? 0.4f : *std::max_element(fR.begin(), fR.end());
  fastjet::JetDefinition mc_jet_def_veto(fastjet::antikt_algorithm, Rcheck);
  fastjet::ClusterSequence mc_cs_veto(MCjetTracks, mc_jet_def_veto);
  std::vector<fastjet::PseudoJet> mcjets_veto =
      sorted_by_pt(mc_cs_veto.inclusive_jets(1.0)); // pT > 1 GeV

  const double ptMaxVeto = 1.5 * fpThatmax;

  for (size_t i = 0; i < mcjets_veto.size(); ++i) {
    if (mcjets_veto[i].perp() > ptMaxVeto) {
      return kStOK;
    }
  }

}
  // RC part
  GetCaloTrackMomentum(mPicoDst, mPrimVtx); // fill array Sump with momenta of tracks which are matched to BEMC

  StEmcGeom *mEmcGeom = StEmcGeom::getEmcGeom("bemc");
  StEmcPosition *mEmcPosition = new StEmcPosition();

//  double TOWE = 0;
for (int iTow = 0; iTow < 4800; iTow++) { // get btow info
  StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
  if (!towHit || towHit->isBad())
    continue;
  int realtowID = towHit->numericIndex2SoftId(iTow);
  if (BadTowerMap[realtowID])
    continue;

  double towE = GetTowerCalibEnergy(iTow + 1);

  if (doTowErrPlus == true)  towE = towE + 0.038 * towE;
  if (doTowErrMinus == true) towE = towE - 0.038 * towE;

  towE -= fHadronCorr * Sump[iTow];
  if (towE < 0) towE = 0;

  float Toweta_tmp = 0, Towphi = 0;
  if (mEmcGeom) mEmcGeom->getEtaPhi(realtowID, Toweta_tmp, Towphi);

  StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(
      StThreeVectorF(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z()), realtowID);

  if (Towphi < 0)              Towphi += 2.0 * TMath::Pi();
  if (Towphi >= 2.0*TMath::Pi()) Towphi -= 2.0 * TMath::Pi();

  float Toweta = towerPosition.pseudoRapidity();
  double ET = towE / cosh(Toweta);

  if (hTowE) hTowE->Fill(towE, fCentralityWeight);
  if (hTowEtaPhi) hTowEtaPhi->Fill(Toweta, Towphi, fCentralityWeight);


  if (ET > 30.0) {
    vetoReco = kTRUE;
    break;   // no need to check further towers
  }

  // no clustering if vetoReco; we will check this later
  double px = ET * cos(Towphi);
  double py = ET * sin(Towphi);
  double pz = towE * tanh(Toweta);

  fastjet::PseudoJet inputTower(px, py, pz, towE);
  if (inputTower.perp() > fETmincut) {
    inputTower.set_user_index(0); // neutral
    int ADC = towHit->adc() >> 4;
    //if (ADC > fTrgthresh) {
      //inputTower.set_user_index(9999); // trigger towers
    //}
    neutraljetTracks.push_back(inputTower);
  }
} // end tower loop

delete mEmcPosition;


// loop over primary tracks
for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
  StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);

  if (doTrackErr) {
    static TRandom3 randGen;
    if (randGen.Rndm() > 0.96) continue;
  }

  const TVector3 p = trk->pMom();
  const double pT  = p.Perp();
  if (!(pT > 0)) continue;

  if (pT > 30.0) {
    vetoReco = kTRUE;
    break;   // no need to look at other tracks
  }

  const float eta = p.PseudoRapidity();
  if (fabs(eta) > 1.0) continue;
  float phi = trk->pMom().Phi();
  if (phi < 0) phi += 2.0*TMath::Pi();/*  */
  if (phi >= 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();

  if (hTrackPt) hTrackPt->Fill(pT, fCentralityWeight);
  if (hTrackEtaPhi) hTrackEtaPhi->Fill(eta, phi, fCentralityWeight);

  float dca = (mPrimVtx - trk->origin()).Mag();

  float charged = trk->charge();
  (void)dca; // to avoid unused variable warning if not used in cuts

  fastjet::PseudoJet pj(p.x(), p.y(), p.z(), p.Mag());
  pj.set_user_index(charged); // for reco tracks, user_index = charge (will be 0 for neutrals from towers)

  /*
  if (mIsEmbedding) {
    if (trk->qaTruth() > 95) pj.set_user_index(trk->idTruth() - 1);
    else                     pj.set_user_index(trk->charge() ? 1 : 0);
  } else {
    pj.set_user_index(trk->charge() ? 1 : 0);
  } 
  */
  
  jetTracks.push_back(pj);

} // end loop over primary tracks

fullTracks.clear();

if (!vetoReco) {
  // build fullTracks from towers + tracks
  fullTracks = neutraljetTracks;
  fullTracks.insert(fullTracks.end(), jetTracks.begin(), jetTracks.end());
}

//==================================================================================//
// Jet part
//==================================================================================//
fastjet::AreaDefinition area_def(
    fastjet::active_area_explicit_ghosts,
    fastjet::GhostedAreaSpec(fGhostMaxrap, 1, 0.001));
//====================background estimate=======================//
float rho_full = 0.0;
float rho_charged = 0.0;

TH2D* hRhoVsRefMultfull = (TH2D*)mOutList->FindObject("hRhoVsRefMultfull");
TH2D* hRhoVsRefMultcharged = (TH2D*)mOutList->FindObject("hRhoVsRefMultcharged");

fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fRBg);
nJetsRemove = (c3 == 1 ? 2 : 1);

fastjet::Selector selector = (!fastjet::SelectorNHardest(nJetsRemove)) *
                             fastjet::SelectorAbsEtaMax(1.0) *
                             fastjet::SelectorPtMin(0.01);

if (!vetoReco && !fullTracks.empty()) { // comment for charged jets, as we want to use only tracks for rho estimation
  fastjet::JetMedianBackgroundEstimator bkgd_estimator(
      selector, jet_def_for_rho, area_def);
  bkgd_estimator.set_particles(fullTracks); // comment for charged jets, as we want to use only tracks for rho estimation
  rho_full = bkgd_estimator.rho();
}

//if (!vetoReco && !fullTracks.empty()) { // comment for charged jets, as we want to use only tracks for rho estimation
if (!vetoReco && !jetTracks.empty()) { //comment for full jets, use jetTracks for rho estimation to be consistent with charged jet case
  fastjet::JetMedianBackgroundEstimator bkgd_estimator(
      selector, jet_def_for_rho, area_def);
  //bkgd_estimator.set_particles(fullTracks); // comment for charged jets, as we want to use only tracks for rho estimation
  bkgd_estimator.set_particles(jetTracks); //comment for full jets, use jetTracks for rho estimation to be consistent with charged jet case
  rho_charged = bkgd_estimator.rho();
}

struct JetConfig {
  TString name;
  std::vector<fastjet::PseudoJet>* tracks;
  float rho;
};

std::vector<JetConfig> jetConfigs = {
  {fConfigNames[0], &fullTracks, rho_full},
  {fConfigNames[1], &jetTracks, rho_charged}
};

if (hRhoVsRefMultfull && !vetoReco && !fullTracks.empty()) { //comment for charged jets, as we want to use only tracks for rho estimation
  hRhoVsRefMultfull->Fill(refMult, rho_full, fCentralityWeight);
}

if (hRhoVsRefMultcharged && !vetoReco && !jetTracks.empty()) {  //comment for full jets, use jetTracks for rho estimation to be consistent with charged jet case
  hRhoVsRefMultcharged->Fill(refMult, rho_charged, fCentralityWeight);
}

//======================================================================//
const double ptMaxVeto = (mIsEmbedding && fpThatmax > 0.0) ? (1.5 * fpThatmax) : -1.0;

for (size_t iConfig = 0; iConfig < jetConfigs.size(); ++iConfig) {
  std::vector<fastjet::PseudoJet>& tracks = *jetConfigs[iConfig].tracks;
  float rho = jetConfigs[iConfig].rho;
  
  for (unsigned int i = 0; i < fR.size(); i++) {
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, fR[i]);
    float maxRapJet = 1 - fR[i];

    //EEC_low = 0.0;
    //EEC_mid = 0.0;
    //EEC_high = 0.0;

    //==============================Reco jets===============================//

    //Create Trees for each R, config, centrality class
    TTree* jetTree = 0;
    const int ciTree = c3 - 1;
    if (i < fTreeRC.size() && iConfig < fTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fTreeRC[i][iConfig].size())
      jetTree = fTreeRC[i][iConfig][ciTree];

    TTree* constituentTree = 0;
    if (i < fConstituentTreeRC.size() && iConfig < fConstituentTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fConstituentTreeRC[i][iConfig].size())
      constituentTree = fConstituentTreeRC[i][iConfig][ciTree];
  
    TTree* try_EECTree = 0;
    if (i < fEECTreeRC.size() && iConfig < fEECTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fEECTreeRC[i][iConfig].size())
    try_EECTree = fEECTreeRC[i][iConfig][ciTree];
    
    //Create EEC histograms for each R, config, centrality class 

    //TH1D* hEEC = 0;
    //if (i < fHistEEC.size() && ciTree >= 0 && ciTree < (int)fHistEEC[i].size())
    //hEEC = fHistEEC[i][ciTree];

    TH1D* hEEC_5_10 = 0;
    if (i < fHistEEC_5_10.size() && iConfig < fHistEEC_5_10[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_5_10[i][iConfig].size())
      hEEC_5_10 = fHistEEC_5_10[i][iConfig][ciTree];

    TH1D* hEEC_10_15 = 0;
    if (i < fHistEEC_10_15.size() && iConfig < fHistEEC_10_15[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_10_15[i][iConfig].size())
      hEEC_10_15 = fHistEEC_10_15[i][iConfig][ciTree];

    TH1D* hEEC_15_20 = 0;
    if (i < fHistEEC_15_20.size() && iConfig < fHistEEC_15_20[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_15_20[i][iConfig].size())
      hEEC_15_20 = fHistEEC_15_20[i][iConfig][ciTree];
  
    TH1D* hEEC_20_30 = 0;
    if (i < fHistEEC_20_30.size() && iConfig < fHistEEC_20_30[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_20_30[i][iConfig].size())
      hEEC_20_30 = fHistEEC_20_30[i][iConfig][ciTree];

    TH1D* hEEC_30_50 = 0;
    if (i < fHistEEC_30_50.size() && iConfig < fHistEEC_30_50[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_30_50[i][iConfig].size())
      hEEC_30_50 = fHistEEC_30_50[i][iConfig][ciTree];
    
    TTree* MCJetTree = 0;
    if (i < MCJetTreeRC.size() && iConfig < MCJetTreeRC[i].size() && ciTree >= 0 && ciTree < (int)MCJetTreeRC[i][iConfig].size())
      MCJetTree = MCJetTreeRC[i][iConfig][ciTree];


    //Create Trees for embedding for R, config, centrality class
    TTree* MatchedConstituentTree = 0;
    if (i < fMatchedJetConstituentTreeRC.size() && iConfig < fMatchedJetConstituentTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fMatchedJetConstituentTreeRC[i][iConfig].size())
      MatchedConstituentTree = fMatchedJetConstituentTreeRC[i][iConfig][ciTree];
    
    TTree* UnmatchedConstituentTree = 0;
    if (i < fUnmatchedJetConstituentTreeRC.size() && iConfig < fUnmatchedJetConstituentTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fUnmatchedJetConstituentTreeRC[i][iConfig].size())
      UnmatchedConstituentTree = fUnmatchedJetConstituentTreeRC[i][iConfig][ciTree];

    TTree* MCConstituentTree = 0;
    if (i < fMCJetConstituentTreeRC.size() && iConfig < fMCJetConstituentTreeRC[i].size() && ciTree >= 0 && ciTree < (int)fMCJetConstituentTreeRC[i][iConfig].size())
      MCConstituentTree = fMCJetConstituentTreeRC[i][iConfig][ciTree];
    
    TTree* EECTree_matched = 0;
    if (i < fEECTreematchedRC.size() && iConfig < fEECTreematchedRC[i].size() && ciTree >= 0 && ciTree < (int)fEECTreematchedRC[i][iConfig].size())
      EECTree_matched = fEECTreematchedRC[i][iConfig][ciTree];

    TTree* EECTree_MC = 0;
    if (i < fEECTree_MC_RC.size() && iConfig < fEECTree_MC_RC[i].size() && ciTree >= 0 && ciTree < (int)fEECTree_MC_RC[i][iConfig].size())
      EECTree_MC = fEECTree_MC_RC[i][iConfig][ciTree];

    TTree* EECTree_unmatched = 0;
    if (i < fEECTreeunmatchedRC.size() && iConfig < fEECTreeunmatchedRC[i].size() && ciTree >= 0 && ciTree < (int)fEECTreeunmatchedRC[i][iConfig].size())
      EECTree_unmatched = fEECTreeunmatchedRC[i][iConfig][ciTree];

    
    TH1D* hEEC_matched_5_10 = 0;
    if (i < fHistEEC_matched_5_10.size() && iConfig < fHistEEC_matched_5_10[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_5_10[i][iConfig].size())
      hEEC_matched_5_10 = fHistEEC_matched_5_10[i][iConfig][ciTree];

    TH1D* hEEC_matched_10_15 = 0;
    if (i < fHistEEC_matched_10_15.size() && iConfig < fHistEEC_matched_10_15[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_10_15[i][iConfig].size())
      hEEC_matched_10_15 = fHistEEC_matched_10_15[i][iConfig][ciTree];

    TH1D* hEEC_matched_15_20 = 0;
    if (i < fHistEEC_matched_15_20.size() && iConfig < fHistEEC_matched_15_20[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_15_20[i][iConfig].size())
      hEEC_matched_15_20 = fHistEEC_matched_15_20[i][iConfig][ciTree];

    TH1D* hEEC_matched_20_30 = 0;
    if (i < fHistEEC_matched_20_30.size() && iConfig < fHistEEC_matched_20_30[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_20_30[i][iConfig].size())
      hEEC_matched_20_30 = fHistEEC_matched_20_30[i][iConfig][ciTree];

    TH1D* hEEC_matched_30_50 = 0;
    if (i < fHistEEC_matched_30_50.size() && iConfig < fHistEEC_matched_30_50[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_30_50[i][iConfig].size())
      hEEC_matched_30_50 = fHistEEC_matched_30_50[i][iConfig][ciTree];

    TH1D* hEEC_matched_50_100 = 0;
    if (i < fHistEEC_matched_50_100.size() && iConfig < fHistEEC_matched_50_100[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_matched_50_100[i][iConfig].size())
      hEEC_matched_50_100 = fHistEEC_matched_50_100[i][iConfig][ciTree];

    TH1D* hEEC_MC_5_10 = 0;
    if (i < fHistEEC_MC_5_10.size() && iConfig < fHistEEC_MC_5_10[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_5_10[i][iConfig].size())
      hEEC_MC_5_10 = fHistEEC_MC_5_10[i][iConfig][ciTree];

    TH1D* hEEC_MC_10_15 = 0;
    if (i < fHistEEC_MC_10_15.size() && iConfig < fHistEEC_MC_10_15[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_10_15[i][iConfig].size())
      hEEC_MC_10_15 = fHistEEC_MC_10_15[i][iConfig][ciTree];

    TH1D* hEEC_MC_15_20 = 0;
    if (i < fHistEEC_MC_15_20.size() && iConfig < fHistEEC_MC_15_20[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_15_20[i][iConfig].size())
      hEEC_MC_15_20 = fHistEEC_MC_15_20[i][iConfig][ciTree];

    TH1D* hEEC_MC_20_30 = 0;
    if (i < fHistEEC_MC_20_30.size() && iConfig < fHistEEC_MC_20_30[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_20_30[i][iConfig].size())
      hEEC_MC_20_30 = fHistEEC_MC_20_30[i][iConfig][ciTree];

    TH1D* hEEC_MC_30_50 = 0;
    if (i < fHistEEC_MC_30_50.size() && iConfig < fHistEEC_MC_30_50[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_30_50[i][iConfig].size())
      hEEC_MC_30_50 = fHistEEC_MC_30_50[i][iConfig][ciTree];

    TH1D* hEEC_MC_50_100 = 0;
    if (i < fHistEEC_MC_50_100.size() && iConfig < fHistEEC_MC_50_100[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_MC_50_100[i][iConfig].size())
      hEEC_MC_50_100 = fHistEEC_MC_50_100[i][iConfig][ciTree];
    
      
    TH1D* hEEC_unmatched_5_10 = 0;
    if (i < fHistEEC_unmatched_5_10.size() && iConfig < fHistEEC_unmatched_5_10[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_5_10[i][iConfig].size())
      hEEC_unmatched_5_10 = fHistEEC_unmatched_5_10[i][iConfig][ciTree];
    
    TH1D* hEEC_unmatched_10_15 = 0;
    if (i < fHistEEC_unmatched_10_15.size() && iConfig < fHistEEC_unmatched_10_15[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_10_15[i][iConfig].size())
      hEEC_unmatched_10_15 = fHistEEC_unmatched_10_15[i][iConfig][ciTree];

    TH1D* hEEC_unmatched_15_20 = 0;
    if (i < fHistEEC_unmatched_15_20.size() && iConfig < fHistEEC_unmatched_15_20[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_15_20[i][iConfig].size())
      hEEC_unmatched_15_20 = fHistEEC_unmatched_15_20[i][iConfig][ciTree];

    TH1D* hEEC_unmatched_20_30 = 0;
    if (i < fHistEEC_unmatched_20_30.size() && iConfig < fHistEEC_unmatched_20_30[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_20_30[i][iConfig].size())
      hEEC_unmatched_20_30 = fHistEEC_unmatched_20_30[i][iConfig][ciTree];

    TH1D* hEEC_unmatched_30_50 = 0;
    if (i < fHistEEC_unmatched_30_50.size() && iConfig < fHistEEC_unmatched_30_50[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_30_50[i][iConfig].size())
      hEEC_unmatched_30_50 = fHistEEC_unmatched_30_50[i][iConfig][ciTree];

    TH1D* hEEC_unmatched_50_100 = 0;
    if (i < fHistEEC_unmatched_50_100.size() && iConfig < fHistEEC_unmatched_50_100[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_50_100[i][iConfig].size())
      hEEC_unmatched_50_100 = fHistEEC_unmatched_50_100[i][iConfig][ciTree];

    TH1D* hEEC_unmatched_all = 0;
    if (i < fHistEEC_unmatched_all.size() && iConfig < fHistEEC_unmatched_all[i].size() && ciTree >= 0 && ciTree < (int)fHistEEC_unmatched_all[i][iConfig].size())
      hEEC_unmatched_all = fHistEEC_unmatched_all[i][iConfig][ciTree];
    

    //Create EEC histograms for embedding for R, config, centrality class
  
    std::vector<MyJet> myRecoJets;
    if (!vetoReco && !tracks.empty()) {
      fastjet::ClusterSequenceArea reco_cluster_seq(tracks, jet_def, area_def); //comment for full jets, use jetTracks for clustering to be consistent with rho estimation
      std::vector<fastjet::PseudoJet> fjets_all =
          sorted_by_pt(reco_cluster_seq.inclusive_jets(fJetPtMin));
    
      fastjet::Selector fiducial_cut_selector = fastjet::SelectorAbsEtaMax(maxRapJet);
      std::vector<fastjet::PseudoJet> RecoJets = fiducial_cut_selector(fjets_all);
      myRecoJets.reserve(RecoJets.size());
      for (size_t j = 0; j < RecoJets.size(); ++j) {
        myRecoJets.push_back(MyJet(RecoJets[j], rho));
        if (!mIsEmbedding){
          auto constituents = RecoJets[j].constituents();
          for (const auto& c : constituents) {
            if (c.pt() < 0.1) continue; // Skip very low pT constituents
            if (c.is_pure_ghost()) continue; //Keeps only real constituents, skip ghosts
            c_runid = fRunNumber;
            c_eventid = fEventId;
            c_ijet = j;
            c_px = c.px();
            c_py = c.py();
            c_pz = c.pz();
            c_pt = c.perp();
            c_E = c.E();
            c_eta = c.eta();
            c_phi = c.phi();
            c_charge = c.user_index(); // Use user_index to store charge information

            constituentTree->Fill();

            if (c_charge != 0){
              phi_vector.push_back(c_phi);
              eta_vector.push_back(c_eta);
              pt_vector.push_back(c_pt);
              E_vector.push_back(c_E);
            }; // Skip neutral constituents
          }
      
          //Calculate EEC and RL
          for (size_t h = 0; h < phi_vector.size(); ++h) {
            for (size_t k = h+1; k < phi_vector.size(); ++k) {
              if (h == k) continue; // Don't compare with itself
              // Calculate RL
              eec_ijet = j;
              eec_eventid = fEventId;
              eec_runid = fRunNumber;
              double delta_phi = phi_vector[h] - phi_vector[k];
              double delta_eta = eta_vector[h] - eta_vector[k];

              RL = sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
              if (RL > 0.8) {
                //cout << "BEFORE____RL: " << RL << "; dEta: " << delta_eta << "; dPhi: " << delta_phi << endl;
                if (delta_phi <= -TMath::Pi()) {
                  delta_phi = delta_phi + TMath::TwoPi();
                  RL = sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
                }
                if (delta_phi >= TMath::Pi()) { 
                  delta_phi = delta_phi - TMath::TwoPi();
                  RL = sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
                }
                //cout << "AFTER_____RL: " << RL << "; dEta: " << delta_eta << "; dPhi: " << delta_phi << endl;
              }
              //eec_data = E_vector[h] * E_vector[k] / (myRecoJets[j].pt_corr * myRecoJets[j].pt_corr);
              eec_data = pt_vector[h] * pt_vector[k] / (myRecoJets[j].pt_corr * myRecoJets[j].pt_corr);
              try_EECTree->Fill();


              //hEEC->Fill(RL, eec_data); // Fill histogram with RL and EEC value
        
              if (myRecoJets[j].pt_corr >= 5 && myRecoJets[j].pt_corr < 10){
                hEEC_5_10->Fill(RL, eec_data); // Fill histogram for 5-10 GeV/c jets
              }
              if (myRecoJets[j].pt_corr >= 10 && myRecoJets[j].pt_corr < 15){
                hEEC_10_15->Fill(RL, eec_data); // Fill histogram for 10-15 GeV/c jets
              }
              if (myRecoJets[j].pt_corr>= 15 && myRecoJets[j].pt_corr < 20){
                hEEC_15_20->Fill(RL, eec_data); // Fill histogram for 15-20 GeV/c jets
                //EEC_low = EEC_low + eec_data;
              }
              if (myRecoJets[j].pt_corr>= 20 && myRecoJets[j].pt_corr < 30){
                hEEC_20_30->Fill(RL, eec_data); // Fill histogram for 20-30 GeV/c jets
                //EEC_mid = EEC_mid + eec_data;
              }
              if (myRecoJets[j].pt_corr >= 30 && myRecoJets[j].pt_corr < 50){
                hEEC_30_50->Fill(RL, eec_data); // Fill histogram for 30-50 GeV/c jets
                //EEC_high = EEC_high + eec_data;
              }
            
            }
          }

          phi_vector.clear();
          eta_vector.clear();
          E_vector.clear();
          pt_vector.clear();
        }
      }
    } 


  if (ptMaxVeto > 0.0 && !vetoReco && !tracks.empty()) {
  for (size_t jr = 0; jr < myRecoJets.size(); ++jr) {
    if (myRecoJets[jr].pt_corr > ptMaxVeto) {
      return kStOK; // veto whole event
    }
  }
}

  //======================================================================//

  // pick tree for this (R, class)

  //==================== Embedding mode ==================//
  if (mIsEmbedding) {
    //============================== MC jets ===============================//
    fastjet::ClusterSequenceArea mc_cluster_seq(MCjetTracks, jet_def, area_def);
    vector<fastjet::PseudoJet> Mcjets_all =
        sorted_by_pt(mc_cluster_seq.inclusive_jets(1.0));

    fastjet::Selector McFiducial_cut_selector =
        fastjet::SelectorAbsEtaMax(maxRapJet) *
        fastjet::SelectorPtMin(0.01);

    vector<fastjet::PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
    vector<MyJet> myMcJets;
    myMcJets.reserve(McJets.size());

    for (size_t j = 0; j < McJets.size(); ++j) {
      myMcJets.push_back(MyJet(McJets[j], 0.0f)); // rho = 0 for truth
    }

    for (size_t j = 0; j < myMcJets.size(); ++j) {
      fMcJet = myMcJets[j];
      if (MCJetTree) MCJetTree->Fill();
    }
    
    //========================= MC–Reco matching ===========================//
    vector<MatchedJetPair> MatchedJets =
        MatchJetsEtaPhi(myMcJets, myRecoJets, fR[i]);

    for (size_t im = 0; im < MatchedJets.size(); ++im) {
      const MatchedJetPair& mp = MatchedJets[im];
      fMcJet   = mp.first;
      fRecoJet = mp.second;
      fDeltaR  = fMcJet.deltaR(fRecoJet);

      const bool haveReco = (fRecoJet.pt >= 0.0);   // with default -999 this is false if no reco
      const bool haveMC   = (fMcJet.pt   >= 0.0);

      if (jetTree && (haveMC || haveReco)) {
          jetTree->Fill();
      }

      //Constituent-level studies for matched jets + EEC calculation
      if (haveMC && haveReco) {
        
        for(size_t k = 0; k < fMcJet.constituents_pt.size(); ++k){
          if (fMcJet.constituents_pt[k] < 0.1) continue; // Skip very low pT constituents
          // Fill MC constituent tree variables here using fMcJet.constituents_px[k], fMcJet.constituents_py[k], etc.

          c_mc_runid = fRunNumber;
          c_mc_eventid = fEventId;
          c_mc_ijet = im;
          c_mc_px = fMcJet.constituents_px[k];
          c_mc_py = fMcJet.constituents_py[k];
          c_mc_pz = fMcJet.constituents_pz[k];
          c_mc_pt = fMcJet.constituents_pt[k];
          c_mc_E = fMcJet.constituents_E[k];
          c_mc_eta = fMcJet.constituents_eta[k];
          c_mc_phi = fMcJet.constituents_phi[k];
          c_mc_charge = fMcJet.constituents_charge[k];
          MCConstituentTree->Fill();

          if (c_mc_charge != 0){
            mc_phi_vector.push_back(c_mc_phi);
            mc_eta_vector.push_back(c_mc_eta);
            mc_pt_vector.push_back(c_mc_pt);
            mc_E_vector.push_back(c_mc_E);
          }; // Skip neutral constituents
        }

        for(size_t k = 0; k < fRecoJet.constituents_pt.size(); ++k){
          if (fRecoJet.constituents_pt[k] < 0.1) continue; // Skip very low pT constituents
          // Fill matched constituent tree variables here using fRecoJet.constituents_px[k], fRecoJet.constituents_py[k], etc.

          c_matched_runid = fRunNumber;
          c_matched_eventid = fEventId;
          c_matched_ijet = im;
          c_matched_px = fRecoJet.constituents_px[k];
          c_matched_py = fRecoJet.constituents_py[k];
          c_matched_pz = fRecoJet.constituents_pz[k];
          c_matched_pt = fRecoJet.constituents_pt[k];
          c_matched_E = fRecoJet.constituents_E[k];
          c_matched_eta = fRecoJet.constituents_eta[k];
          c_matched_phi = fRecoJet.constituents_phi[k];
          c_matched_charge = fRecoJet.constituents_charge[k];

          MatchedConstituentTree->Fill();

          if(c_matched_charge != 0){
            matched_phi_vector.push_back(c_matched_phi);
            matched_eta_vector.push_back(c_matched_eta);
            matched_pt_vector.push_back(c_matched_pt);
            matched_E_vector.push_back(c_matched_E);
          }; // Skip neutral constituents
        }

        for (size_t k = 0; k < mc_phi_vector.size(); ++k) {
          for(size_t l = k+1; l<mc_phi_vector.size(); ++l){
            // Calculate RL for MC constituents
            eec_mc_ijet = im;
            eec_mc_eventid = fEventId;
            eec_mc_runid = fRunNumber;
            
            double mc_delta_phi = mc_phi_vector[k] - mc_phi_vector[l];
            double mc_delta_eta = mc_eta_vector[k] - mc_eta_vector[l];
            RL_mc = sqrt(mc_delta_phi * mc_delta_phi + mc_delta_eta * mc_delta_eta);
            if (RL_mc > 0.8) {
              if (mc_delta_phi <= -TMath::Pi()) {
                mc_delta_phi = mc_delta_phi + TMath::TwoPi();
                RL_mc = sqrt(mc_delta_phi * mc_delta_phi + mc_delta_eta * mc_delta_eta);
              }
              if (mc_delta_phi >= TMath::Pi()) { 
                mc_delta_phi = mc_delta_phi - TMath::TwoPi();
                RL_mc = sqrt(mc_delta_phi * mc_delta_phi + mc_delta_eta * mc_delta_eta);
              }
            }
            eec_mc = mc_pt_vector[k] * mc_pt_vector[l] / (fMcJet.pt * fMcJet.pt);
            // Fill MC EEC tree or histogram here using RL_mc and eec_mc
            EECTree_MC->Fill();

            if (fMcJet.pt >= 5 && fMcJet.pt < 10){
              hEEC_MC_5_10->Fill(RL_mc, eec_mc); // Fill histogram for 5-10 GeV/c jets
            }

            if (fMcJet.pt >= 10 && fMcJet.pt < 15){
              hEEC_MC_10_15->Fill(RL_mc, eec_mc); // Fill histogram for 10-15 GeV/c jets
            }

            if (fMcJet.pt >= 15 && fMcJet.pt < 20){
              hEEC_MC_15_20->Fill(RL_mc, eec_mc); // Fill histogram for 15-20 GeV/c jets
            }

            if (fMcJet.pt >= 20 && fMcJet.pt < 30){
              hEEC_MC_20_30->Fill(RL_mc, eec_mc); // Fill histogram for 20-30 GeV/c jets
            }

            if (fMcJet.pt >= 30 && fMcJet.pt < 50){
              hEEC_MC_30_50->Fill(RL_mc, eec_mc); // Fill histogram for 30-50 GeV/c jets
            }

            if (fMcJet.pt >= 50 && fMcJet.pt < 100){
              hEEC_MC_50_100->Fill(RL_mc, eec_mc); // Fill histogram for 50-100 GeV/c jets
            }
          }
        }

        for(size_t k = 0; k < matched_phi_vector.size(); ++k){
          for(size_t l = k+1; l<matched_phi_vector.size(); ++l){
            // Calculate RL for matched constituents
            eec_matched_ijet = im;
            eec_matched_eventid = fEventId;
            eec_matched_runid = fRunNumber;
            
            double matched_delta_phi = matched_phi_vector[k] - matched_phi_vector[l];
            double matched_delta_eta = matched_eta_vector[k] - matched_eta_vector[l];
            RL_matched = sqrt(matched_delta_phi * matched_delta_phi + matched_delta_eta * matched_delta_eta);
            if (RL_matched > 0.8) {
              if (matched_delta_phi <= -TMath::Pi()) {
                matched_delta_phi = matched_delta_phi + TMath::TwoPi();
                RL_matched = sqrt(matched_delta_phi * matched_delta_phi + matched_delta_eta * matched_delta_eta);
              }
              if (matched_delta_phi >= TMath::Pi()) { 
                matched_delta_phi = matched_delta_phi - TMath::TwoPi();
                RL_matched = sqrt(matched_delta_phi * matched_delta_phi + matched_delta_eta * matched_delta_eta);
              }
            }
            eec_matched = matched_pt_vector[k] * matched_pt_vector[l] / (fRecoJet.pt_corr * fRecoJet.pt_corr);
            // Fill Matched EEC tree or histogram here using RL_matched and eec_matched
            EECTree_matched->Fill();

            if (fRecoJet.pt_corr >= 5 && fRecoJet.pt_corr < 10){
              hEEC_matched_5_10->Fill(RL_matched, eec_matched); // Fill histogram for 5-10 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 10 && fRecoJet.pt_corr < 15){
              hEEC_matched_10_15->Fill(RL_matched, eec_matched); // Fill histogram for 10-15 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 15 && fRecoJet.pt_corr < 20){
              hEEC_matched_15_20->Fill(RL_matched, eec_matched); // Fill histogram for 15-20 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 20 && fRecoJet.pt_corr < 30){
              hEEC_matched_20_30->Fill(RL_matched, eec_matched); // Fill histogram for 20-30 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 30 && fRecoJet.pt_corr < 50){
              hEEC_matched_30_50->Fill(RL_matched, eec_matched); // Fill histogram for 30-50 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 50 && fRecoJet.pt_corr < 100){
              hEEC_matched_50_100->Fill(RL_matched, eec_matched); // Fill histogram for 50-100 GeV/c jets
            }
          }
        }
        mc_phi_vector.clear();
        mc_eta_vector.clear();
        mc_pt_vector.clear();
        mc_E_vector.clear();

        matched_phi_vector.clear();
        matched_eta_vector.clear();
        matched_pt_vector.clear();
        matched_E_vector.clear(); 
      }
      if (!haveMC && haveReco){
        for (size_t k = 0; k < fRecoJet.constituents_pt.size(); ++k) {
          if (fRecoJet.constituents_pt[k] < 0.1) continue; // Skip very low pT constituents
          // Fill unmatched constituent tree variables here using fRecoJet.constituents_px[k], fRecoJet.constituents_py[k], etc.

          c_unmatched_runid = fRunNumber;
          c_unmatched_eventid = fEventId;
          c_unmatched_ijet = im;
          c_unmatched_px = fRecoJet.constituents_px[k];
          c_unmatched_py = fRecoJet.constituents_py[k];
          c_unmatched_pz = fRecoJet.constituents_pz[k];
          c_unmatched_pt = fRecoJet.constituents_pt[k];
          c_unmatched_E = fRecoJet.constituents_E[k];
          c_unmatched_eta = fRecoJet.constituents_eta[k];
          c_unmatched_phi = fRecoJet.constituents_phi[k];
          c_unmatched_charge = fRecoJet.constituents_charge[k];

          UnmatchedConstituentTree->Fill();

          if(c_unmatched_charge != 0){
            unmatched_phi_vector.push_back(c_unmatched_phi);
            unmatched_eta_vector.push_back(c_unmatched_eta);
            unmatched_pt_vector.push_back(c_unmatched_pt);
            unmatched_E_vector.push_back(c_unmatched_E);
          }; // Skip neutral constituents
        }
        
        for (size_t k = 0; k < unmatched_phi_vector.size(); ++k) {
          for (size_t l = k+1; l < unmatched_phi_vector.size(); ++l) {
            // Calculate RL for unmatched constituents
            eec_unmatched_ijet = im;
            eec_unmatched_eventid = fEventId;
            eec_unmatched_runid = fRunNumber;

            double unmatched_delta_phi = unmatched_phi_vector[k] - unmatched_phi_vector[l];
            double unmatched_delta_eta = unmatched_eta_vector[k] - unmatched_eta_vector[l];
            RL_unmatched = sqrt(unmatched_delta_phi * unmatched_delta_phi + unmatched_delta_eta * unmatched_delta_eta);
            if (RL_unmatched > 0.8) {
              if (unmatched_delta_phi <= -TMath::Pi()) {
                unmatched_delta_phi = unmatched_delta_phi + TMath::TwoPi();
                RL_unmatched = sqrt(unmatched_delta_phi * unmatched_delta_phi + unmatched_delta_eta * unmatched_delta_eta);
              }
              if (unmatched_delta_phi >= TMath::Pi()) { 
                unmatched_delta_phi = unmatched_delta_phi - TMath::TwoPi();
                RL_unmatched = sqrt(unmatched_delta_phi * unmatched_delta_phi + unmatched_delta_eta * unmatched_delta_eta);
              }
            }
            eec_unmatched = unmatched_pt_vector[k] * unmatched_pt_vector[l] / (fRecoJet.pt_corr * fRecoJet.pt_corr);
            // Fill Unmatched EEC tree or histogram here using RL_unmatched and eec_unmatched
            // You can create a separate tree or histogram for unmatched jets if needed


            EECTree_unmatched->Fill();
            hEEC_unmatched_all->Fill(RL_unmatched, eec_unmatched); // Fill histogram for all unmatched jets

            if (fRecoJet.pt_corr >= 5 && fRecoJet.pt_corr < 10){
              hEEC_unmatched_5_10->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 5-10 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 10 && fRecoJet.pt_corr < 15){
              hEEC_unmatched_10_15->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 10-15 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 15 && fRecoJet.pt_corr < 20){
              hEEC_unmatched_15_20->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 15-20 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 20 && fRecoJet.pt_corr < 30){
              hEEC_unmatched_20_30->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 20-30 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 30 && fRecoJet.pt_corr < 50){
              hEEC_unmatched_30_50->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 30-50 GeV/c jets
            }

            if (fRecoJet.pt_corr >= 50 && fRecoJet.pt_corr < 100){
              hEEC_unmatched_50_100->Fill(RL_unmatched, eec_unmatched); // Fill histogram for 50-100 GeV/c jets
            }
          }
        }

        unmatched_phi_vector.clear();
        unmatched_eta_vector.clear();
        unmatched_pt_vector.clear();
        unmatched_E_vector.clear();
      }
      

    } // end loop over MatchedJets

  } else {
    //==================== Data mode ===========================//
    if (!vetoReco) {
      for (size_t jr = 0; jr < myRecoJets.size(); ++jr) {
        fRecoJet = myRecoJets[jr];
        fMcJet   = MyJet();   // dummy
        fDeltaR  = -1.0;

        if (!jetTree) continue;
          jetTree->Fill(); 
      }     
    } // if (!vetoReco)
  } // end embedding/data
  
  } // end loop over R
  } // end loop over configs
  return kStOK;
}



//-----------------------------------------------------------------------------
////Correct tower energy
//-----------------------------------------------------------------------------
Double_t StPicoHFJetMaker::GetTowerCalibEnergy(Int_t TowerId) {
  StPicoBTowHit *tower =
      static_cast<StPicoBTowHit *>(mPicoDst->btowHit(TowerId - 1));
  Float_t pedestal, rms;
  Int_t status;
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);
  Double_t *TowerCoeff;
  if (fRunNumber <= 15094020)
    TowerCoeff = CPre;
  else
    TowerCoeff = CLowMidHigh;
  Double_t calibEnergy = TowerCoeff[TowerId - 1] * (tower->adc() - pedestal);
  return calibEnergy;
}

//-----------------------------------------------------------------------------
////Correct tower eta for Vz position //// Not used anymore
//-----------------------------------------------------------------------------
Double_t StPicoHFJetMaker::vertexCorrectedEta(double eta, double vz) {
  double tower_theta = 2.0 * atan(exp(-eta));
  double z = 0.0;
  if (eta != 0.0)
    z = mBarrelRadius / tan(tower_theta);
  double z_diff = z - vz;
  double theta_corr = atan2(mBarrelRadius, z_diff);
  double eta_corr = -log(tan(theta_corr / 2.0));
  return eta_corr;
}

//-----------------------------------------------------------------------------
// Fill array with momentum of BEMC-matched tracks
//-----------------------------------------------------------------------------
Bool_t StPicoHFJetMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst,
                                              TVector3 mPrimVtx) {
  // loop over global tracks  - towers
  UInt_t nTracks = mPicoDst->numberOfTracks();
  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
    StPicoTrack *trk = mPicoDst->track(itrack);
    TVector3 gMom = trk->gMom();
    // using global tracks
    double pT = gMom.Perp();
    if (pT != pT || pT < 0.2)
      continue;
    float eta = gMom.PseudoRapidity();
    if (fabs(eta) > 1)
      continue;
  //  float phi = gMom.Phi();

    float nHitsFit = trk->nHitsFit();
    float nHitsMax = trk->nHitsMax();
    if (nHitsFit < 15 || nHitsFit / nHitsMax < 0.52)
      continue; // some basic QA cuts
    double Bfield = mPicoDst->event()->bField();

    StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
  //  float vtx_x = mPrimVtx.x();
  //  float vtx_y = mPrimVtx.y();
  //  float vtx_z = mPrimVtx.z();

    float dca_z = abs(trk->gDCAz(mPicoDst->event()->primaryVertex().z()));
    if (fabs(dca_z) > maxdcazhadroncorr)
      continue;
    int TowIndex = -99999;
    TowIndex = trk->bemcTowerIndex();
    float p = 0;
    if (TowIndex >= 0) {
      p = gMom.Mag();
      Sump[TowIndex] += p;
    }
  } // END global track loop
  return true;
}

//-----------------------------------------------------------------------------
// Jet matching function using only eta-phi criteria
//-----------------------------------------------------------------------------
vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R) {
  const double matchRadius = 0.6*R; 
  vector<char> recoUsed(RecoJets.size(), 0);
  vector<MatchedJetPair> matchedJets;
  matchedJets.reserve(McJets.size() + RecoJets.size());

  // For each MC jet, find the closest unused reco jet within matchRadius
  for (const auto &mcJet : McJets) {
    int bestIdx = -1;
    double bestDr = matchRadius;
    for (size_t j = 0; j < RecoJets.size(); ++j) {
      if (recoUsed[j]) continue;
      double dr = mcJet.deltaR(RecoJets[j]); 
      if (dr < bestDr) {
        bestDr = dr;
        bestIdx = (int)j;
      }
    }
    if (bestIdx >= 0) {
      matchedJets.emplace_back(mcJet, RecoJets[bestIdx]);
      recoUsed[bestIdx] = 1;
    } else {
      matchedJets.emplace_back(mcJet, MyJet()); // unmatched MC -> dummy reco
    }
  }

  // Append unmatched reco jets
  for (size_t j = 0; j < RecoJets.size(); ++j) {
    if (!recoUsed[j]) matchedJets.emplace_back(MyJet(), RecoJets[j]);
  }

  return matchedJets;
}