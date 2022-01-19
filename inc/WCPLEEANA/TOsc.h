#ifndef TOsc_ana
#define TOsc_ana

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

struct EventInfo {
  int           e2e_pdg;
  int           e2e_flag_FC;
  double        e2e_Etrue;
  double        e2e_Ereco;
  double        e2e_weight_xs;
  double        e2e_baseline;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////// TOsc

class TOsc {
 public:
  TOsc() {
    cout<<endl<<" ---> Hello TOsc"<<endl<<endl;

    default_oldworld_rows = 0;
    default_newworld_rows = 0;
  
    flag_NuMI_nue2nue   = 0;
    flag_NuMI_numu2numu = 0;
    flag_NuMI_numu2nue  = 0;  
    flag_NuMI_nue2numu  = 0;
    flag_NuMI_NC_1minus_nue2sterile  = 0;
    flag_NuMI_NC_1minus_numu2sterile = 0;
  
    flag_BNB_nue2nue   = 0;
    flag_BNB_numu2numu = 0;
    flag_BNB_numu2nue  = 0;  
    flag_BNB_nue2numu  = 0;
    flag_BNB_NC_1minus_nue2sterile  = 0;
    flag_BNB_NC_1minus_numu2sterile = 0;
    
  }

  ////////////////////////////////////////////////////// data members

  TMatrixD matrix_transform;

  map<int, TH1D*>map_default_h1d_meas;
  map<int, int>map_default_h1d_meas_bins;
  map<int, double>map_default_h1d_meas_xlow;
  map<int, double>map_default_h1d_meas_xhgh;  
  vector<double> vector_default_newworld_meas;
  TMatrixD matrix_default_newworld_meas;

  map<int, TH1D*>map_default_h1d_pred;
  map<int, int>map_default_h1d_pred_bins;
  map<int, double>map_default_h1d_pred_xlow;
  map<int, double>map_default_h1d_pred_xhgh;  
  vector<double> vector_default_oldworld_pred;
  TMatrixD matrix_default_oldworld_pred;
  TMatrixD matrix_default_newworld_pred;

  TMatrixD matrix_default_oldworld_abs_syst_addi;// for dirt additional syst
  TMatrixD matrix_default_oldworld_abs_syst_mcstat;// only newworld
  TMatrixD matrix_default_oldworld_abs_syst_flux;
  TMatrixD matrix_default_oldworld_abs_syst_geant;
  TMatrixD matrix_default_oldworld_abs_syst_Xs;
  TMatrixD matrix_default_oldworld_abs_syst_det;
  
  TMatrixD matrix_default_newworld_abs_syst_addi;
  TMatrixD matrix_default_newworld_abs_syst_mcstat;// only newworld
  TMatrixD matrix_default_newworld_abs_syst_flux;
  TMatrixD matrix_default_newworld_abs_syst_geant;
  TMatrixD matrix_default_newworld_abs_syst_Xs;
  TMatrixD matrix_default_newworld_abs_syst_det;

  ///////
  
  bool flag_NuMI_nue2nue;
  bool flag_NuMI_numu2numu;
  bool flag_NuMI_numu2nue;  
  bool flag_NuMI_nue2numu;
  bool flag_NuMI_NC_1minus_nue2sterile;
  bool flag_NuMI_NC_1minus_numu2sterile;
  
  bool flag_BNB_nue2nue;
  bool flag_BNB_numu2numu;
  bool flag_BNB_numu2nue;  
  bool flag_BNB_nue2numu;
  bool flag_BNB_NC_1minus_nue2sterile;
  bool flag_BNB_NC_1minus_numu2sterile;
  
  ///////
  TMatrixD matrix_oscillation_base;
  
  vector<double>vector_oscillation_base_NuMI_nueCC_scaleFPOT; vector< vector<EventInfo> >vector_vector_oscillation_base_NuMI_nueCC;
  vector<double>vector_oscillation_base_NuMI_numuCC_scaleFPOT;
  vector<double>vector_oscillation_base_NuMI_nueNC_scaleFPOT;
  vector<double>vector_oscillation_base_NuMI_numuNC_scaleFPOT;
  
  vector<double>vector_oscillation_base_NuMIdirt_nueCC_scaleFPOT;
  vector<double>vector_oscillation_base_NuMIdirt_numuCC_scaleFPOT;
  vector<double>vector_oscillation_base_NuMIdirt_nueNC_scaleFPOT;
  vector<double>vector_oscillation_base_NuMIdirt_numuNC_scaleFPOT;
  
  vector<double>vector_oscillation_base_BNB_nueCC_scaleFPOT;
  vector<double>vector_oscillation_base_BNB_numuCC_scaleFPOT;
  vector<double>vector_oscillation_base_BNB_nueNC_scaleFPOT;
  vector<double>vector_oscillation_base_BNB_numuNC_scaleFPOT;
  
  vector<double>vector_oscillation_base_BNBdirt_nueCC_scaleFPOT;
  vector<double>vector_oscillation_base_BNBdirt_numuCC_scaleFPOT;
  vector<double>vector_oscillation_base_BNBdirt_nueNC_scaleFPOT;
  vector<double>vector_oscillation_base_BNBdirt_numuNC_scaleFPOT;
  
  ////////////////////////////////////////////////////// member functions

  void Set_default_cv_cov(TString default_cv_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir);

  void Set_oscillation_base();
  void Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, vector< vector<EventInfo> > *vec_vec_eventinfo);
  
  ////////////////////////////////////////////////////// data members

 private:
  int default_oldworld_rows;
  int default_newworld_rows;
  
};

#endif
