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

    scaleF_POT_BNB  = 1;
    scaleF_POT_NuMI = 1;
    
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
    
    dm2_41 = 0;
    sin2_2theta_14 = 0;
    sin2_theta_24  = 0;
    sin2_theta_34  = 0;

    minimization_status     = -1;
    minimization_chi2       = -1;
    minimization_dm2_41_val = -1;
    minimization_sin2_2theta_14_val = -1;
    minimization_sin2_theta_24_val  = -1;
    minimization_sin2_theta_34_val  = -1;
    minimization_dm2_41_err         = -1;
    minimization_sin2_2theta_14_err = -1;
    minimization_sin2_theta_24_err  = -1;
    minimization_sin2_theta_34_err  = -1;
  }

  ////////////////////////////////////////////////////// data members
  
  TMatrixD matrix_transform;
  TMatrixD matrix_transform_wiPOT;

  double scaleF_POT_BNB;
  double scaleF_POT_NuMI;
  
  map<int, TH1D*>map_default_h1d_meas;
  map<int, int>map_default_h1d_meas_bins;
  map<int, double>map_default_h1d_meas_xlow;
  map<int, double>map_default_h1d_meas_xhgh;  
  vector<double> vector_default_newworld_meas;
  TMatrixD matrix_default_newworld_meas;
  TMatrixD matrix_default_newworld_meas_POT;

  map<int, TH1D*>map_default_h1d_pred;
  map<int, int>map_default_h1d_pred_bins;
  map<int, double>map_default_h1d_pred_xlow;
  map<int, double>map_default_h1d_pred_xhgh;  
  vector<double> vector_default_oldworld_pred;
  TMatrixD matrix_default_oldworld_pred;
  TMatrixD matrix_default_newworld_pred;

  TMatrixD matrix_default_oldworld_abs_syst_addi;  // for dirt additional syst, approximation: always use the same absolute cov
  TMatrixD matrix_default_oldworld_abs_syst_addi_POT;
  TMatrixD matrix_default_oldworld_abs_syst_mcstat;// only newworld, iteration
  TMatrixD matrix_default_oldworld_abs_syst_flux;
  TMatrixD matrix_default_oldworld_abs_syst_geant;
  TMatrixD matrix_default_oldworld_abs_syst_Xs;
  TMatrixD matrix_default_oldworld_abs_syst_det;
  
  TMatrixD matrix_default_oldworld_rel_syst_addi;  //
  TMatrixD matrix_default_oldworld_rel_syst_mcstat;//
  TMatrixD matrix_default_oldworld_rel_syst_flux;  // initialized
  TMatrixD matrix_default_oldworld_rel_syst_geant; // initialized
  TMatrixD matrix_default_oldworld_rel_syst_Xs;    // initialized
  TMatrixD matrix_default_oldworld_rel_syst_det;   // initialized
  
  TMatrixD matrix_default_newworld_abs_syst_addi;
  TMatrixD matrix_default_newworld_abs_syst_mcstat;// only newworld
  TMatrixD matrix_default_newworld_abs_syst_mcstat_POT;
  TMatrixD matrix_default_newworld_abs_syst_flux;
  TMatrixD matrix_default_newworld_abs_syst_geant;
  TMatrixD matrix_default_newworld_abs_syst_Xs;
  TMatrixD matrix_default_newworld_abs_syst_det;

  TMatrixD matrix_default_newworld_rel_syst_addi;
  TMatrixD matrix_default_newworld_rel_syst_mcstat;
  TMatrixD matrix_default_newworld_rel_syst_flux;
  TMatrixD matrix_default_newworld_rel_syst_geant;
  TMatrixD matrix_default_newworld_rel_syst_Xs;
  TMatrixD matrix_default_newworld_rel_syst_det;

  TMatrixD matrix_default_newworld_abs_syst_total;
  
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
  TMatrixD matrix_oscillation_effect_result_oldworld;

  vector<double>vector_oscillation_base_NuMI_nueCC_scaleFPOT; vector< vector<EventInfo> >vector_vector_oscillation_base_NuMI_nueCC_info;
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
  
  double dm2_41;
  double sin2_2theta_14;
  double sin2_theta_24;
  double sin2_theta_34;

  TMatrixD matrix_meas2fitdata_newworld;
  
  int    minimization_status;
  double minimization_chi2;
  double minimization_dm2_41_val;
  double minimization_sin2_2theta_14_val;
  double minimization_sin2_theta_24_val;
  double minimization_sin2_theta_34_val;
  double minimization_dm2_41_err;
  double minimization_sin2_2theta_14_err;
  double minimization_sin2_theta_24_err;
  double minimization_sin2_theta_34_err;
  
  ////////////////////////////////////////////////////// member functions

  void Set_default_cv_cov(TString default_cv_file, TString default_dirtadd_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir);
 
  void Set_oscillation_base();
  void Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, vector< vector<EventInfo> > *vec_vec_eventinfo);

  /// matrix_oscillation_effect_result_oldworld = matrix_oscillation_base + matrix_oscillation_effect;

  void Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34) {
    dm2_41 = val_dm2_41;
    sin2_2theta_14 = val_sin2_2theta_14;
    sin2_theta_24  = val_sin2_theta_24;
    sin2_theta_34  = val_sin2_theta_34;
  }
  
  void Apply_oscillation();

  double Prob_oscillaion(double Etrue, double baseline, TString strflag_osc);
  
  void Set_apply_POT();

  void Set_meas2fitdata() {
    //matrix_meas2fitdata_newworld = matrix_default_newworld_meas_POT;
    matrix_meas2fitdata_newworld = matrix_default_newworld_pred;
  }

  int Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar);
  
  ////////////////////////////////////////////////////// data members

 private:
  int default_oldworld_rows;
  int default_newworld_rows;
  
};

#endif
