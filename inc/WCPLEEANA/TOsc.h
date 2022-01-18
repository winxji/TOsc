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

///////////////////////////////////////////////////////////////////////////////////////////////////////// TOsc

class TOsc {
 public:
  TOsc() {
    cout<<endl<<" ---> Hello TOsc"<<endl<<endl;

    default_oldworld_rows = 0;
    default_newworld_rows = 0;
  }

  ////////////////////////////////////////////////////// data members

  TMatrixD matrix_transform;

  map<int, TH1D*>map_default_h1d_meas;
  vector<double> vector_newworld_meas;
  TMatrixD matrix_newworld_meas;

  map<int, TH1D*>map_default_h1d_pred;
  vector<double> vector_oldworld_pred;
  TMatrixD matrix_oldworld_pred;
  TMatrixD matrix_newworld_pred;

  TMatrixD matrix_default_oldworld_abs_syst_addi;
  TMatrixD matrix_default_oldworld_abs_syst_mcstat;
  TMatrixD matrix_default_oldworld_abs_syst_flux;
  TMatrixD matrix_default_oldworld_abs_syst_geant;
  TMatrixD matrix_default_oldworld_abs_syst_Xs;
  TMatrixD matrix_default_oldworld_abs_syst_det;
  
  TMatrixD matrix_default_newworld_abs_syst_addi;
  TMatrixD matrix_default_newworld_abs_syst_mcstat;
  TMatrixD matrix_default_newworld_abs_syst_flux;
  TMatrixD matrix_default_newworld_abs_syst_geant;
  TMatrixD matrix_default_newworld_abs_syst_Xs;
  TMatrixD matrix_default_newworld_abs_syst_det;
  
  ////////////////////////////////////////////////////// member functions

  void Set_default_cv_cov(TString default_cv_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir);
  
  ////////////////////////////////////////////////////// data members

 private:
  int default_oldworld_rows;
  int default_newworld_rows;
  
};

#endif
