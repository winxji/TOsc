#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

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

#include "TImage.h"

#include "./src/draw.icc"

////////////////////////////////////////////////////////////////////// main

void draw_CLs_edit()
{

  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";
       
  ////////////////////////////////////////////////////////////////////////////////////////

  TFile *roofile_BNBNuMI_both = new TFile("./zm_nueapp/za_roofile_CL_both_right_fixt14_0d99.root", "read");
  TH2D *h2_basic_pars = (TH2D*)roofile_BNBNuMI_both->Get("h2_basic_CLs_data");
  TGraph *gh_CLs_data_1c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_CLs_data_1");
  TGraph *gh_CLs_pred_1c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_CLs_pred_1");
  TGraph *gh_CLs_pred_1sigma_1c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_contour_1sigma_1");
  TGraph *gh_CLs_pred_2sigma_1c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_contour_2sigma_1");
  TGraph *gh_wilk_CL_pred_1c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_wilk_CL_pred_1");
  TGraph *gh_CLs_data_2c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_CLs_data_2");
  TGraph *gh_CLs_pred_2c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_CLs_pred_2");
  TGraph *gh_CLs_pred_1sigma_2c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_contour_1sigma_2");
  TGraph *gh_CLs_pred_2sigma_2c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_contour_2sigma_2");
  TGraph *gh_wilk_CL_pred_2c_both = (TGraph*)roofile_BNBNuMI_both->Get("gh_wilk_CL_pred_2");
  gh_CLs_pred_1c_both->SetLineStyle(1);
  gh_CLs_pred_2c_both->SetLineStyle(1);
  gh_wilk_CL_pred_1c_both->SetLineStyle(1);
  gh_wilk_CL_pred_2c_both->SetLineStyle(1);  
  
  TFile *roofile_BNBNuMI_BNBonly = new TFile("./zm_nueapp/za_roofile_CL_both_right_fixt14_0d90.root", "read");
  //TH2D *h2_basic_pars = (TH2D*)roofile_BNBNuMI_BNBonly->Get("h2_basic_CLs_data");
  TGraph *gh_CLs_data_1c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_CLs_data_1");
  TGraph *gh_CLs_pred_1c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_CLs_pred_1");
  TGraph *gh_CLs_pred_1sigma_1c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_contour_1sigma_1");
  TGraph *gh_CLs_pred_2sigma_1c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_contour_2sigma_1");
  TGraph *gh_wilk_CL_pred_1c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_wilk_CL_pred_1");
  TGraph *gh_CLs_data_2c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_CLs_data_2");
  TGraph *gh_CLs_pred_2c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_CLs_pred_2");
  TGraph *gh_CLs_pred_1sigma_2c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_contour_1sigma_2");
  TGraph *gh_CLs_pred_2sigma_2c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_contour_2sigma_2");
  TGraph *gh_wilk_CL_pred_2c_BNBonly = (TGraph*)roofile_BNBNuMI_BNBonly->Get("gh_wilk_CL_pred_2");
  gh_CLs_pred_1c_BNBonly->SetLineStyle(1);
  gh_CLs_pred_2c_BNBonly->SetLineStyle(1);
  gh_wilk_CL_pred_1c_BNBonly->SetLineStyle(1);
  gh_wilk_CL_pred_2c_BNBonly->SetLineStyle(1);  
  
  TFile *roofile_BNBNuMI_NuMIonly = new TFile("./zm_nueapp/za_roofile_CL_both_right_fixt14_0d50.root", "read");
  //TH2D *h2_basic_pars = (TH2D*)roofile_BNBNuMI_NuMIonly->Get("h2_basic_CLs_data");
  TGraph *gh_CLs_data_1c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_CLs_data_1");
  TGraph *gh_CLs_pred_1c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_CLs_pred_1");
  TGraph *gh_CLs_pred_1sigma_1c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_contour_1sigma_1");
  TGraph *gh_CLs_pred_2sigma_1c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_contour_2sigma_1");
  TGraph *gh_wilk_CL_pred_1c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_wilk_CL_pred_1");
  TGraph *gh_CLs_data_2c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_CLs_data_2");
  TGraph *gh_CLs_pred_2c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_CLs_pred_2");
  TGraph *gh_CLs_pred_1sigma_2c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_contour_1sigma_2");
  TGraph *gh_CLs_pred_2sigma_2c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_contour_2sigma_2");
  TGraph *gh_wilk_CL_pred_2c_NuMIonly = (TGraph*)roofile_BNBNuMI_NuMIonly->Get("gh_wilk_CL_pred_2");
  gh_CLs_pred_1c_NuMIonly->SetLineStyle(1);
  gh_CLs_pred_2c_NuMIonly->SetLineStyle(1);
  gh_wilk_CL_pred_1c_NuMIonly->SetLineStyle(1);
  gh_wilk_CL_pred_2c_NuMIonly->SetLineStyle(1);  
       
  ////////////////////////////////////////////////////////////////////////////////////////
  
  TGraph *gh_n4new = new TGraph();
  gh_n4new->SetPoint(0, 0.36, 7.3);
  gh_n4new->SetMarkerStyle(20);
  gh_n4new->SetMarkerSize(1.8);
  gh_n4new->SetMarkerColor(kRed);

  roostr = "canv_h2_basic_CLs_AA";
  TCanvas *canv_h2_basic_CLs_AA = new TCanvas(roostr, roostr, 800, 600);  
  func_canv_margin(canv_h2_basic_CLs_AA, 0.15, 0.1, 0.1, 0.18);
  canv_h2_basic_CLs_AA->SetLogx(); canv_h2_basic_CLs_AA->SetLogy();
  h2_basic_pars->Draw();
  h2_basic_pars->SetXTitle("sin^{2}2#theta_{#mue}");

  if( 0 ) {
    gh_CLs_pred_2sigma_1c_NuMIonly->Draw("same f");
    gh_CLs_pred_1sigma_1c_NuMIonly->Draw("same f");
    gh_CLs_pred_1c_NuMIonly->Draw("same l");
    //gh_CLs_data_1c_NuMIonly->Draw("same l");
    gh_wilk_CL_pred_1c_NuMIonly->Draw("same l");
  }

  if( 0 ) {
    //gh_CLs_pred_2sigma_1c_BNBonly->Draw("same f");
    //gh_CLs_pred_1sigma_1c_BNBonly->Draw("same f");
    gh_CLs_pred_1c_BNBonly->Draw("same l");
    //gh_CLs_data_1c_BNBonly->Draw("same l");
    gh_wilk_CL_pred_1c_BNBonly->Draw("same l"); gh_wilk_CL_pred_1c_BNBonly->SetLineColor(kRed);

    h2_basic_pars->GetXaxis()->SetLabelSize(0);
    h2_basic_pars->GetXaxis()->SetTitleSize(0);
    h2_basic_pars->GetYaxis()->SetLabelSize(0);
    h2_basic_pars->GetYaxis()->SetTitleSize(0);
    
  }

  if( 0 ) {
    gh_CLs_pred_2sigma_1c_both->Draw("same f");
    gh_CLs_pred_1sigma_1c_both->Draw("same f");
    gh_CLs_pred_1c_both->Draw("same l");
    //gh_CLs_data_1c_both->Draw("same l");
    gh_wilk_CL_pred_1c_both->Draw("same l");
  }
  
  if( 0 ) {
    gh_CLs_pred_2sigma_2c_both->Draw("same f");
    gh_CLs_pred_1sigma_2c_both->Draw("same f");
    gh_CLs_pred_2c_both->Draw("same l");
    //gh_CLs_data_2c_both->Draw("same l");    
  }

  if( 0 ) {
    gh_CLs_pred_1c_BNBonly->Draw("same l");
    gh_CLs_pred_1c_NuMIonly->Draw("same l");
    gh_CLs_pred_1c_both->Draw("same l");
  }

  if( 0 ) {
    gh_CLs_pred_1c_both->Draw("same l");          gh_CLs_pred_1c_both->SetLineColor(kBlue);
    //gh_CLs_pred_1c_both_run123->Draw("same l");   gh_CLs_pred_1c_both_run123->SetLineColor(kRed);
    //gh_CLs_pred_1c_both_runwhole->Draw("same l"); gh_CLs_pred_1c_both_runwhole->SetLineColor(kGreen+1);

    //gh_wilk_CL_pred_1c_both->Draw("same l");        gh_wilk_CL_pred_1c_both->SetLineColor(kBlue);       gh_wilk_CL_pred_1c_both->SetLineStyle(7);
    //gh_wilk_CL_pred_1c_both_run123->Draw("same l"); gh_wilk_CL_pred_1c_both_run123->SetLineColor(kRed); gh_wilk_CL_pred_1c_both_run123->SetLineStyle(7);    
  }


  if( 1 ) {// za_roofile_CL_both_t24_fix
    //gh_CLs_pred_2c_BNBonly->Draw("same l");  gh_CLs_pred_2c_BNBonly->SetLineColor(kBlue);  // fix 0.01
    gh_CLs_pred_2c_both->Draw("same l");     gh_CLs_pred_2c_both->SetLineColor(kRed);      // fix 0    
    //gh_CLs_pred_2c_NuMIonly->Draw("same l"); gh_CLs_pred_2c_NuMIonly->SetLineColor(kGreen);// fix 0.1
  }
  
  
  //////////////
  
  //gh_n4new->Draw("same p");

  //TLegend *lg_AA = new TLegend(0.18, 0.2, 0.5, 0.5);
  TLegend *lg_AA = new TLegend(0.5, 0.7, 0.75, 0.85);
  
  lg_AA->SetBorderSize(0); lg_AA->SetFillStyle(0); lg_AA->SetTextSize(0.05);
  
  //lg_AA->AddEntry(gh_CLs_pred_1c_BNBonly, "CLs Sensitivity, 95% CL", "l");
  //lg_AA->AddEntry(gh_CLs_pred_1sigma_1c_BNBonly, "1sigma", "f");
  //lg_AA->AddEntry(gh_CLs_pred_2sigma_1c_BNBonly, "2sigma", "f");
  //lg_AA->AddEntry(gh_wilk_CL_pred_1c_BNBonly, "Wilks Sensitivity, 95% CL", "l");
  //lg_AA->AddEntry(gh_CLs_data_1c_BNBonly, "CLs Exclusion, 95% CL", "l");
  
  //lg_AA->AddEntry("", "CLs Sensitivity, 95% CL", "");
  //lg_AA->AddEntry(gh_CLs_pred_1c_BNBonly, "BNB only", "l");     gh_CLs_pred_1c_BNBonly->SetLineColor(kGreen+1);
  //lg_AA->AddEntry(gh_CLs_pred_1c_NuMIonly, "NuMI only", "l");   gh_CLs_pred_1c_NuMIonly->SetLineColor(kBlue);
  //lg_AA->AddEntry(gh_CLs_pred_1c_both, "Both", "l");            gh_CLs_pred_1c_both->SetLineColor(kRed);

  lg_AA->AddEntry("", "CLs Sensitivity, 99.75% CL", "");  
  lg_AA->AddEntry(gh_CLs_pred_2c_both, "sin^{2}2#theta_{14} = 0.04", "l");
  //lg_AA->AddEntry(gh_CLs_pred_2c_BNBonly, "sin^{2}2#theta_{14} = 0.36", "l");
  //lg_AA->AddEntry(gh_CLs_pred_2c_NuMIonly, "sin^{2}#theta_{24} = 0.1", "l");
  
  
  lg_AA->Draw();
  
  canv_h2_basic_CLs_AA->SaveAs("canv_fix_t14.png");


  
  /*
  TFile *outfile = new TFile("outfile_gg.root", "recreate");
  gh_CLs_pred_1c_both->SetName("gh_CLs_pred_1c_both"); gh_CLs_pred_1c_both->Write();
  gh_CLs_pred_1c_both_run123->SetName("gh_CLs_pred_1c_both_run123"); gh_CLs_pred_1c_both_run123->Write();
  outfile->Close();
  */

  
}
