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

#include "./src/draw.icc"

////////////////////////////////////////////////////////////////////// main

void draw_check()
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

  TFile *file_noosc = new TFile("roofile_check_pred_noosc.root", "read");
  TMatrixD *matrix_pred_noosc = (TMatrixD*)file_noosc->Get("matrix_pred");
  
  TFile *file_n4 = new TFile("roofile_check_pred_dm2_7d3_s22th_0d36.root", "read");
  TMatrixD *matrix_pred_n4 = (TMatrixD*)file_n4->Get("matrix_pred");
  
  TH1D *h1_pred_noosc = new TH1D("h1_pred_noosc", "", 100, 0, 100);
  TH1D *h1_pred_n4 = new TH1D("h1_pred_n4", "", 100, 0, 100);
  TH1D *h1_ratio_wi2no = new TH1D("h1_ratio_wi2no", "", 100, 0, 100);

  for(int ibin=1; ibin<=25; ibin++) {
    double BNB_FC  = (*matrix_pred_noosc)(0,      ibin-1);
    double BNB_PC  = (*matrix_pred_noosc)(0, 26  +ibin-1);
    double NuMI_FC = (*matrix_pred_noosc)(0, 26*7+ibin-1);
    double NuMI_PC = (*matrix_pred_noosc)(0, 26*8+ibin-1);
    h1_pred_noosc->SetBinContent(     ibin, BNB_FC);
    h1_pred_noosc->SetBinContent(25  +ibin, BNB_PC);
    h1_pred_noosc->SetBinContent(25*2+ibin, NuMI_FC);
    h1_pred_noosc->SetBinContent(25*3+ibin, NuMI_PC);    
  }

  for(int ibin=1; ibin<=25; ibin++) {
    double BNB_FC  = (*matrix_pred_n4)(0,      ibin-1);
    double BNB_PC  = (*matrix_pred_n4)(0, 26  +ibin-1);
    double NuMI_FC = (*matrix_pred_n4)(0, 26*7+ibin-1);
    double NuMI_PC = (*matrix_pred_n4)(0, 26*8+ibin-1);
    h1_pred_n4->SetBinContent(     ibin, BNB_FC);
    h1_pred_n4->SetBinContent(25  +ibin, BNB_PC);
    h1_pred_n4->SetBinContent(25*2+ibin, NuMI_FC);
    h1_pred_n4->SetBinContent(25*3+ibin, NuMI_PC);    
  }

  for(int ibin=1; ibin<=100; ibin++) {
    double val_noosc = h1_pred_noosc->GetBinContent(ibin);
    double val_n4 = h1_pred_n4->GetBinContent(ibin);
    double ratio = 0;
    if( val_noosc ) ratio = val_n4/val_noosc;
    h1_ratio_wi2no->SetBinContent(ibin, ratio);
  }
  
  
  TLine *line_sub[3];
  for(int idx=0; idx<3; idx++) {
    double xx = (idx+1)*25;
    line_sub[idx] = new TLine( xx, 0, xx, 1e3);
    line_sub[idx]->SetLineStyle(7);
    line_sub[idx]->SetLineWidth(1);
  }
  
  roostr = "canv_cv_oscillation";
  TCanvas *canv_cv_oscillation = new TCanvas(roostr, roostr, 1000, 900);
  canv_cv_oscillation->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  //pad_top_no->SetLogy();
  pad_top_no->Draw(); pad_top_no->cd();
  h1_pred_noosc->Draw("hist");
  h1_pred_noosc->SetLineColor(kBlack);
  func_title_size(h1_pred_noosc, 0.055, 0.055, 0.055, 0.055);
  h1_pred_noosc->GetYaxis()->CenterTitle();
  h1_pred_noosc->SetYTitle("Events");
  h1_pred_n4->Draw("same hist");
  h1_pred_n4->SetLineColor(kRed);
  h1_pred_n4->SetLineWidth(2);;
  canv_cv_oscillation->cd(); pad_top_no->cd(); pad_top_no->Update(); double y2_pad_top_no = gPad->GetUymax();
  for(int idx=0; idx<3; idx++) {
    line_sub[idx]->SetY2( y2_pad_top_no );
    line_sub[idx]->Draw("same");
  }
  h1_pred_noosc->Draw("same axis");
  
  canv_cv_oscillation->cd();
  TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_no->Draw(); pad_bot_no->cd();
  h1_ratio_wi2no->SetMinimum(0.6);
  h1_ratio_wi2no->SetMaximum(1);
  h1_ratio_wi2no->Draw("hist");
  h1_ratio_wi2no->SetLineColor(kRed);
  func_title_size(h1_ratio_wi2no, 0.07, 0.07, 0.07, 0.07);
  func_xy_title(h1_ratio_wi2no, "Bin index", "Osc./No Osc.");
  h1_ratio_wi2no->GetXaxis()->CenterTitle();
  h1_ratio_wi2no->GetYaxis()->CenterTitle();
  h1_ratio_wi2no->GetYaxis()->SetTitleOffset(0.74);
  h1_ratio_wi2no->GetYaxis()->SetLabelOffset(0.008);
  h1_ratio_wi2no->GetYaxis()->SetNdivisions(506);
  h1_ratio_wi2no->Draw("same axis");

  canv_cv_oscillation->SaveAs("canv_cv_oscillation.png");
}
