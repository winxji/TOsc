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

////////////////////////////////////////////////////////////////////// ccc

void func_get_contours(TH2D *h2_CL, TGraph *gh_CL[3], int index)
{
  gROOT->SetBatch( 1 );
  
  TString roostr = "";
  
  const int Ncontour = 3;
  double contours[Ncontour] = {0};
  contours[0] = 0.9;
  contours[1] = 0.95;
  contours[2] = 0.9973;
  
  ///////
  roostr = TString::Format("canv_h2_CL_%d", index);
  TCanvas *canv_h2_CL = new TCanvas(roostr, roostr, 800, 600);
  h2_CL->SetStats(0);
  h2_CL->SetContour(Ncontour, contours);
  h2_CL->Draw("cont z list");
  canv_h2_CL->Update(); // Needed to force the plotting and retrieve the contours in TGraphs
     
  // Get Contours
  TObjArray *conts_mc = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TList* contLevel_mc = NULL;
  TGraph* gh_curv_mc[10] = {0};

  Int_t nGraphs_mc    = 0;
  Int_t TotalConts_mc = 0;

  if (conts_mc == NULL){
    printf("*** No Contours Were Extracted!\n");
    TotalConts_mc = 0;
    exit(1);
  } else {
    TotalConts_mc = conts_mc->GetSize();
  }

  printf("TotalConts_mc = %d\n", TotalConts_mc);

  for(int i = 0; i < TotalConts_mc; i++){
    contLevel_mc = (TList*)conts_mc->At(i);
    printf("Contour %d has %d Graphs\n", i, contLevel_mc->GetSize());
    nGraphs_mc += contLevel_mc->GetSize();
  }

  nGraphs_mc = 0;
  for(int i = 0; i < TotalConts_mc; i++){
    contLevel_mc = (TList*)conts_mc->At(i);

    // Get first graph from list on curves on this level
    gh_curv_mc[i] = (TGraph*)contLevel_mc->First();     
  }

  ///////  
  for(int ic=0; ic<Ncontour; ic++) {
    int np_curv = gh_curv_mc[ic]->GetN();
    for(int idx=0; idx<np_curv; idx++) {
      double t14 = 0;
      double m41 = 0;
      gh_curv_mc[ic]->GetPoint(idx, t14, m41);
      gh_CL[ic]->SetPoint(idx, pow(10., t14), pow(10., m41) );      
    }    
  }

  delete canv_h2_CL;
}

//////////////////////////////////////////////////////////////////
///////////////////////// MAIN ///////////////////////////////////
//////////////////////////////////////////////////////////////////

void plot_CLs()
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
       
  ///////
  const int Ncontour = 3;
  //int colors[Ncontour] = {kBlue, kGreen+3, kRed};
  int colors[Ncontour] = {kBlue, kBlack, kRed};
  
  ////////////////////////////////////////////////////////////////////////////////////////

  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  TString file_roostr = "zz_NuMI_20by20_all.dat";
  
  int bins_theta = 20;
  int bins_dm2   = 20;
  
  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103

  roostr = "h2_space_data";
  TH2D *h2_space_data = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_pred";
  TH2D *h2_space_pred = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_pred_1sigma_plus";
  TH2D *h2_space_pred_1sigma_plus = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_pred_1sigma_minus";
  TH2D *h2_space_pred_1sigma_minus = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_pred_2sigma_plus";
  TH2D *h2_space_pred_2sigma_plus = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_pred_2sigma_minus";
  TH2D *h2_space_pred_2sigma_minus = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  ifstream InputFile(file_roostr, ios::in);
  if(!InputFile) { cerr<<" No input-list"<<endl; exit(1); }

  for(int idx=1; idx<=bins_theta*bins_dm2; idx++) {
    int theta(0), dm2(0);    
    double chi2_4v_4vAsimov(0), chi2_3v_4vAsimov(0),
      chi2_4v_3vAsimov(0), chi2_3v_3vAsimov(0),
      chi2_4v_data(0), chi2_3v_data(0),
      CL_data(0), CL_pred(0),
      CL_pred_1sigma_plus(0), CL_pred_1sigma_minus(0),
      CL_pred_2sigma_plus(0), CL_pred_2sigma_minus(0);

    InputFile >> theta >> dm2
	      >> chi2_4v_4vAsimov >> chi2_3v_4vAsimov
	      >> chi2_4v_3vAsimov >> chi2_3v_3vAsimov
	      >> chi2_4v_data >> chi2_3v_data
	      >> CL_data >> CL_pred
	      >> CL_pred_1sigma_plus >> CL_pred_1sigma_minus
	      >> CL_pred_2sigma_plus >> CL_pred_2sigma_minus;

    CL_data *= 0.01;
    CL_pred *= 0.01;
    CL_pred_1sigma_plus  *= 0.01;
    CL_pred_1sigma_minus *= 0.01;
    CL_pred_2sigma_plus  *= 0.01;
    CL_pred_2sigma_minus *= 0.01;
    
    h2_space_data->SetBinContent( theta, dm2, CL_data );
    h2_space_pred->SetBinContent( theta, dm2, CL_pred );
    h2_space_pred_1sigma_plus->SetBinContent( theta, dm2, CL_pred_1sigma_plus );
    h2_space_pred_1sigma_minus->SetBinContent( theta, dm2, CL_pred_1sigma_minus );
    h2_space_pred_2sigma_plus->SetBinContent( theta, dm2, CL_pred_2sigma_plus );
    h2_space_pred_2sigma_minus->SetBinContent( theta, dm2, CL_pred_2sigma_minus );
  }

  ////////////////////////////////////////////////////////////

  int index = 0;
  
  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_data[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_data_%d", idx);
    gh_CLs_data[idx] = new TGraph();
    gh_CLs_data[idx]->SetName(roostr);
    gh_CLs_data[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_data, gh_CLs_data, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_%d", idx);
    gh_CLs_pred[idx] = new TGraph();
    gh_CLs_pred[idx]->SetName(roostr);
    gh_CLs_pred[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred, gh_CLs_pred, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred_1sigma_plus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_1sigma_plus_%d", idx);
    gh_CLs_pred_1sigma_plus[idx] = new TGraph();
    gh_CLs_pred_1sigma_plus[idx]->SetName(roostr);
    gh_CLs_pred_1sigma_plus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_1sigma_plus, gh_CLs_pred_1sigma_plus, index);

  index++;
  TGraph *gh_CLs_pred_1sigma_minus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_1sigma_minus_%d", idx);
    gh_CLs_pred_1sigma_minus[idx] = new TGraph();
    gh_CLs_pred_1sigma_minus[idx]->SetName(roostr);
    gh_CLs_pred_1sigma_minus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_1sigma_minus, gh_CLs_pred_1sigma_minus, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred_2sigma_plus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_2sigma_plus_%d", idx);
    gh_CLs_pred_2sigma_plus[idx] = new TGraph();
    gh_CLs_pred_2sigma_plus[idx]->SetName(roostr);
    gh_CLs_pred_2sigma_plus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_2sigma_plus, gh_CLs_pred_2sigma_plus, index);

  index++;
  TGraph *gh_CLs_pred_2sigma_minus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_2sigma_minus_%d", idx);
    gh_CLs_pred_2sigma_minus[idx] = new TGraph();
    gh_CLs_pred_2sigma_minus[idx]->SetName(roostr);
    gh_CLs_pred_2sigma_minus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_2sigma_minus, gh_CLs_pred_2sigma_minus, index);

  ////////////////////////////////////////////////////////////

  TGraph *gh_contour_1sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_contour_1sigma_%d", idx);
    gh_contour_1sigma[idx] = new TGraph();
    gh_contour_1sigma[idx]->SetName(roostr);
    gh_contour_1sigma[idx]->SetFillColor(kGreen);
    
    int num_plus = gh_CLs_pred_1sigma_plus[idx]->GetN();
    for(int ip=0; ip<num_plus; ip++) {
      double xx(0), yy(0);
      gh_CLs_pred_1sigma_plus[idx]->GetPoint(ip, xx, yy);
      gh_contour_1sigma[idx]->SetPoint( gh_contour_1sigma[idx]->GetN(), xx, yy );
    }// ip
    
    int num_minus = gh_CLs_pred_1sigma_minus[idx]->GetN();
    for(int ip=num_minus-1; ip>=0; ip--) {
      double xx(0), yy(0);
      gh_CLs_pred_1sigma_minus[idx]->GetPoint(ip, xx, yy);
      gh_contour_1sigma[idx]->SetPoint( gh_contour_1sigma[idx]->GetN(), xx, yy );
    }// ip          
  }// idx
  
  TGraph *gh_contour_2sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_contour_2sigma_%d", idx);
    gh_contour_2sigma[idx] = new TGraph();
    gh_contour_2sigma[idx]->SetName(roostr);
    gh_contour_2sigma[idx]->SetFillColor(kYellow);
    
    int num_plus = gh_CLs_pred_2sigma_plus[idx]->GetN();
    for(int ip=0; ip<num_plus; ip++) {
      double xx(0), yy(0);
      gh_CLs_pred_2sigma_plus[idx]->GetPoint(ip, xx, yy);
      gh_contour_2sigma[idx]->SetPoint( gh_contour_2sigma[idx]->GetN(), xx, yy );
    }// ip
    
    int num_minus = gh_CLs_pred_2sigma_minus[idx]->GetN();
    for(int ip=num_minus-1; ip>=0; ip--) {
      double xx(0), yy(0);
      gh_CLs_pred_2sigma_minus[idx]->GetPoint(ip, xx, yy);
      gh_contour_2sigma[idx]->SetPoint( gh_contour_2sigma[idx]->GetN(), xx, yy );
    }// ip          
  }// idx
  
  //////////////////////////////////////////////////////////// plotting  
  //////////////////////////////////////////////////////////// plotting
  
  gROOT->SetBatch( 0 );  
  
  int index_90 = 0;
  int index_95 = 1;
  int index_99 = 2;
  
  double xxlow = h2_space_data->GetXaxis()->GetBinLowEdge(1);
  double xxhgh = h2_space_data->GetXaxis()->GetBinUpEdge(bins_theta);
    
  double yylow = h2_space_data->GetYaxis()->GetBinLowEdge(1);
  double yyhgh = h2_space_data->GetYaxis()->GetBinUpEdge(bins_dm2);

  
  TGraph *gh_n4 = new TGraph();
  gh_n4->SetPoint(0, 0.26, 7.2);
  gh_n4->SetMarkerStyle(22);
  gh_n4->SetMarkerSize(1.8);
  gh_n4->SetMarkerColor(kRed);
 
  TGraph *gh_n4new = new TGraph();
  gh_n4new->SetPoint(0, 0.36, 7.3);
  gh_n4new->SetMarkerStyle(22);
  gh_n4new->SetMarkerSize(1.8);
  gh_n4new->SetMarkerColor(kBlue);

  
  roostr = "h2_basic_CLs_data";
  TH2D *h2_basic_CLs_data = new TH2D(roostr, "",
				     bins_theta, pow(10, xxlow), pow(10, xxhgh),
				     bins_dm2, pow(10, yylow), pow(10, yyhgh));


  /////////////////////////////////////////////////////// 95
  
  roostr = "canv_h2_basic_CLs_data_95";
  TCanvas *canv_h2_basic_CLs_data_95 = new TCanvas(roostr, roostr, 800, 600);
  canv_h2_basic_CLs_data_95->SetLeftMargin(0.15);
  canv_h2_basic_CLs_data_95->SetRightMargin(0.1);
  canv_h2_basic_CLs_data_95->SetBottomMargin(0.18);
  canv_h2_basic_CLs_data_95->SetLogx();
  canv_h2_basic_CLs_data_95->SetLogy();
  
  h2_basic_CLs_data->Draw();
  h2_basic_CLs_data->SetXTitle("sin^{2}2#theta_{14}");
  h2_basic_CLs_data->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  h2_basic_CLs_data->GetXaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleSize(0.045);    
  h2_basic_CLs_data->GetXaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetXaxis()->SetTitleSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.4);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.4);

  //////////////
  
  gh_contour_2sigma[index_95]->Draw("same f");
  gh_contour_1sigma[index_95]->Draw("same f");
  
  gh_CLs_data[index_95]->Draw("same l");
  gh_CLs_data[index_95]->SetLineColor(kBlack);
    
  gh_CLs_pred[index_95]->Draw("same l");
  gh_CLs_pred[index_95]->SetLineStyle(7);
  
  gh_n4->Draw("same p");
  gh_n4new->Draw("same p");
 
  /////////////////////////////////////////////////////// 99
  
  roostr = "canv_h2_basic_CLs_data_99";
  TCanvas *canv_h2_basic_CLs_data_99 = new TCanvas(roostr, roostr, 800, 600);
  canv_h2_basic_CLs_data_99->SetLeftMargin(0.15);
  canv_h2_basic_CLs_data_99->SetRightMargin(0.1);
  canv_h2_basic_CLs_data_99->SetBottomMargin(0.18);
  canv_h2_basic_CLs_data_99->SetLogx();
  canv_h2_basic_CLs_data_99->SetLogy();
  
  h2_basic_CLs_data->Draw();
  h2_basic_CLs_data->SetXTitle("sin^{2}2#theta_{14}");
  h2_basic_CLs_data->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  h2_basic_CLs_data->GetXaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleSize(0.045);    
  h2_basic_CLs_data->GetXaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetXaxis()->SetTitleSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.4);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.4);

  //////////////
  
  gh_contour_2sigma[index_99]->Draw("same f");
  gh_contour_1sigma[index_99]->Draw("same f");
  
  gh_CLs_data[index_99]->Draw("same l");
  gh_CLs_data[index_99]->SetLineColor(kBlack);
    
  gh_CLs_pred[index_99]->Draw("same l");
  gh_CLs_pred[index_99]->SetLineStyle(7);
  
  gh_n4->Draw("same p");
  gh_n4new->Draw("same p");
 
}
