#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Plot_user()
{
  if( 0 ) {
    gStyle->SetPalette(kTemperatureMap);
    gStyle->SetNumberContours(100);
      
    TString roostr = "";
    roostr = "h2_corr_newworld";
    TH2D *h2_corr_newworld = new TH2D(roostr, "", default_newworld_rows, 0, default_newworld_rows, default_newworld_rows, 0, default_newworld_rows);
    
    for(int idx=0; idx<default_newworld_rows; idx++) {
      for(int jdx=0; jdx<default_newworld_rows; jdx++) {
	double cov = matrix_eff_newworld_abs_syst_total(idx, jdx);
	double si = sqrt( matrix_eff_newworld_abs_syst_total(idx, idx) );
	double sj = sqrt( matrix_eff_newworld_abs_syst_total(jdx, jdx) );
	double corr = 0;
	if( si!=0 && sj!=0 ) corr = cov/si/sj;
	if( idx==jdx ) corr = 1;
	h2_corr_newworld->SetBinContent(idx, jdx, corr);
      }
    }

    TLine *line_sub_XX[13];
    for(int idx=0; idx<13; idx++) {
      double xx = 0;
      for(int jdx=0; jdx<=idx; jdx++) xx += (map_default_h1d_meas_bins[idx+1]+1);// hack
      line_sub_XX[idx] = new TLine( xx, 0, xx, default_newworld_rows);
      line_sub_XX[idx]->SetLineStyle(7);
      line_sub_XX[idx]->SetLineWidth(1);
    }    
    
    TLine *line_sub_YY[13];
    for(int idx=0; idx<13; idx++) {
      double xx = 0;
      for(int jdx=0; jdx<=idx; jdx++) xx += (map_default_h1d_meas_bins[idx+1]+1);// hack
      line_sub_YY[idx] = new TLine( 0, xx, default_newworld_rows, xx);
      line_sub_YY[idx]->SetLineStyle(7);
      line_sub_YY[idx]->SetLineWidth(1);
    }    
    
    TCanvas *canv_h2_corr_newworld = new TCanvas(roostr, roostr, 800, 700);
    func_canv_margin(canv_h2_corr_newworld, 0.15, 0.15, 0.1, 0.15);
    h2_corr_newworld->Draw("colz");
    h2_corr_newworld->GetZaxis()->SetRangeUser(-1,1);   
    func_title_size(h2_corr_newworld, 0.045, 0.045, 0.045, 0.045);
    h2_corr_newworld->GetZaxis()->SetLabelSize(0.04);
    h2_corr_newworld->GetZaxis()->SetRangeUser(-1, 1);
    h2_corr_newworld->GetXaxis()->SetTickLength(0);
    h2_corr_newworld->GetYaxis()->SetTickLength(0);
    func_xy_title(h2_corr_newworld, "Bin index", "Bin index");
    h2_corr_newworld->GetXaxis()->CenterTitle(); h2_corr_newworld->GetYaxis()->CenterTitle();
    h2_corr_newworld->GetXaxis()->SetTitleOffset(1.1); h2_corr_newworld->GetYaxis()->SetTitleOffset(1.3);
    
    for(int idx=0; idx<13; idx++) {
      line_sub_XX[idx]->Draw("same");
      line_sub_YY[idx]->Draw("same");
    }

    canv_h2_corr_newworld->SaveAs("canv_h2_corr_newworld.png");
  }

  if( 1 ) {
    TFile *roofile_check = new TFile("roofile_check.root", "recreate");
    matrix_eff_newworld_pred.Write("matrix_pred");
    roofile_check->Close();
  }
  
  if( 0 ) {
    TString roostr = "";

    TF1 *f1_one = new TF1("f1_one", "1", 0, 1e6);
    f1_one->SetLineColor(kBlack);
    f1_one->SetLineStyle(7);
    
    TF1 *f1_zero = new TF1("f1_zero", "0", 0, 1e6);
    f1_zero->SetLineColor(kBlack);
    f1_zero->SetLineStyle(7);
    
    roostr = "h1_newworld_pred_test";
    TH1D *h1_newworld_pred_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_pred_relerr_test";
    TH1D *h1_newworld_pred_relerr_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_meas_test";
    TH1D *h1_newworld_meas_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_meas2pred_test";
    TH1D *h1_newworld_meas2pred_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    for(int ibin=1; ibin<=default_newworld_rows; ibin++) {
      double content = matrix_eff_newworld_pred(0, ibin-1);
      double abs_err = sqrt( matrix_eff_newworld_abs_syst_total(ibin-1, ibin-1) );
      h1_newworld_pred_test->SetBinContent(ibin, content);
      h1_newworld_pred_test->SetBinError(ibin, abs_err);

      double relerr = 0;
      if(content!=0) relerr = abs_err/content;
      h1_newworld_pred_relerr_test->SetBinError(ibin, relerr);
      h1_newworld_pred_relerr_test->SetBinContent(ibin, 0);

      double meas = matrix_eff_newworld_meas(0, ibin-1);
      h1_newworld_meas_test->SetBinContent(ibin, meas );

      double ratio_val = 0;
      double ratio_err = 0;
      if(content!=0) {
	ratio_val = meas/content;
	ratio_err = sqrt(meas)/content;
      }
      h1_newworld_meas2pred_test->SetBinContent(ibin, ratio_val);
      h1_newworld_meas2pred_test->SetBinError(ibin, ratio_err);
    }

    roostr = "canv_h1_newworld_pred_test";
    TCanvas *canv_h1_newworld_pred_test = new TCanvas(roostr, roostr, 1000, 600);
    func_canv_margin(canv_h1_newworld_pred_test, 0.1, 0.1, 0.1, 0.15);
    canv_h1_newworld_pred_test->SetLogy();
    TH1D *h1_newworld_pred_test_err = (TH1D*)h1_newworld_pred_test->Clone("h1_newworld_pred_test_err");
    h1_newworld_pred_test_err->Draw("e2");
    h1_newworld_pred_test_err->SetMinimum(1e-1);
    h1_newworld_pred_test_err->SetFillColor(kRed-9);
    h1_newworld_pred_test_err->SetMarkerSize(0);
    h1_newworld_pred_test->Draw("same hist");
    h1_newworld_pred_test->SetLineColor(kBlack);
    h1_newworld_meas_test->Draw("same e1");
    h1_newworld_meas_test->SetLineColor(kBlue);
    h1_newworld_meas_test->SetMarkerColor(kBlue);
    canv_h1_newworld_pred_test->SaveAs("canv_h1_newworld_pred_test.root");
    canv_h1_newworld_pred_test->SaveAs("canv_h1_newworld_pred_test.png");
      
    roostr = "canv_h1_newworld_pred_relerr_test";
    TCanvas *canv_h1_newworld_pred_relerr_test = new TCanvas(roostr, roostr, 1000, 600);
    func_canv_margin(canv_h1_newworld_pred_relerr_test, 0.1, 0.1, 0.1, 0.15);
    h1_newworld_pred_relerr_test->SetMinimum(0);
    h1_newworld_pred_relerr_test->SetMaximum(2);
    //h1_newworld_pred_relerr_test->SetFillColor(kRed-9);
    //h1_newworld_pred_relerr_test->SetMarkerSize(0);
    h1_newworld_pred_relerr_test->Draw("e2");
    h1_newworld_pred_relerr_test->SetFillColor(kRed-9);
    h1_newworld_pred_relerr_test->SetLineColor(kRed);
    h1_newworld_pred_relerr_test->SetLineStyle(7);
    h1_newworld_pred_relerr_test->SetMarkerSize(0);
    h1_newworld_meas2pred_test->Draw("same e1");
    h1_newworld_meas2pred_test->SetLineColor(kBlue);
    h1_newworld_meas2pred_test->SetMarkerColor(kBlue);
    f1_one->Draw("same");
    canv_h1_newworld_pred_relerr_test->SaveAs("canv_h1_newworld_pred_relerr_test.root");
    canv_h1_newworld_pred_relerr_test->SaveAs("canv_h1_newworld_pred_relerr_test.png");

    /////////////

    TLine *line_sub[13];
    for(int idx=0; idx<13; idx++) {
      double xx = 0;
      for(int jdx=0; jdx<=idx; jdx++) xx += (map_default_h1d_meas_bins[idx+1]+1);// hack
      line_sub[idx] = new TLine( xx, 0, xx, 1e3);
      line_sub[idx]->SetLineStyle(7);
      line_sub[idx]->SetLineWidth(1);
    }
    
    roostr = "canv_cv_and_err";
    TCanvas *canv_cv_and_err = new TCanvas(roostr, roostr, 1000, 900);
    canv_cv_and_err->cd();
    TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
    func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
    pad_top_no->SetLogy();
    pad_top_no->Draw(); pad_top_no->cd();
    h1_newworld_pred_test_err->Draw("e2");
    func_title_size(h1_newworld_pred_test_err, 0.055, 0.055, 0.055, 0.055);
    h1_newworld_pred_test_err->GetYaxis()->CenterTitle();
    h1_newworld_pred_test_err->SetYTitle("Events");
    h1_newworld_pred_test->Draw("same hist");
    h1_newworld_pred_test->SetLineWidth(2);;
    canv_cv_and_err->cd(); pad_top_no->cd(); pad_top_no->Update(); double y2_pad_top_no = gPad->GetUymax();
    for(int idx=0; idx<13; idx++) {
      line_sub[idx]->SetY2( pow(10,y2_pad_top_no) );
      line_sub[idx]->Draw("same");
    }
	
    canv_cv_and_err->cd();
    TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
    func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
    pad_bot_no->Draw(); pad_bot_no->cd();
    h1_newworld_pred_relerr_test->SetMinimum(-1);
    h1_newworld_pred_relerr_test->SetMaximum(1);
    h1_newworld_pred_relerr_test->Draw("e2");
    func_title_size(h1_newworld_pred_relerr_test, 0.07, 0.07, 0.07, 0.07);
    func_xy_title(h1_newworld_pred_relerr_test, "Bin index", "Rel. Uncertainty");
    h1_newworld_pred_relerr_test->GetXaxis()->CenterTitle();
    h1_newworld_pred_relerr_test->GetYaxis()->CenterTitle();
    h1_newworld_pred_relerr_test->GetYaxis()->SetTitleOffset(0.84);
    h1_newworld_pred_relerr_test->GetYaxis()->SetLabelOffset(0.008);
    f1_zero->Draw("same");
    
    canv_cv_and_err->SaveAs("canv_cv_and_err.png");
    canv_cv_and_err->SaveAs("canv_cv_and_err.root");
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

double TOsc::FCN(const double *par)
{
  double chi2_final = 0;
  double fit_dm2_41         = par[0];
  double fit_sin2_2theta_14 = par[1];
        
  /////// standard order
  Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, 0, 0);  
  Apply_oscillation();
  Set_apply_POT();// meas, CV, COV: all ready

  ///////
  
  TMatrixD matrix_data_total = matrix_fitdata_newworld;
  TMatrixD matrix_pred_total = matrix_eff_newworld_pred;        
  TMatrixD matrix_cov_syst_total = matrix_eff_newworld_abs_syst_total;
  int rows = matrix_cov_syst_total.GetNrows();
  
  /////// modify
  /*
  TMatrixD matrix_data_total = matrix_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
  TMatrixD matrix_pred_total = matrix_eff_newworld_pred.GetSub(0,0, 0, 26*7-1);
  TMatrixD matrix_cov_syst_total = matrix_eff_newworld_abs_syst_total.GetSub(0, 26*7-1, 0, 26*7-1);
  int rows = matrix_cov_syst_total.GetNrows();
  */
  /*
  TMatrixD matrix_data_total = matrix_fitdata_newworld.GetSub(0,0, 26*7, 26*14-1);
  TMatrixD matrix_pred_total = matrix_eff_newworld_pred.GetSub(0,0, 26*7, 26*14-1);
  TMatrixD matrix_cov_syst_total = matrix_eff_newworld_abs_syst_total.GetSub(26*7, 26*14-1, 26*7, 26*14-1);
  int rows = matrix_cov_syst_total.GetNrows();
  */
  ///////  

  TMatrixD matrix_cov_stat_total(rows, rows);
  TMatrixD matrix_cov_total(rows, rows);

  for(int idx=0; idx<rows; idx++) {
    double val_stat_cov = 0;        
    double val_data = matrix_data_total(0, idx);
    double val_pred = matrix_pred_total(0, idx);
        
    if( val_data==0 ) { val_stat_cov = val_pred/2; }
    else {
      if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
      else val_stat_cov = val_data;
    }

    if( val_stat_cov==0 ) val_stat_cov = 1e-6;
    if( matrix_cov_syst_total(idx, idx)==0 ) matrix_cov_syst_total(idx, idx) = 1e-6;
    
    matrix_cov_stat_total(idx, idx) = val_stat_cov;    
  }
  
  matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
  
  ///////
	
  TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
  TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;
  TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();
   
  TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
  chi2_final = matrix_chi2(0,0);           

  return chi2_final;
}

///////

void TOsc::Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar)
{
  ROOT::Minuit2::Minuit2Minimizer min_osc( ROOT::Minuit2::kMigrad );
  min_osc.SetPrintLevel(2);
  min_osc.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min_osc.SetMaxFunctionCalls(50000);
  min_osc.SetMaxIterations(50000);
  min_osc.SetTolerance(1e-5); // tolerance*2e-3 = edm precision
  min_osc.SetPrecision(1e-18); //precision in the target function
    
  /// set fitting parameters
  ROOT::Math::Functor Chi2Functor_osc(
				      [&](const double *par) {return FCN( par );},// FCN
				      2// number of fitting parameters
				      );
    
  min_osc.SetFunction(Chi2Functor_osc);
    
  min_osc.SetVariable( 0, "dm2_41", init_dm2_41, 1e-2);
  min_osc.SetVariable( 1, "sin2_2theta_14", init_sin2_2theta_14, 1e-2);

  min_osc.SetLowerLimitedVariable(0, "dm2_41", init_dm2_41, 1e-3, 0);  
  min_osc.SetLimitedVariable(1, "sin2_2theta_14", init_sin2_2theta_14, 1e-3, 0, 1);
    
  min_osc.Minimize();

  ///////
    
  minimization_status = min_osc.Status();
  minimization_chi2   = min_osc.MinValue();
    
  const double *par_val = min_osc.X();
  const double *par_err = min_osc.Errors();

  minimization_dm2_41_val = par_val[0];
  minimization_dm2_41_err = par_err[0];
    
  minimization_sin2_2theta_14_val = par_val[1];
  minimization_sin2_2theta_14_err = par_err[1];

  if(1) {
    cout<<endl;
    cout<<TString::Format(" ---> minimization, status %2d, chi2 %6.4f, dm2 %4.2f +/- %4.2f, s22t14 %5.3f +/- %5.3f",
			  minimization_status, minimization_chi2,
			  minimization_dm2_41_val, minimization_dm2_41_err,
			  minimization_sin2_2theta_14_val, minimization_sin2_2theta_14_err
			  )<<endl;
  }
   
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_toy_variations(int num_toys)
{
  //cout<<endl<<" ---> Set_toy_variations "<<endl;

  map_matrix_toy_pred.clear();
  NUM_TOYS = num_toys;

  ///////
  TMatrixDSym DSmatrix_cov(default_newworld_rows);
  for(int idx=0; idx<default_newworld_rows; idx++) {
    for(int jdx=0; jdx<default_newworld_rows; jdx++) {
      DSmatrix_cov(idx, jdx) = matrix_eff_newworld_abs_syst_total(idx, jdx);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

  ///////
  for(int itoy=1; itoy<=num_toys; itoy++) {
    TMatrixD matrix_gaus_rand(default_newworld_rows, 1);
    int eff_count = 0;

  RANDOM_AGAIN:
    eff_count++;
    for(int idx=0; idx<default_newworld_rows; idx++) {
      if( matrix_eigenvalue(idx)>=0 ) {
	matrix_gaus_rand(idx, 0) = rand->Gaus( 0, sqrt(matrix_eigenvalue(idx)) );// systematics
      }
      else {
	matrix_gaus_rand(idx, 0) = 0;
      }
    }// for(int idx=0; idx<default_newworld_rows; idx++)

    TMatrixD matrix_variation = (matrix_eigenvector*matrix_gaus_rand).T();

    bool flag_negative = 0;

    for(int idx=0; idx<default_newworld_rows; idx++) {
      double val_with_syst = matrix_variation(0, idx) + matrix_eff_newworld_pred(0, idx);// key point
      if(val_with_syst<0) {
	if( matrix_eff_newworld_pred(0, idx)<0.8 ) {// hack for very low statistics
	  matrix_variation(0, idx) *= (-1);
	}
	else {
	  flag_negative = 1;
	  break;
	}
      }
    }// for(int idx=0; idx<default_newworld_rows; idx++)

    if( flag_negative ) goto RANDOM_AGAIN;

    TMatrixD matrix_temp_pred(1, default_newworld_rows);
    for(int idx=0; idx<default_newworld_rows; idx++) {
      double val_with_syst = matrix_variation(0, idx) + matrix_eff_newworld_pred(0, idx);// key point
      val_with_syst = rand->PoissonD(val_with_syst);// statistics
      matrix_temp_pred(0, idx) = val_with_syst;
    }// for(int idx=0; idx<default_newworld_rows; idx++)
    map_matrix_toy_pred[itoy].Clear(); map_matrix_toy_pred[itoy].ResizeTo(1, default_newworld_rows);
    map_matrix_toy_pred[itoy] = matrix_temp_pred;
    
  }// for(int itoy=1; itoy<=num_toys; itoy++)

  ///////
  
  if( 0 ) {// validation
    TPrincipal principal_cov(default_newworld_rows, "ND");
    double *array_pseudo_expt = new double[default_newworld_rows];
    for(int itoy=1; itoy<=num_toys; itoy++) {
      if(itoy%1000==0) cout<<TString::Format(" ---> Processing %6d", itoy)<<endl;
      for(int idx=0; idx<default_newworld_rows; idx++) {
	array_pseudo_expt[idx] = map_matrix_toy_pred[itoy](0, idx);
      }
      principal_cov.AddRow( array_pseudo_expt );
    }

    TMatrixD *matrix_principal_cov = (TMatrixD*)principal_cov.GetCovarianceMatrix();
    for(int idx=0; idx<default_newworld_rows; idx++) {
      for(int jdx=0; jdx<idx; jdx++) {
	(*matrix_principal_cov)(jdx,idx) = (*matrix_principal_cov)(idx,jdx);
      }      
    }

    for(int idx=0; idx<default_newworld_rows; idx++) {
      cout<<TString::Format("%3d cv %f, old,new: %f %f", idx+1, matrix_eff_newworld_pred(0, idx),
			    sqrt(matrix_eff_newworld_abs_syst_total(idx, idx)),
			    sqrt((*matrix_principal_cov)(idx,idx)))<<endl;
    }
  }// if( 0 ) {// validation
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_apply_POT()
{  
  ////////////////////////// hack for POT scale, should know which part is BNB, and which part is NuMI
  
  matrix_eff_newworld_meas = matrix_default_newworld_meas;
  matrix_eff_newworld_pred *= 0;
  matrix_eff_newworld_abs_syst_total *= 0;

  /////// hack, oldworld
  TMatrixD matrix_temp_oscillation_oldworld_pred = matrix_oscillation_oldworld_pred;
  TMatrixD matrix_temp_oldworld_abs_syst_addi = matrix_default_oldworld_abs_syst_addi;

  // for(int idx=0; idx<default_oldworld_rows; idx++) {
  //   matrix_temp_oscillation_oldworld_pred(0, idx) *= scaleF_POT_NuMI;
  //   matrix_temp_oldworld_abs_syst_addi(idx, idx) *= (scaleF_POT_NuMI*scaleF_POT_NuMI);
  // }

  for(int idx=0; idx<default_oldworld_rows; idx++) {
    if( idx<26*14 ) {
      matrix_temp_oscillation_oldworld_pred(0, idx) *= scaleF_POT_BNB;
      matrix_temp_oldworld_abs_syst_addi(idx, idx)  *= (scaleF_POT_BNB*scaleF_POT_BNB);
    }
    else {
      matrix_temp_oscillation_oldworld_pred(0, idx) *= scaleF_POT_NuMI;
      matrix_temp_oldworld_abs_syst_addi(idx, idx)  *= (scaleF_POT_NuMI*scaleF_POT_NuMI);
    }
  }

  /////// hack, newworld
  TMatrixD matrix_temp_newworld_abs_syst_mcstat = matrix_default_newworld_abs_syst_mcstat;
  for(int idx=0; idx<default_newworld_rows; idx++) {
    if( idx<26*7 ) {
      matrix_eff_newworld_meas(0, idx) *= scaleF_POT_BNB;
      matrix_temp_newworld_abs_syst_mcstat(idx, idx) *= (scaleF_POT_BNB*scaleF_POT_BNB);
    }
    else {
      matrix_eff_newworld_meas(0, idx) *= scaleF_POT_NuMI;
      matrix_temp_newworld_abs_syst_mcstat(idx, idx) *= (scaleF_POT_NuMI*scaleF_POT_NuMI);
    }
  }

  matrix_eff_newworld_pred = matrix_temp_oscillation_oldworld_pred * matrix_transform;
  
  ///////
  ///////
  
  TMatrixD matrix_temp_oldworld_abs_syst_flux(default_oldworld_rows, default_oldworld_rows);
  TMatrixD matrix_temp_oldworld_abs_syst_geant(default_oldworld_rows, default_oldworld_rows);
  TMatrixD matrix_temp_oldworld_abs_syst_Xs(default_oldworld_rows, default_oldworld_rows);
  TMatrixD matrix_temp_oldworld_abs_syst_det(default_oldworld_rows, default_oldworld_rows);

  for(int idx=0; idx<default_oldworld_rows; idx++) {        
    for(int jdx=0; jdx<default_oldworld_rows; jdx++) {
      double cv_i = matrix_temp_oscillation_oldworld_pred(0, idx);
      double cv_j = matrix_temp_oscillation_oldworld_pred(0, jdx);
      
      matrix_temp_oldworld_abs_syst_flux(idx, jdx)  = cv_i*cv_j*matrix_default_oldworld_rel_syst_flux(idx, jdx);
      matrix_temp_oldworld_abs_syst_geant(idx, jdx) = cv_i*cv_j*matrix_default_oldworld_rel_syst_geant(idx, jdx);
      matrix_temp_oldworld_abs_syst_Xs(idx, jdx)    = cv_i*cv_j*matrix_default_oldworld_rel_syst_Xs(idx, jdx);
      matrix_temp_oldworld_abs_syst_det(idx, jdx)   = cv_i*cv_j*matrix_default_oldworld_rel_syst_det(idx, jdx);      
    }
  }
  
  TMatrixD matrix_transform_T =  matrix_transform.T();  matrix_transform.T();

  matrix_eff_newworld_abs_syst_addi   = matrix_transform_T * matrix_temp_oldworld_abs_syst_addi * matrix_transform;
  matrix_eff_newworld_abs_syst_mcstat = matrix_temp_newworld_abs_syst_mcstat;
  matrix_eff_newworld_abs_syst_flux   = matrix_transform_T * matrix_temp_oldworld_abs_syst_flux * matrix_transform;
  matrix_eff_newworld_abs_syst_geant  = matrix_transform_T * matrix_temp_oldworld_abs_syst_geant * matrix_transform;
  matrix_eff_newworld_abs_syst_Xs     = matrix_transform_T * matrix_temp_oldworld_abs_syst_Xs * matrix_transform;
  matrix_eff_newworld_abs_syst_det    = matrix_transform_T * matrix_temp_oldworld_abs_syst_det * matrix_transform;

  if( flag_syst_dirt )   matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_addi;
  if( flag_syst_mcstat ) matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_mcstat;
  if( flag_syst_flux )   matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_flux;
  if( flag_syst_geant )  matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_geant;
  if( flag_syst_Xs )     matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_Xs;
  if( flag_syst_det )    matrix_eff_newworld_abs_syst_total += matrix_eff_newworld_abs_syst_det;  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

double TOsc::Prob_oscillaion(double Etrue, double baseline, TString strflag_osc)
{  
  double prob = 0;
  int flag_osc = -1;
    
  if( strflag_osc=="nue2nue" )        flag_osc = 1;
  else if( strflag_osc=="numu2numu" ) flag_osc = 2;
  else if( strflag_osc=="numu2nue" )  flag_osc = 3;
  else if( strflag_osc=="nue2numu" )  flag_osc = 4;
  else if( strflag_osc=="nueNC")      flag_osc = 5;
  else if( strflag_osc=="numuNC")     flag_osc = 6;
  else flag_osc = -1;

  if( flag_osc==-1 ) {
    cerr<<" ERROR: flag_osc==-1, strflag_osc=="<<strflag_osc<<" ---> Candidates: nue2nue numu2numu numu2nue nue2numu"<<endl;
    exit(1);
  }

  const int nue2nue   = 1;
  const int numu2numu = 2;
  const int numu2nue  = 3;
  const int nue2numu  = 4;
  const int nueNC     = 5;
  const int numuNC    = 6;

  /// dm2*L/4E = 1.267 * dm2(eV2) * L(m) / E(MeV)
    
  switch( flag_osc ) {
  case nue2nue:
    prob = 1 - sin2_2theta_14 * pow(TMath::Sin(1.267 * dm2_41 * baseline/Etrue), 2);
    break;
  case numu2numu:
    prob = 1 - sin2_2theta_14 * pow(TMath::Sin(1.267 * dm2_41 * baseline/Etrue), 2);// only numu2numu
    break;
  case numu2nue:
    break;
  case nue2numu:
    break;
  case nueNC:
    break;
  case numuNC:
    break;
  default:
    cerr<<"ERROR: NAN flag_osc"<<endl; exit(1);
  }
    
  return prob;
}

///////////////////

void TOsc::Set_oscillation_base_minus(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode)
{
  int total_pred_chs = map_default_h1d_pred.size();
  if( pred_channel_index > total_pred_chs ) { cerr<<TString::Format(" ERROR: pred_channel_index(%d) > total_pred_chs(%d)", pred_channel_index, total_pred_chs)<<endl; exit(1); }  
  cout<<TString::Format("            ---> (self-check) work on pred_channel: %3d,  oscillation mode: %-10s", pred_channel_index, str_osc_mode.Data())<<endl;

  
  for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ ) {
    TH1D *h1d_temp = (TH1D*)map_default_h1d_pred[pred_channel_index]->Clone("h1d_temp"); h1d_temp->Reset();// define and clear
    
    for(int ievent=0; ievent<(int)vec_vec_eventinfo->at(isize).size(); ievent++) {
      EventInfo info = vec_vec_eventinfo->at(isize).at(ievent);
      h1d_temp->Fill(info.e2e_Ereco, info.e2e_weight_xs);
    }

    h1d_temp->Scale( vec_ratioPOT->at(isize) );

    ///////    
    TMatrixD matrix_temp(1, default_oldworld_rows);        
    int bin_index_base = 0;
    if( pred_channel_index > 1 ) { for(int ich=1; ich<pred_channel_index; ich++) bin_index_base += ( map_default_h1d_pred[ich]->GetNbinsX()+1 ); }
    for(int ibin=1; ibin<=h1d_temp->GetNbinsX()+1; ibin++) { matrix_temp(0, bin_index_base + ibin -1) = h1d_temp->GetBinContent(ibin); }
    matrix_oscillation_base_oldworld_pred -= matrix_temp;    
  }// for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ )  
}

void TOsc::Set_oscillation_base_added(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode)
{
  int total_pred_chs = map_default_h1d_pred.size();
  if( pred_channel_index > total_pred_chs ) { cerr<<TString::Format(" ERROR: pred_channel_index(%d) > total_pred_chs(%d)", pred_channel_index, total_pred_chs)<<endl; exit(1); }  
  
  for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ ) {
    TH1D *h1d_temp = (TH1D*)map_default_h1d_pred[pred_channel_index]->Clone("h1d_temp"); h1d_temp->Reset();// define and clear
    
    for(int ievent=0; ievent<(int)vec_vec_eventinfo->at(isize).size(); ievent++) {
      EventInfo info = vec_vec_eventinfo->at(isize).at(ievent);
      double prob = Prob_oscillaion(info.e2e_Etrue, info.e2e_baseline, str_osc_mode);
      h1d_temp->Fill(info.e2e_Ereco, prob * info.e2e_weight_xs);
    }

    h1d_temp->Scale( vec_ratioPOT->at(isize) );

    ///////    
    TMatrixD matrix_temp(1, default_oldworld_rows);        
    int bin_index_base = 0;
    if( pred_channel_index > 1 ) { for(int ich=1; ich<pred_channel_index; ich++) bin_index_base += ( map_default_h1d_pred[ich]->GetNbinsX()+1 ); }
    for(int ibin=1; ibin<=h1d_temp->GetNbinsX()+1; ibin++) { matrix_temp(0, bin_index_base + ibin -1) = h1d_temp->GetBinContent(ibin); }
    matrix_oscillation_oldworld_pred += matrix_temp;    
  }// for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ )  
}

///////////////////
  
void TOsc::Apply_oscillation()
{
  //cout<<endl<<" ---> Apply_oscillation"<<endl;

  matrix_oscillation_oldworld_pred = matrix_oscillation_base_oldworld_pred;

  /////////////////// same as the ones in void TOsc::Set_oscillation_base(), but use "added" instead of "minus"

  if( flag_NuMI_nueCC_from_intnue ) {    
    Set_oscillation_base_added(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo, 15, "nue2nue");// hack
    Set_oscillation_base_added(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo, 16, "nue2nue");// hack
  }

  if( flag_NuMI_nueCC_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo, 15, "numu2numu");// hack
    Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo, 16, "numu2numu");// hack
  }
  
  if( flag_NuMI_numuCC_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo, 17, "numu2numu");// hack
    Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo, 18, "numu2numu");// hack
  }
  
  if( flag_NuMI_CCpi0_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo, 19, "numu2numu");// hack
    Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo, 20, "numu2numu");// hack
  }

  if( flag_NuMI_NCpi0_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo, 21, "numu2numu");// hack
  }
  
  /////////
  /////////
  /////////
  
  if( flag_BNB_nueCC_from_intnue ) {
    Set_oscillation_base_added(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo, 1, "nue2nue");// hack
    Set_oscillation_base_added(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo, 2, "nue2nue");// hack
  }
  
  if( flag_BNB_nueCC_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo, 1, "numu2numu");// hack
    Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo, 2, "numu2numu");// hack
  }
    
  if( flag_BNB_numuCC_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo, 3, "numu2numu");// hack
    Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo, 4, "numu2numu");// hack    
  }

  if( flag_BNB_CCpi0_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo, 5, "numu2numu");// hack
    Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo, 6, "numu2numu");// hack
  }

  if( flag_BNB_NCpi0_from_overlaynumu ) {
    Set_oscillation_base_added(&vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo, 7, "numu2numu");// hack
  }
  
  /////// winxp check with results from framework
  if( 0 ) {// self-check
    int idx_aa = 26*20;
    int idx_bb = 26*20;
    
    for(int idx=1; idx<=26; idx++) {
      cout<<TString::Format("%3d   %15.6f ---> (origin)%15.6f,    %15.6f  ---> (origin)%15.6f", idx,
			    matrix_oscillation_oldworld_pred(0, idx_aa  + idx-1), matrix_default_oldworld_pred(0, idx_aa  + idx-1),
			    matrix_oscillation_oldworld_pred(0, idx_bb  + idx-1), matrix_default_oldworld_pred(0, idx_bb  + idx-1)
			    )<<endl;
    }
  }// if( 0 ) {// self-check
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, TString str_treename, vector< vector<EventInfo> > *vec_vec_eventinfo)
{
  // Declaration of leaf types
  Int_t           e2e_pdg;
  //Int_t           e2e_flag_FC;
  Double_t        e2e_Etrue;
  Double_t        e2e_Ereco;
  Double_t        e2e_weight_xs;
  Double_t        e2e_baseline;

  // List of branches
  TBranch        *b_e2e_pdg;   //!
  //TBranch        *b_e2e_flag_FC;   //!
  TBranch        *b_e2e_Etrue;   //!
  TBranch        *b_e2e_Ereco;   //!
  TBranch        *b_e2e_weight_xs;   //!
  TBranch        *b_e2e_baseline;   //!
  
  TFile *roofile_obj = new TFile(strfile_mc_e2e, "read");
  TTree *tree_obj = (TTree*)roofile_obj->Get(str_treename);
  if(!tree_obj) { cerr<<" ERROR: Set_oscillation_base_subfunc, no treename: "<<str_treename<<endl; exit(1); }
  
  // Set branch addresses and branch pointers
  tree_obj->SetBranchAddress("e2e_pdg", &e2e_pdg, &b_e2e_pdg);
  //tree_obj->SetBranchAddress("e2e_flag_FC", &e2e_flag_FC, &b_e2e_flag_FC);
  tree_obj->SetBranchAddress("e2e_Etrue", &e2e_Etrue, &b_e2e_Etrue);
  tree_obj->SetBranchAddress("e2e_Ereco", &e2e_Ereco, &b_e2e_Ereco);
  tree_obj->SetBranchAddress("e2e_weight_xs", &e2e_weight_xs, &b_e2e_weight_xs);
  tree_obj->SetBranchAddress("e2e_baseline", &e2e_baseline, &b_e2e_baseline);
      
  int entries = tree_obj->GetEntries();

  if( vec_ratioPOT!=NULL ) cout<<endl;
  cout<<TString::Format("            ---> entries %10d     %50s   --> %20s", entries, strfile_mc_e2e.Data(), str_treename.Data())<<endl;

  vector<EventInfo>vector_eventinfo;
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_obj->GetEntry(ientry);
    EventInfo eventinfo;
    eventinfo.e2e_pdg = e2e_pdg;
    //eventinfo.e2e_flag_FC = e2e_flag_FC;
    eventinfo.e2e_Etrue = e2e_Etrue;
    eventinfo.e2e_Ereco = e2e_Ereco;
    eventinfo.e2e_weight_xs = e2e_weight_xs;
    eventinfo.e2e_baseline = e2e_baseline;
    vector_eventinfo.push_back(eventinfo);
  }
    
  delete tree_obj;
  delete roofile_obj;
  
  vec_vec_eventinfo->push_back( vector_eventinfo );

  if( vec_ratioPOT!=NULL ) {
    Double_t        pot;
    TBranch        *b_pot;   //!  
    double mc_pot = 1; double data_pot = 0;

    TFile *roofile_mc = new TFile(strfile_mcPOT, "read");
    TTree *tree_mc = (TTree*)roofile_mc->Get("T");
    tree_mc->SetBranchAddress("pot", &pot, &b_pot);
    tree_mc->GetEntry(0); mc_pot = pot;
    delete tree_mc;
    delete roofile_mc;
      
    TFile *roofile_data = new TFile(strfile_dataPOT, "read");
    TTree *tree_data = (TTree*)roofile_data->Get("T");
    tree_data->SetBranchAddress("pot", &pot, &b_pot);
    tree_data->GetEntry(0); data_pot = pot;
    delete tree_data;
    delete roofile_data;
    
    vec_ratioPOT->push_back( data_pot/mc_pot );// from the output of ./convert_histo.pl
    
    cout<<"            ---> MC POT "<<mc_pot<<"\t"<<strfile_mcPOT<<endl;
    cout<<"            ---> DD POT "<<data_pot<<"\t"<<strfile_dataPOT<<endl;      
  }// if( vec_ratioPOT!=NULL )
  
}

///////////////////
  
void TOsc::Set_oscillation_base()
{
  /// Eventlist generator
  /// /home/xji/data0/work/work_oscillation/301_Framework_for_Osc/wcp-uboone-bdt/apps/convert_checkout_hist.cxx, winxp
  
  cout<<endl;
  cout<<" ---> Set_oscillation_base"<<endl;

  matrix_oscillation_base_oldworld_pred = matrix_default_oldworld_pred;

  TString str_dirbase = "./data_inputs/newdir_BNBNuMI_cv_cov_list/";
  
  ///////////////////
  
  if( flag_NuMI_nueCC_from_intnue ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_intnue"<<endl;        
    {// run1 FHC
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo, 15, "nue2nue");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo, 16, "nue2nue");// hack
  
  }// if( flag_NuMI_nueCC_from_intnue )

  ///////////////////
  ///////////////////  

  if( flag_NuMI_nueCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_overlaynumu"<<endl;        
    {// run1 FHC
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_numi_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo, 15, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo, 16, "numu2numu");// hack
  
  }// if( flag_NuMI_nueCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_NuMI_numuCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_numuCC_from_overlaynumu"<<endl;        
    {// run1 FHC
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_numi_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo, 17, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo, 18, "numu2numu");// hack
  
  }// if( flag_NuMI_numuCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_NuMI_CCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_CCpi0_from_overlaynumu"<<endl;        
    {// run1 FHC
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_numi_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo, 19, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo, 20, "numu2numu");// hack
  
  }// if( flag_NuMI_CCpi0_from_overlaynumu )
  
  ///////////////////  

  if( flag_NuMI_NCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_NCpi0_from_overlaynumu"<<endl;        
    {// run1 FHC
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_numi_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo, 21, "numu2numu");// hack
  
  }// if( flag_NuMI_NCpi0_from_overlaynumu )

  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    

  if( flag_BNB_nueCC_from_intnue ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_intnue"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo, 1, "nue2nue");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo, 2, "nue2nue");// hack
    
  }// if( flag_BNB_nueCC_from_intnue )

  ///////////////////
  ///////////////////  

  if( flag_BNB_nueCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo, 1, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo, 2, "numu2numu");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_BNB_numuCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_numuCC_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo, 3, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo, 4, "numu2numu");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynumu )

  ///////////////////  

  if( flag_BNB_CCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_CCpi0_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo, 5, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo, 6, "numu2numu");// hack
    
  }// if( flag_BNB_CCpi0_from_overlaynumu )
  
  ///////////////////  

  if( flag_BNB_NCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_NCpi0_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo, 7, "numu2numu");// hack
    
  }// if( flag_BNB_NCpi0_from_overlaynumu )
  
  cout<<endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_default_cv_cov(TString default_cv_file, TString default_dirtadd_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir)
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_default_cv_cov"<<endl;
  
  cout<<endl;
  cout<<TString::Format("      ---> default_cv_file       %10s", default_cv_file.Data())<<endl;
  cout<<TString::Format("      ---> default_dirtadd_file  %10s", default_dirtadd_file.Data())<<endl;
  cout<<TString::Format("      ---> default_mcstat_file   %10s", default_mcstat_file.Data())<<endl;
  cout<<TString::Format("      ---> default_fluxXs_dir    %10s", default_fluxXs_dir.Data())<<endl;
  cout<<TString::Format("      ---> default_detector_dir  %10s", default_detector_dir.Data())<<endl;

  //////////////////////////////////////

  {
    TFile *roofile_default_cv_file = new TFile(default_cv_file, "read");

    ///////
    cout<<endl;
    cout<<"      ---> matrix_transform"<<endl;
    TMatrixD *temp_mat_collapse = (TMatrixD*)roofile_default_cv_file->Get("mat_collapse");  
    default_oldworld_rows = temp_mat_collapse->GetNrows();
    default_newworld_rows = temp_mat_collapse->GetNcols();
    matrix_transform.Clear(); matrix_transform.ResizeTo(default_oldworld_rows, default_newworld_rows);
    matrix_transform += (*temp_mat_collapse);
    delete temp_mat_collapse;
    
    /// hack, set LEE channel to 0
    // for(int idx=0; idx<default_oldworld_rows; idx++) {
    //   for(int jdx=0; jdx<default_newworld_rows; jdx++) {
    // 	if(idx>=26*4+11*3 && idx<26*4+11*3+26*2) matrix_transform(idx, jdx) = 0;
    //   }
    // }
    
    cout<<"      ---> default_oldworld_rows  "<<default_oldworld_rows<<endl;
    cout<<"      ---> default_newworld_rows  "<<default_newworld_rows<<endl;

    ///////

    matrix_default_newworld_meas.Clear(); matrix_default_newworld_meas.ResizeTo(1, default_newworld_rows);
    matrix_default_oldworld_pred.Clear(); matrix_default_oldworld_pred.ResizeTo(1, default_oldworld_rows);

    matrix_oscillation_base_oldworld_pred.Clear(); matrix_oscillation_base_oldworld_pred.ResizeTo(1, default_oldworld_rows);
    matrix_oscillation_oldworld_pred.Clear();      matrix_oscillation_oldworld_pred.ResizeTo(1, default_oldworld_rows);

    ///////
    
    matrix_default_oldworld_abs_syst_addi.Clear();   matrix_default_oldworld_abs_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  
    matrix_default_oldworld_rel_syst_flux.Clear();   matrix_default_oldworld_rel_syst_flux.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_geant.Clear();  matrix_default_oldworld_rel_syst_geant.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_Xs.Clear();     matrix_default_oldworld_rel_syst_Xs.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_det.Clear();    matrix_default_oldworld_rel_syst_det.ResizeTo(default_oldworld_rows, default_oldworld_rows);

    matrix_default_newworld_abs_syst_mcstat.Clear(); matrix_default_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    
    //////
    
    matrix_eff_newworld_abs_syst_addi.Clear();   matrix_eff_newworld_abs_syst_addi.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_mcstat.Clear(); matrix_eff_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_flux.Clear();   matrix_eff_newworld_abs_syst_flux.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_geant.Clear();  matrix_eff_newworld_abs_syst_geant.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_Xs.Clear();     matrix_eff_newworld_abs_syst_Xs.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_det.Clear();    matrix_eff_newworld_abs_syst_det.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_eff_newworld_abs_syst_total.Clear();  matrix_eff_newworld_abs_syst_total.ResizeTo(default_newworld_rows, default_newworld_rows);
	  
    matrix_eff_newworld_meas.Clear();  matrix_eff_newworld_meas.ResizeTo(1, default_newworld_rows);
    matrix_eff_newworld_pred.Clear();  matrix_eff_newworld_pred.ResizeTo(1, default_newworld_rows);
    matrix_eff_newworld_noosc.Clear(); matrix_eff_newworld_noosc.ResizeTo(1, default_newworld_rows);

    matrix_fitdata_newworld.Clear(); matrix_fitdata_newworld.ResizeTo(1, default_newworld_rows);
    
    ///////
    cout<<endl;
    cout<<"      ---> measurement"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("hdata_obsch_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %3d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      int bins = h1f_temp->GetNbinsX(); double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1); double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_default_h1d_meas_bins[idx] = bins; map_default_h1d_meas_xlow[idx] = xlow; map_default_h1d_meas_xhgh[idx] = xhgh;                 
      for(int ibin=1; ibin<=bins+1; ibin++) vector_default_newworld_meas.push_back( h1f_temp->GetBinContent(ibin) );      
      delete h1f_temp;    
    }
    
    if( default_newworld_rows!=(int)(vector_default_newworld_meas.size()) ) { cerr<<" ---> ERROR: default_newworld_rows!=vector_default_newworld_meas"<<endl; exit(1); }
    for(int idx=0; idx<(int)(vector_default_newworld_meas.size()); idx++)  matrix_default_newworld_meas(0, idx) = vector_default_newworld_meas.at(idx);
        
    ///////
    cout<<endl;
    cout<<"      ---> prediction"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("histo_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %3d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      int bins = h1f_temp->GetNbinsX(); double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1); double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_default_h1d_pred_bins[idx] = bins; map_default_h1d_pred_xlow[idx] = xlow; map_default_h1d_pred_xhgh[idx] = xhgh;                 
      for(int ibin=1; ibin<=bins+1; ibin++) vector_default_oldworld_pred.push_back( h1f_temp->GetBinContent(ibin) );      
      delete h1f_temp;    
    }
    
    if( default_oldworld_rows!=(int)(vector_default_oldworld_pred.size()) ) { cerr<<" ---> ERROR: default_oldworld_rows!=vector_default_oldworld_pred"<<endl; exit(1); }
    for(int idx=0; idx<(int)(vector_default_oldworld_pred.size()); idx++)  matrix_default_oldworld_pred(0, idx) = vector_default_oldworld_pred.at(idx);
    
    ///////
    
    delete roofile_default_cv_file;
  }

  
  {
    for(auto it=map_default_h1d_meas_bins.begin(); it!=map_default_h1d_meas_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("default_h1d_meas_%d", idx);
      map_default_h1d_meas[idx] = new TH1D(roostr, roostr, it->second, map_default_h1d_meas_xlow[idx], map_default_h1d_meas_xhgh[idx]);
    }
    
    for(auto it=map_default_h1d_pred_bins.begin(); it!=map_default_h1d_pred_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("default_h1d_pred_%d", idx);
      map_default_h1d_pred[idx] = new TH1D(roostr, roostr, it->second, map_default_h1d_pred_xlow[idx], map_default_h1d_pred_xhgh[idx]);      
    }

    TFile *roofile_default_cv_file = new TFile(default_cv_file, "read");

    for(auto it=map_default_h1d_meas_bins.begin(); it!=map_default_h1d_meas_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("hdata_obsch_%d", idx);
      int bins = it->second;
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_default_h1d_meas[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    for(auto it=map_default_h1d_pred_bins.begin(); it!=map_default_h1d_pred_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("histo_%d", idx);
      int bins = it->second;
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_default_h1d_pred[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    delete roofile_default_cv_file;
  }
  
  //////////////////////////////////////

  if( flag_syst_dirt ) {
    cout<<endl;
    cout<<"      ---> Dirt: additional uncertainty, Yes"<<endl;
    
    TFile *roofile_default_dirtadd_file = new TFile(default_dirtadd_file, "read");
    TMatrixD *matrix_temp = (TMatrixD*)roofile_default_dirtadd_file->Get("cov_mat_add");
    matrix_default_oldworld_abs_syst_addi = (*matrix_temp);
    delete matrix_temp;
    delete roofile_default_dirtadd_file;
  }
  else {
    cout<<endl;
    cout<<"      ---> Dirt: additional uncertainty, No"<<endl;
  }
  
  //////////////////////////////////////

  if( flag_syst_mcstat ) {
    cout<<endl;
    cout<<"      ---> MCstat, Yes"<<endl;
    
    ifstream file_mcstat_aa(default_mcstat_file, ios::in);
    if(!file_mcstat_aa) { cerr<<" Error: No file_mcstat_aa"<<endl; exit(1); }
    int count_aa = 0;
    string str_count_aa;    
    ifstream file_check_aa(default_mcstat_file);
    while( getline(file_check_aa, str_count_aa) ) count_aa++;
    if( count_aa-1 != default_newworld_rows ) { cerr<<" Error: mcstat != default_newworld_rows"<<endl; exit(1);  }

    ifstream file_mcstat_bb(default_mcstat_file, ios::in);
    double Lee = 1; double run = 1;
    file_mcstat_bb >> Lee >> run;
    int gch = 0; int lbin = 0; double val_pred = 0; double mc_stat = 0; double nn_stat = 0;
    for(int idx=1; idx<=default_newworld_rows; idx++) {
      file_mcstat_bb >> gch >> lbin >> val_pred >> mc_stat >> nn_stat;
      matrix_default_newworld_abs_syst_mcstat(idx-1, idx-1) = mc_stat;
    }
  }
  else {
    cout<<endl;
    cout<<"      ---> MCstat, No"<<endl;
  }
  
  ////////////////////////////////////// flux, geant, Xs

  if( flag_syst_flux ) {
    cout<<endl;
    cout<<"      ---> flux, Yes"<<endl;
    
    for(int idx=1; idx<=13; idx++) {
      TFile *roofile_syst_flux = new TFile(default_fluxXs_dir+TString::Format("cov_%d.root", idx), "read");
      TMatrixD *matrix_temp_flux = (TMatrixD*)roofile_syst_flux->Get(TString::Format("frac_cov_xf_mat_%d", idx) );
      matrix_default_oldworld_rel_syst_flux += (*matrix_temp_flux);
      delete matrix_temp_flux;
      delete roofile_syst_flux;
    }
	
  }// if( flag_syst_flux )
  else {
    cout<<endl;
    cout<<"      ---> flux, No"<<endl;
  }
  
  if( flag_syst_geant ) {
    cout<<endl;
    cout<<"      ---> geant, Yes"<<endl;

    for(int idx=14; idx<=16; idx++) {
      TFile *roofile_syst_geant = new TFile(default_fluxXs_dir+TString::Format("cov_%d.root", idx), "read");
      TMatrixD *matrix_temp_geant = (TMatrixD*)roofile_syst_geant->Get(TString::Format("frac_cov_xf_mat_%d", idx) );
      matrix_default_oldworld_rel_syst_geant += (*matrix_temp_geant);
      delete matrix_temp_geant;
      delete roofile_syst_geant;
    }   
  }// if( flag_syst_geant )
  else {
    cout<<endl;
    cout<<"      ---> geant, No"<<endl;
  }
  
  if( flag_syst_Xs ) {
    cout<<endl;
    cout<<"      ---> Xs, Yes"<<endl;
    
    TFile *roofile_syst_Xs = new TFile(default_fluxXs_dir+"cov_17.root", "read");
    TMatrixD *matrix_temp_Xs = (TMatrixD*)roofile_syst_Xs->Get("frac_cov_xf_mat_17");
    matrix_default_oldworld_rel_syst_Xs = (*matrix_temp_Xs);
    delete matrix_temp_Xs;
    delete roofile_syst_Xs; 
  }// if( flag_syst_Xs )
  else {
    cout<<endl;
    cout<<"      ---> Xs, No"<<endl;
  }
  
  ////////////////////////////////////// detector

  if( flag_syst_det ) {
    cout<<endl;
    cout<<"      ---> detector, Yes"<<endl;
 
    map<int, TString>map_detectorfile_str;    
    map_detectorfile_str[1] = "cov_LYDown.root";
    map_detectorfile_str[2] = "cov_LYRayleigh.root";
    map_detectorfile_str[3] = "cov_Recomb2.root";
    map_detectorfile_str[4] = "cov_SCE.root";
    //map_detectorfile_str[5] = "cov_WMdEdx.root";
    map_detectorfile_str[6] = "cov_WMThetaXZ.root";
    map_detectorfile_str[7] = "cov_WMThetaYZ.root";
    map_detectorfile_str[8] = "cov_WMX.root";
    map_detectorfile_str[9] = "cov_WMYZ.root";
    map_detectorfile_str[10]= "cov_LYatt.root";
    
    for(auto it_map=map_detectorfile_str.begin(); it_map!=map_detectorfile_str.end(); it_map++) {
      int idx = it_map->first;
      TFile *roofile_det = new TFile(default_detector_dir + map_detectorfile_str[idx], "read");
      TMatrixD *matrix_temp_det = (TMatrixD*)roofile_det->Get(TString::Format("frac_cov_det_mat_%d", idx) );
      matrix_default_oldworld_rel_syst_det += (*matrix_temp_det);
      delete matrix_temp_det;
      delete roofile_det;
    }
  }// if( flag_syst_det )
  else {
    cout<<endl;
    cout<<"      ---> detector, No"<<endl;
  }
  
}
  
