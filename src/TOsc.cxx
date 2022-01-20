#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

int TOsc::Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar)
{
  int flag_minimization = 0;

  /////////////////////////////////////////// ddd

  {
    ROOT::Minuit2::Minuit2Minimizer min_osc( ROOT::Minuit2::kMigrad );
    min_osc.SetPrintLevel(2);
    min_osc.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
    min_osc.SetMaxFunctionCalls(500000);
    min_osc.SetMaxIterations(500000);
    min_osc.SetTolerance(1e-4); // tolerance*2e-3 = edm precision
    min_osc.SetPrecision(1e-18); //precision in the target function

    /// set fitting parameters
    ROOT::Math::Functor Chi2Functor_osc( [&](const double *par) {// FCN
	
	double chi2_final = 0;
	double fit_dm2_41         = par[0];
	double fit_sin2_2theta_14 = par[1];
	
	/////// standard order
	Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, init_sin2_theta_24, init_sin2_theta_34);  
	Apply_oscillation();
	Set_apply_POT();// meas, CV, COV: all ready

	///////
	
	TMatrixD matrix_data_total = matrix_meas2fitdata_newworld;
	TMatrixD matrix_pred_total = matrix_default_newworld_pred;
	
	TMatrixD matrix_cov_syst_total = matrix_default_newworld_abs_syst_total;
	int rows = matrix_cov_syst_total.GetNrows();

	TMatrixD matrix_cov_stat_total(rows, rows);
	TMatrixD matrix_cov_total(rows, rows);

	for(int idx=0; idx<=6; idx++) {
	  double val_data = matrix_data_total(0, idx);
	  double val_pred = matrix_pred_total(0, idx);       
	  //cout<<TString::Format(" %3d  %8.4f %8.4f", idx+1, val_data, val_pred)<<endl;
	}
	
	for(int idx=0; idx<rows; idx++) {
	  double val_stat_cov = 0;        
	  double val_data = matrix_data_total(0, idx);
	  double val_pred = matrix_pred_total(0, idx);
        
	  if( val_data==0 ) { val_stat_cov = val_pred/2; }
	  else {
	    if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
	    else val_stat_cov = val_data;
	  }
	
	  matrix_cov_stat_total(idx, idx) += val_stat_cov;
	  if( matrix_cov_syst_total(idx, idx)==0 ) matrix_cov_syst_total(idx, idx) = 1e-6;
	}

	matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
	
	///////
	
	TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
	for(int idx=0; idx<rows; idx++) {
	  if( matrix_cov_total(idx, idx)<=2e-6 ) cout<<idx<<"\t"<<matrix_cov_total(idx, idx)<<endl;
	}
      
	TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;
	TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();
   
	TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
	chi2_final = matrix_chi2(0,0);           
      
	///////

	return chi2_final;
	
      },// end of FCN
      2// number of fitting parameters
      );


    
    min_osc.SetFunction(Chi2Functor_osc);
    
    min_osc.SetVariable( 0, "dm2_41", init_dm2_41, 1e-2);
    min_osc.SetVariable( 1, "sin2_2theta_14", init_sin2_2theta_14, 1e-2);

    //min_osc.SetLowerLimitedVariable(0, "dm2_41", init_dm2_41, 1e-3, 0);  
    //min_osc.SetLimitedVariable(1, "sin2_2theta_14", init_sin2_2theta_14, 1e-3, 0, 1);
    
    min_osc.Minimize();
    
  }
  
  ///////////////////////////////////////////
  
  return flag_minimization;
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_apply_POT()
{
  ////////////////////////// hack for POT scale, should know which part is BNB, and which part is NuMI
  
  matrix_default_oldworld_abs_syst_addi_POT = matrix_default_oldworld_abs_syst_addi;
  matrix_default_newworld_abs_syst_mcstat_POT = matrix_default_newworld_abs_syst_mcstat;
  matrix_default_newworld_meas_POT = matrix_default_newworld_meas;
  
  /////// hack
  
  for(int idx=0; idx<default_oldworld_rows; idx++) {
    matrix_oscillation_effect_result_oldworld(0, idx) *= scaleF_POT_BNB;
    matrix_default_oldworld_abs_syst_addi_POT(idx, idx) *= (scaleF_POT_BNB*scaleF_POT_BNB);
  }

  /////// hack
    
  for(int idx=0; idx<default_newworld_rows; idx++) {
    matrix_default_newworld_abs_syst_mcstat_POT(idx, idx) *= (scaleF_POT_BNB*scaleF_POT_BNB);
    matrix_default_newworld_meas_POT(0, idx) *= scaleF_POT_BNB;
  }
  
  matrix_default_newworld_pred = matrix_oscillation_effect_result_oldworld * matrix_transform;
  
  //////////////////////////
  
  for(int idx=0; idx<default_oldworld_rows; idx++) {        
    for(int jdx=0; jdx<default_oldworld_rows; jdx++) {
      double cv_i = matrix_oscillation_effect_result_oldworld(0, idx);
      double cv_j = matrix_oscillation_effect_result_oldworld(0, jdx);

      /// flux
      matrix_default_oldworld_abs_syst_flux(idx,jdx) = cv_i*cv_j*matrix_default_oldworld_rel_syst_flux(idx,jdx);

      /// geant
      matrix_default_oldworld_abs_syst_geant(idx,jdx) = cv_i*cv_j*matrix_default_oldworld_rel_syst_geant(idx,jdx);
      
      /// Xs
      matrix_default_oldworld_abs_syst_Xs(idx,jdx) = cv_i*cv_j*matrix_default_oldworld_rel_syst_Xs(idx,jdx);
            
      /// det
      matrix_default_oldworld_abs_syst_det(idx,jdx) = cv_i*cv_j*matrix_default_oldworld_rel_syst_det(idx,jdx);      
    }
  }

  //////////////////////////

  TMatrixD matrix_transform_T =  matrix_transform.T();  matrix_transform.T();
  
  /// dirt addi
  matrix_default_newworld_abs_syst_addi = matrix_transform_T * matrix_default_oldworld_abs_syst_addi_POT * matrix_transform;

  /// flux
  matrix_default_newworld_abs_syst_flux = matrix_transform_T * matrix_default_oldworld_abs_syst_flux * matrix_transform;
 
  /// geant
  matrix_default_newworld_abs_syst_geant = matrix_transform_T * matrix_default_oldworld_abs_syst_geant * matrix_transform;
  
  /// Xs
  matrix_default_newworld_abs_syst_Xs = matrix_transform_T * matrix_default_oldworld_abs_syst_Xs * matrix_transform;
  
  /// det
  matrix_default_newworld_abs_syst_det = matrix_transform_T * matrix_default_oldworld_abs_syst_det * matrix_transform;

  ///
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_addi;
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_flux;
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_geant;
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_Xs;
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_det;
  matrix_default_newworld_abs_syst_total += matrix_default_newworld_abs_syst_mcstat_POT;

  ///
  if( 0 ) {
    TString roostr = "";

    TF1 *f1_one = new TF1("f1_one", "1", 0, 1e6);
    f1_one->SetLineColor(kBlack);
    f1_one->SetLineStyle(7);
    
    roostr = "h1_newworld_pred_test";
    TH1D *h1_newworld_pred_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_pred_relerr_test";
    TH1D *h1_newworld_pred_relerr_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_meas_test";
    TH1D *h1_newworld_meas_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    roostr = "h1_newworld_meas2pred_test";
    TH1D *h1_newworld_meas2pred_test = new TH1D(roostr, "", default_newworld_rows, 0, default_newworld_rows);
    
    for(int ibin=1; ibin<default_newworld_rows; ibin++) {
      double content = matrix_default_newworld_pred(0, ibin-1);
      double abs_err = sqrt( matrix_default_newworld_abs_syst_total(ibin, ibin) );
      h1_newworld_pred_test->SetBinContent(ibin, content);
      h1_newworld_pred_test->SetBinError(ibin, abs_err);

      double relerr = 0;
      if(content!=0) relerr = abs_err/content;
      h1_newworld_pred_relerr_test->SetBinError(ibin, relerr);
      h1_newworld_pred_relerr_test->SetBinContent(ibin, 1);

      double meas = matrix_default_newworld_meas_POT(0, ibin-1);
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
  }
}
 
//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

double TOsc::Prob_oscillaion(double Etrue, double baseline, TString strflag_osc) {
  int flag_osc = -1;
  double prob = 0;
    
  if( strflag_osc=="nue2nue" )        flag_osc = 1;
  else if( strflag_osc=="numu2numu" ) flag_osc = 2;
  else if( strflag_osc=="numu2nue" )  flag_osc = 3;
  else if( strflag_osc=="nue2numu" )  flag_osc = 4;
  else flag_osc = -1;

  if( flag_osc==-1 ) {
    cerr<<" ERROR: flag_osc==-1, strflag_osc=="<<strflag_osc<<" ---> Candidates: nue2nue numu2numu numu2nue nue2numu"<<endl;
    exit(1);
  }

  const int nue2nue   = 1;
  const int numu2numu = 2;
  const int numu2nue  = 3;
  const int nue2numu  = 4;

  /// dm2*L/4E = 1.267 * dm2(eV2) * L(m) / E(MeV)
    
  switch( flag_osc ) {
  case nue2nue:
    prob = 1 - sin2_2theta_14 * pow(TMath::Sin(1.267 * dm2_41 * baseline/Etrue), 2);
    break;
  case numu2numu:
    break;
  case numu2nue:
    break;
  case nue2numu:
    break;
  default:
    cerr<<"ERROR: NAN flag_osc"<<endl; exit(1);
  }
    
  return prob;
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Apply_oscillation()
{
  //cout<<endl<<" ---> Apply_oscillation"<<endl;  
  
  matrix_oscillation_effect_result_oldworld.Clear(); matrix_oscillation_effect_result_oldworld.ResizeTo(1, default_oldworld_rows);
  matrix_oscillation_effect_result_oldworld = matrix_oscillation_base;
    
  //////////////////////////////////////

  if( flag_NuMI_nue2nue ) {
            
    int size = vector_vector_oscillation_base_NuMI_nueCC_info.size();
    for(int isize=0; isize<size; isize++) {
      
      TMatrixD matrix_oscillation_base_add(1, default_oldworld_rows);
      
      TH1D *h1_FC = (TH1D*)map_default_h1d_pred[1]->Clone("h1_FC"); h1_FC->Reset();// define and clear
      TH1D *h1_PC = (TH1D*)map_default_h1d_pred[2]->Clone("h1_PC"); h1_PC->Reset();// define and clear
      
      double scaleFPOT = vector_oscillation_base_NuMI_nueCC_scaleFPOT.at(isize);
      int num_events = vector_vector_oscillation_base_NuMI_nueCC_info.at(isize).size();
      
      for(int ievent=0; ievent<num_events; ievent++) {
	EventInfo info = vector_vector_oscillation_base_NuMI_nueCC_info.at(isize).at(ievent);

	double prob = Prob_oscillaion(info.e2e_Etrue, info.e2e_baseline, "nue2nue");
	
	if( info.e2e_flag_FC ) h1_FC->Fill(info.e2e_Ereco, prob * info.e2e_weight_xs);
	else h1_PC->Fill(info.e2e_Ereco, prob * info.e2e_weight_xs);
      }
      
      h1_FC->Scale( scaleFPOT );
      h1_PC->Scale( scaleFPOT );
      
      /// hack for NuMI nueCC
      for(int ibin=1; ibin<=h1_FC->GetNbinsX()+1; ibin++) matrix_oscillation_base_add(0, ibin-1) = h1_FC->GetBinContent(ibin);
      for(int ibin=1; ibin<=h1_PC->GetNbinsX()+1; ibin++) matrix_oscillation_base_add(0, (h1_FC->GetNbinsX()+1) +ibin-1) = h1_PC->GetBinContent(ibin);
      matrix_oscillation_effect_result_oldworld += matrix_oscillation_base_add;
      
      delete h1_FC;
      delete h1_PC;      
    }// for(int isize=0; isize<size; isize++)
        
  }
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, vector< vector<EventInfo> > *vec_vec_eventinfo)
{
  // Declaration of leaf types
  Int_t           e2e_pdg;
  Int_t           e2e_flag_FC;
  Double_t        e2e_Etrue;
  Double_t        e2e_Ereco;
  Double_t        e2e_weight_xs;
  Double_t        e2e_baseline;

  // List of branches
  TBranch        *b_e2e_pdg;   //!
  TBranch        *b_e2e_flag_FC;   //!
  TBranch        *b_e2e_Etrue;   //!
  TBranch        *b_e2e_Ereco;   //!
  TBranch        *b_e2e_weight_xs;   //!
  TBranch        *b_e2e_baseline;   //!
  
  TFile *roofile_obj = new TFile(strfile_mc_e2e, "read");
  TTree *tree_obj = (TTree*)roofile_obj->Get("tree_obj");
  
  // Set branch addresses and branch pointers
  tree_obj->SetBranchAddress("e2e_pdg", &e2e_pdg, &b_e2e_pdg);
  tree_obj->SetBranchAddress("e2e_flag_FC", &e2e_flag_FC, &b_e2e_flag_FC);
  tree_obj->SetBranchAddress("e2e_Etrue", &e2e_Etrue, &b_e2e_Etrue);
  tree_obj->SetBranchAddress("e2e_Ereco", &e2e_Ereco, &b_e2e_Ereco);
  tree_obj->SetBranchAddress("e2e_weight_xs", &e2e_weight_xs, &b_e2e_weight_xs);
  tree_obj->SetBranchAddress("e2e_baseline", &e2e_baseline, &b_e2e_baseline);
      
  int entries = tree_obj->GetEntries();

  cout<<endl;
  cout<<TString::Format("      ---> entries %10d   %10s", entries, strfile_mc_e2e.Data())<<endl;

  vector<EventInfo>vector_eventinfo;
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_obj->GetEntry(ientry);
    EventInfo eventinfo;
    eventinfo.e2e_pdg = e2e_pdg;
    eventinfo.e2e_flag_FC = e2e_flag_FC;
    eventinfo.e2e_Etrue = e2e_Etrue;
    eventinfo.e2e_Ereco = e2e_Ereco;
    eventinfo.e2e_weight_xs = e2e_weight_xs;
    eventinfo.e2e_baseline = e2e_baseline;
    vector_eventinfo.push_back(eventinfo);
  }
    
  delete tree_obj;
  delete roofile_obj;

  //
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

  cout<<"      ---> POT "<<mc_pot<<"\t"<<strfile_mcPOT<<endl;
  cout<<"      ---> POT "<<data_pot<<"\t"<<strfile_dataPOT<<endl;
  
  vec_ratioPOT->push_back( data_pot/mc_pot );// from the output of ./convert_histo.pl
  vec_vec_eventinfo->push_back( vector_eventinfo );
  
}
 
/////////////////////////////////////////////// ccc

void TOsc::Set_oscillation_base()
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_oscillation_base: must appear only one time"<<endl;
  
  //////////////////////////////////////

  matrix_oscillation_base.Clear();  matrix_oscillation_base.ResizeTo(1, default_oldworld_rows);
  matrix_oscillation_base = matrix_default_oldworld_pred;
  
  //////////////////////////////////////

  if( flag_NuMI_nue2nue ) {
    
    {// intrinsic run1     
      TString strfile_mcPOT = "./data_inputs/hist_rootfiles_default_noosc/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root";
      TString strfile_dataPOT = "./data_inputs/hist_rootfiles_default_noosc/run1_data_bnb.root";
      TString strfile_mc_e2e = "./data_inputs/hist_rootfiles_default_noosc/roofile_obj_NuMI_run1_intrinsic_nue.root";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_oscillation_base_NuMI_nueCC_scaleFPOT, strfile_mc_e2e, &vector_vector_oscillation_base_NuMI_nueCC_info);
    }

    ///////////////////
            
    int size = vector_vector_oscillation_base_NuMI_nueCC_info.size();
    for(int isize=0; isize<size; isize++) {
      
      TMatrixD matrix_oscillation_base_subtract(1, default_oldworld_rows);
      
      TH1D *h1_FC = (TH1D*)map_default_h1d_pred[1]->Clone("h1_FC"); h1_FC->Reset();// define and clear
      TH1D *h1_PC = (TH1D*)map_default_h1d_pred[2]->Clone("h1_PC"); h1_PC->Reset();// define and clear
    
      double scaleFPOT = vector_oscillation_base_NuMI_nueCC_scaleFPOT.at(isize);
      int num_events = vector_vector_oscillation_base_NuMI_nueCC_info.at(isize).size();

      for(int ievent=0; ievent<num_events; ievent++) {
	EventInfo info = vector_vector_oscillation_base_NuMI_nueCC_info.at(isize).at(ievent);
	if( info.e2e_flag_FC ) h1_FC->Fill(info.e2e_Ereco, info.e2e_weight_xs);
	else h1_PC->Fill(info.e2e_Ereco, info.e2e_weight_xs);
      }

      h1_FC->Scale( scaleFPOT );
      h1_PC->Scale( scaleFPOT );

      /// hack for NuMI nueCC
      for(int ibin=1; ibin<=h1_FC->GetNbinsX()+1; ibin++) matrix_oscillation_base_subtract(0, ibin-1) = h1_FC->GetBinContent(ibin);
      for(int ibin=1; ibin<=h1_PC->GetNbinsX()+1; ibin++) matrix_oscillation_base_subtract(0, (h1_FC->GetNbinsX()+1) +ibin-1) = h1_PC->GetBinContent(ibin);
      matrix_oscillation_base -= matrix_oscillation_base_subtract;
      
      delete h1_FC;
      delete h1_PC;      
    }// for(int isize=0; isize<size; isize++)
    
  }// if( flag_NuMI_nue2nue )
  
  if( flag_NuMI_numu2numu ) {
    
  }// if( flag_NuMI_numu2numu )

  if( flag_NuMI_numu2nue ) {
    
  }// if( flag_NuMI_numu2nue )

  if( flag_NuMI_nue2numu ) {
    
  }// if( flag_NuMI_nue2numu )

  if( flag_NuMI_NC_1minus_nue2sterile ) {

  }// if( flag_NuMI_NC_1minus_nue2sterile )

  if( flag_NuMI_NC_1minus_numu2sterile ) {

  }// if( flag_NuMI_NC_1minus_numu2sterile )
  
  //////////////////////////////////////

  if( flag_BNB_nue2nue ) {
    
  }// if( flag_BNB_nue2nue )

  if( flag_BNB_numu2numu ) {
    
  }// if( flag_BNB_numu2numu )

  if( flag_BNB_numu2nue ) {
    
  }// if( flag_BNB_numu2nue )

  if( flag_BNB_nue2numu ) {
    
  }// if( flag_BNB_nue2numu )

  if( flag_BNB_NC_1minus_nue2sterile ) {

  }// if( flag_BNB_NC_1minus_nue2sterile )

  if( flag_BNB_NC_1minus_numu2sterile ) {

  }// if( flag_BNB_NC_1minus_numu2sterile )
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_default_cv_cov(TString default_cv_file, TString default_dirtadd_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir)
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_default_cv_cov"<<endl;

  cout<<endl;
  cout<<TString::Format("      ---> default_cv_file       %50s", default_cv_file.Data())<<endl;
  cout<<TString::Format("      ---> default_dirtadd_file  %50s", default_dirtadd_file.Data())<<endl;
  cout<<TString::Format("      ---> default_mcstat_file   %50s", default_mcstat_file.Data())<<endl;
  cout<<TString::Format("      ---> default_fluxXs_dir    %50s", default_fluxXs_dir.Data())<<endl;
  cout<<TString::Format("      ---> default_detector_dir  %50s", default_detector_dir.Data())<<endl;

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
    for(int idx=0; idx<default_oldworld_rows; idx++) {
      for(int jdx=0; jdx<default_newworld_rows; jdx++) {
	if(idx>=26*4+11*3 && idx<26*4+11*3+26*2) matrix_transform(idx, jdx) = 0;
      }
    }
    
    ///////
    matrix_default_oldworld_abs_syst_addi.Clear();   matrix_default_oldworld_abs_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_addi_POT.Clear();   matrix_default_oldworld_abs_syst_addi_POT.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_mcstat.Clear(); matrix_default_oldworld_abs_syst_mcstat.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_flux.Clear();   matrix_default_oldworld_abs_syst_flux.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_geant.Clear();  matrix_default_oldworld_abs_syst_geant.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_Xs.Clear();     matrix_default_oldworld_abs_syst_Xs.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_det.Clear();    matrix_default_oldworld_abs_syst_det.ResizeTo(default_oldworld_rows, default_oldworld_rows);

    matrix_default_newworld_abs_syst_addi.Clear();   matrix_default_newworld_abs_syst_addi.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_mcstat.Clear(); matrix_default_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_mcstat_POT.Clear(); matrix_default_newworld_abs_syst_mcstat_POT.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_flux.Clear();   matrix_default_newworld_abs_syst_flux.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_geant.Clear();  matrix_default_newworld_abs_syst_geant.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_Xs.Clear();     matrix_default_newworld_abs_syst_Xs.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_det.Clear();    matrix_default_newworld_abs_syst_det.ResizeTo(default_newworld_rows, default_newworld_rows); 
  
    ///////
    matrix_default_oldworld_rel_syst_addi.Clear();   matrix_default_oldworld_rel_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_mcstat.Clear(); matrix_default_oldworld_rel_syst_mcstat.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_flux.Clear();   matrix_default_oldworld_rel_syst_flux.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_geant.Clear();  matrix_default_oldworld_rel_syst_geant.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_Xs.Clear();     matrix_default_oldworld_rel_syst_Xs.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_rel_syst_det.Clear();    matrix_default_oldworld_rel_syst_det.ResizeTo(default_oldworld_rows, default_oldworld_rows);

    matrix_default_newworld_rel_syst_addi.Clear();   matrix_default_newworld_rel_syst_addi.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_rel_syst_mcstat.Clear(); matrix_default_newworld_rel_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_rel_syst_flux.Clear();   matrix_default_newworld_rel_syst_flux.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_rel_syst_geant.Clear();  matrix_default_newworld_rel_syst_geant.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_rel_syst_Xs.Clear();     matrix_default_newworld_rel_syst_Xs.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_rel_syst_det.Clear();    matrix_default_newworld_rel_syst_det.ResizeTo(default_newworld_rows, default_newworld_rows); 

    ///////
    matrix_default_newworld_abs_syst_total.Clear();  matrix_default_newworld_abs_syst_total.ResizeTo(default_newworld_rows, default_newworld_rows); 

    matrix_meas2fitdata_newworld.Clear(); matrix_meas2fitdata_newworld.ResizeTo(1, default_newworld_rows);
    
    ///////
    cout<<endl;
    cout<<"      ---> measurement"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("hdata_obsch_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %2d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      roostr = TString::Format("default_h1d_meas_%d", idx);
      int bins = h1f_temp->GetNbinsX();
      double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1);
      double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_default_h1d_meas_bins[idx] = bins;
      map_default_h1d_meas_xlow[idx] = xlow;
      map_default_h1d_meas_xhgh[idx] = xhgh;      
           
      for(int ibin=1; ibin<=bins+1; ibin++) {
      	double content = h1f_temp->GetBinContent(ibin);
      	vector_default_newworld_meas.push_back( content );      
      }
      delete h1f_temp;    
    }

    if( default_newworld_rows!=(int)(vector_default_newworld_meas.size()) ) { cerr<<" ---> ERROR: default_newworld_rows!=vector_default_newworld_meas"<<endl; exit(1); }
    matrix_default_newworld_meas.Clear(); matrix_default_newworld_meas.ResizeTo(1, default_newworld_rows);
    matrix_default_newworld_meas_POT.Clear(); matrix_default_newworld_meas_POT.ResizeTo(1, default_newworld_rows);
    for(int idx=0; idx<(int)(vector_default_newworld_meas.size()); idx++)  matrix_default_newworld_meas(0, idx) = vector_default_newworld_meas.at(idx);
  
    ///////
    cout<<endl;
    cout<<"      ---> prediction"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("histo_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %2d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      roostr = TString::Format("default_h1d_pred_%d", idx);
      int bins = h1f_temp->GetNbinsX();
      double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1);
      double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_default_h1d_pred_bins[idx] = bins;
      map_default_h1d_pred_xlow[idx] = xlow;
      map_default_h1d_pred_xhgh[idx] = xhgh;      
      
      for(int ibin=1; ibin<=bins+1; ibin++) {
	double content = h1f_temp->GetBinContent(ibin);
	vector_default_oldworld_pred.push_back( content );            
      }
      delete h1f_temp;     
    }
  
    if( default_oldworld_rows!=(int)(vector_default_oldworld_pred.size()) ) { cerr<<" ---> ERROR: default_oldworld_rows!=vector_default_oldworld_pred"<<endl; exit(1); }
    matrix_default_newworld_pred.Clear(); matrix_default_newworld_pred.ResizeTo(1, default_newworld_rows);
    matrix_default_oldworld_pred.Clear(); matrix_default_oldworld_pred.ResizeTo(1, default_oldworld_rows);
    for(int idx=0; idx<(int)(vector_default_oldworld_pred.size()); idx++)  matrix_default_oldworld_pred(0, idx) = vector_default_oldworld_pred.at(idx);
  
    delete roofile_default_cv_file;
  }

  
  {
    for(auto it=map_default_h1d_meas_bins.begin(); it!=map_default_h1d_meas_bins.end(); it++) {
      int idx = it->first;
      int bins = it->second;
      double xlow = map_default_h1d_meas_xlow[idx];
      double xhgh = map_default_h1d_meas_xhgh[idx];
      roostr = TString::Format("default_h1d_meas_%d", idx);
      map_default_h1d_meas[idx] = new TH1D(roostr, roostr, bins, xlow, xhgh);
    }
    
    for(auto it=map_default_h1d_pred_bins.begin(); it!=map_default_h1d_pred_bins.end(); it++) {
      int idx = it->first;
      int bins = it->second;
      double xlow = map_default_h1d_pred_xlow[idx];
      double xhgh = map_default_h1d_pred_xhgh[idx];
      roostr = TString::Format("default_h1d_pred_%d", idx);
      map_default_h1d_pred[idx] = new TH1D(roostr, roostr, bins, xlow, xhgh);
    }
    
    ///////
    
    TFile *roofile_default_cv_file = new TFile(default_cv_file, "read");

    for(auto it=map_default_h1d_meas_bins.begin(); it!=map_default_h1d_meas_bins.end(); it++) {
      int idx = it->first;
      int bins = it->second;
      roostr = TString::Format("hdata_obsch_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_default_h1d_meas[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    for(auto it=map_default_h1d_pred_bins.begin(); it!=map_default_h1d_pred_bins.end(); it++) {
      int idx = it->first;
      int bins = it->second;
      roostr = TString::Format("histo_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_default_h1d_pred[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    delete roofile_default_cv_file;
  }
  
  //////////////////////////////////////
   
  cout<<endl;
  cout<<"      ---> Dirt: additional uncertainty"<<endl;
  {
    TFile *roofile_default_dirtadd_file = new TFile(default_dirtadd_file, "read");
    TMatrixD *matrix_temp = (TMatrixD*)roofile_default_dirtadd_file->Get("cov_mat_add");
    matrix_default_oldworld_abs_syst_addi = (*matrix_temp);
    delete matrix_temp;
    delete roofile_default_dirtadd_file;
  }
  
  //////////////////////////////////////
  
  cout<<endl;
  cout<<"      ---> MCstat"<<endl;
  {
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

  //////////////////////////////////////
  
  cout<<endl;
  cout<<"      ---> flux, geant, and Xs"<<endl;
  {
    TFile *roofile_syst_flux = new TFile(default_fluxXs_dir+"cov_3.root", "read");
    TMatrixD *matrix_temp_flux = (TMatrixD*)roofile_syst_flux->Get("frac_cov_xf_mat_3");
    matrix_default_oldworld_abs_syst_flux = (*matrix_temp_flux);
    matrix_default_oldworld_rel_syst_flux = (*matrix_temp_flux);
    delete matrix_temp_flux;
    delete roofile_syst_flux;
    for(int ibin=0; ibin<default_oldworld_rows; ibin++) {
      for(int jbin=0; jbin<default_oldworld_rows; jbin++) {
	double rel_cov = matrix_default_oldworld_abs_syst_flux(ibin, jbin);
	double cv_i = matrix_default_oldworld_pred(0, ibin);
	double cv_j = matrix_default_oldworld_pred(0, jbin);
	double abs_cov = rel_cov * cv_i * cv_j;
	matrix_default_oldworld_abs_syst_flux(ibin, jbin) = abs_cov;
      }
    }
  }
  
  
  {
    for(int idx=14; idx<=16; idx++) {
      TFile *roofile_syst_geant = new TFile(default_fluxXs_dir+TString::Format("cov_%d.root", idx), "read");
      TMatrixD *matrix_temp_geant = (TMatrixD*)roofile_syst_geant->Get(TString::Format("frac_cov_xf_mat_%d", idx) );
      matrix_default_oldworld_abs_syst_geant += (*matrix_temp_geant);
      matrix_default_oldworld_rel_syst_geant += (*matrix_temp_geant);
      delete matrix_temp_geant;
      delete roofile_syst_geant;
    }
    for(int ibin=0; ibin<default_oldworld_rows; ibin++) {
      for(int jbin=0; jbin<default_oldworld_rows; jbin++) {
	double rel_cov = matrix_default_oldworld_abs_syst_geant(ibin, jbin);
	double cv_i = matrix_default_oldworld_pred(0, ibin);
	double cv_j = matrix_default_oldworld_pred(0, jbin);
	double abs_cov = rel_cov * cv_i * cv_j;
	matrix_default_oldworld_abs_syst_geant(ibin, jbin) = abs_cov;
      }
    }
  }

  
  {
    TFile *roofile_syst_Xs = new TFile(default_fluxXs_dir+"cov_17.root", "read");
    TMatrixD *matrix_temp_Xs = (TMatrixD*)roofile_syst_Xs->Get("frac_cov_xf_mat_17");
    matrix_default_oldworld_abs_syst_Xs = (*matrix_temp_Xs);
    matrix_default_oldworld_rel_syst_Xs = (*matrix_temp_Xs);
    delete matrix_temp_Xs;
    delete roofile_syst_Xs;
    for(int ibin=0; ibin<default_oldworld_rows; ibin++) {
      for(int jbin=0; jbin<default_oldworld_rows; jbin++) {
	double rel_cov = matrix_default_oldworld_abs_syst_Xs(ibin, jbin);
	double cv_i = matrix_default_oldworld_pred(0, ibin);
	double cv_j = matrix_default_oldworld_pred(0, jbin);
	double abs_cov = rel_cov * cv_i * cv_j;
	matrix_default_oldworld_abs_syst_Xs(ibin, jbin) = abs_cov;
      }
    }
  }
    
  //////////////////////////////////////
  
  cout<<endl;
  cout<<"      ---> detector"<<endl;
  {
    map<int, TString>map_detectorfile_str;    
    map_detectorfile_str[1] = default_detector_dir+"cov_LYDown.root";
    map_detectorfile_str[2] = default_detector_dir+"cov_LYRayleigh.root";
    map_detectorfile_str[3] = default_detector_dir+"cov_Recomb2.root";
    map_detectorfile_str[4] = default_detector_dir+"cov_SCE.root";
    //map_detectorfile_str[5] = default_detector_dir+"cov_WMdEdx.root";
    map_detectorfile_str[6] = default_detector_dir+"cov_WMThetaXZ.root";
    map_detectorfile_str[7] = default_detector_dir+"cov_WMThetaYZ.root";
    map_detectorfile_str[8] = default_detector_dir+"cov_WMX.root";
    map_detectorfile_str[9] = default_detector_dir+"cov_WMYZ.root";
    map_detectorfile_str[10]= default_detector_dir+"cov_LYatt.root";

    for(auto it_map=map_detectorfile_str.begin(); it_map!=map_detectorfile_str.end(); it_map++) {
      int idx = it_map->first;
      TFile *roofile_det = new TFile(map_detectorfile_str[idx], "read");
      TMatrixD *matrix_temp_det = (TMatrixD*)roofile_det->Get(TString::Format("frac_cov_det_mat_%d", idx) );
      matrix_default_oldworld_abs_syst_det += (*matrix_temp_det);
      matrix_default_oldworld_rel_syst_det += (*matrix_temp_det);
      delete matrix_temp_det;
      delete roofile_det;
    }

    for(int ibin=0; ibin<default_oldworld_rows; ibin++) {
      for(int jbin=0; jbin<default_oldworld_rows; jbin++) {
	double rel_cov = matrix_default_oldworld_abs_syst_det(ibin, jbin);
	double cv_i = matrix_default_oldworld_pred(0, ibin);
	double cv_j = matrix_default_oldworld_pred(0, jbin);
	double abs_cov = rel_cov * cv_i * cv_j;
	matrix_default_oldworld_abs_syst_det(ibin, jbin) = abs_cov;
      }
    }      
  }
  
  //////////////////////////////////////
    
  cout<<endl;
}
