#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

//#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//milliseconds, minutes

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A story ..."<<endl<<endl;

  int ifile = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;
  
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pbnb")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_BNB ) ) { cerr<<" ---> Error scaleF_POT_BNB !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pnumi")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_NuMI ) ) { cerr<<" ---> Error scaleF_POT_NuMI !"<<endl; exit(1); }
    }    
    if( strcmp(argv[i],"-d")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>display ) ) { cerr<<" ---> Error display !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it24")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it24 ) ) { cerr<<" ---> Error it24 !"<<endl; exit(1); }
    }    
  }

  ///////////////////////////////////////////////////////////
  
  if( !display ) {
    gROOT->SetBatch( 1 );
  }
  
  TApplication theApp("theApp",&argc,argv);
  
  /////////////////////////////////////////////////////////// Draw style

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

  ///////////////////////////////////////////////////////////

  TOsc *osc_test = new TOsc();

  ///////

  osc_test->scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->scaleF_POT_NuMI = scaleF_POT_NuMI;
     
  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
  ///////

  osc_test->flag_NuMI_nueCC_from_intnue       = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue       = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu      = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC= Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
  osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC= Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
  osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
  /////// set only one time
  
  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
			       Configure_Osc::default_dirtadd_file,
			       Configure_Osc::default_mcstat_file,
			       Configure_Osc::default_fluxXs_dir,
			       Configure_Osc::default_detector_dir);
  
  osc_test->Set_oscillation_base();
  
  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
  double val_dm2_41         = 7.3;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  val_dm2_41         = 7.3;
  val_sin2_2theta_14 = 0.1;
  val_sin2_theta_24  = 0.01;  
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  //osc_test->Set_meas2fitdata();
  osc_test->Set_asimov2fitdata();
  
  ///////
  osc_test->Plot_user();
  //osc_test->Minimization_OscPars_FullCov(7.3, 0.1, 0.0048, 0, "str_flag_fixpar");

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if( 0 ) {
    TH1D *h1d_nue = new TH1D("h1d_nue", "h1d_nue", 25, 0, 2500);
    for(int idx=1; idx<=25; idx++) {
      h1d_nue->SetBinContent(idx, osc_test->matrix_eff_newworld_pred(0, 26*7 + idx-1));
    }
    h1d_nue->SaveAs("file_numi_nue.root");

    
    cout<<" ---> test "<<h1d_nue->Integral()<<endl;    
    //h1d_nue->Draw("hist");
    
  }

  /////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////// 
  
  if( 0 ) {
    cout<<endl<<" ---> Asimov scan"<<endl<<endl;
    
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.01;

    int    min_status             = 10;
    int    flag_negative          = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;
    double min_dm2_41_err         = 0;
    double min_sin2_2theta_14_err = 0;
    double min_sin2_theta_24_err  = 0;
    double min_sin2_theta_34_err  = 0;
    double true_dm2_41_val         = val_dm2_41;
    double true_sin2_2theta_14_val = val_sin2_2theta_14;
    double true_sin2_theta_24_val  = val_sin2_theta_24;
    double true_sin2_theta_34_val  = 0;
    
    vector<double> vec_asimov_data;
    vector<double> vec_pseudo_data;
    vector<double> vec_bestfit_pred;
    vector<double> vec_syst_rand;

    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready   
    for(int idx=0; idx<osc_test->matrix_eff_newworld_pred.GetNcols(); idx++) {
      vec_asimov_data.push_back( osc_test->matrix_eff_newworld_pred(0,idx) );
    }
    

    roostr = TString::Format("sub_fit_%06d.root", ifile);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "flag_negative",          &flag_negative,          "flag_negative/I" );
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );
 
    double dm2_low  = 0;
    double dm2_hgh  = 16;
    double dm2_step = 0.1;
    int num_dm2 = (dm2_hgh-dm2_low)/dm2_step;
    //num_dm2 = 2;
    
    double t14_low  = 0;
    double t14_hgh  = 0.5;
    double t14_step = 0.02;
    int num_t14 = (t14_hgh-t14_low)/t14_step;
    //num_t14 = 2;
    
    double t24_low  = 0;
    double t24_hgh  = 0.04;
    double t24_step = 0.002;
    int num_t24 = (t24_hgh-t24_low)/t24_step;
    //num_t24 = 2;
    
    for(int itoy=1; itoy<=1; itoy++) {

      cout<<endl<<" ---> proceesing toy "<<itoy<<endl;
      
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_asimov2fitdata();
      //osc_test->Set_toy_variations( 1 );
      //osc_test->Set_toy2fitdata( 1 );

      flag_negative = 0;
      
      double obj_dm2  = 0;
      double obj_t14  = 0;
      double obj_t24  = 0;
      double obj_chi2 = 1e6;
      
      double pars_4v[4] = {0};

      
      for(int idm2=ifile; idm2<=ifile; idm2++ ) {
	double dm2_val = dm2_low + (ifile-1)*dm2_step;
	cout<<" ---> processing "<<idm2<<endl;
	
	for(int it14=0; it14<num_t14; it14++) {
	  cout<<TString::Format("      idm2 %3d, it14 %3d", idm2, it14)<<endl;
	  for(int it24=0; it24<num_t24; it24++) {
	    //double dm2_val = dm2_low + idm2*dm2_step;
	    double t14_val = t14_low + it14*t14_step;
	    double t24_val = t24_low + it24*t24_step;

	    pars_4v[0] = dm2_val;
	    pars_4v[1] = t14_val;
	    pars_4v[2] = t24_val;
	    
	    double chi2_val = osc_test->FCN( pars_4v );// 1000 times ~ 3 min
	    
	    min_status            = 0;
	    min_chi2              = chi2_val;
	    min_dm2_41_val        = dm2_val;
	    min_sin2_2theta_14_val= t14_val;
	    min_sin2_theta_24_val = t24_val;
	    min_sin2_theta_34_val = 0;

	    cout<<" ---> "<<t14_val<<"\t"<<t24_val<<"\t   "<<chi2_val<<endl;
	    
	    tree->Fill();	   
	  }
	}
      }

    }

    tree->Write();
    subroofile->Close();
  }
  

  
  /////////////////////////////////////////////////////
  
  if( 0 ) {


    // 7.49894 0.0926119 0.00501187
    // 8.42076e-05 // new: no fix
    // 0.000254071 // fixed 680
    // 0.000334495 7.49894 0.0926119 0.00501187 // old: no fix
    
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.005;

    val_dm2_41         = 7.49894;
    val_sin2_2theta_14 = 0.0926119;
    val_sin2_theta_24  = 0.00501187;

    if( 1 ) {      
      double pars_4v[4] = {val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0};
      double pars_3v[4] = {0};

      //////
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_asimov2fitdata();
      
      double chi2_4v_4v = osc_test->FCN( pars_4v );
      double chi2_3v_4v = osc_test->FCN( pars_3v );

      /////
      osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_asimov2fitdata();
      
      double chi2_4v_3v = osc_test->FCN( pars_4v );
      double chi2_3v_3v = osc_test->FCN( pars_3v );

      //////
      cout<<endl;
      cout<<" chi2_4v_4v "<<chi2_4v_4v<<endl;
      cout<<" chi2_3v_4v "<<chi2_3v_4v<<endl;
      cout<<" chi2_4v_3v "<<chi2_4v_3v<<endl;
      cout<<" chi2_3v_3v "<<chi2_3v_3v<<endl;
      cout<<" dchi2_4v "<<chi2_4v_4v - chi2_3v_4v<<endl;
      cout<<" dchi2_3v "<<chi2_4v_3v - chi2_3v_3v<<endl;      
      cout<<endl;

      cout<<osc_test->func_CLs(chi2_4v_4v - chi2_3v_4v, chi2_4v_3v - chi2_3v_3v, chi2_4v_3v - chi2_3v_3v)<<endl;
      
      return 0;
    }
  }
  

  if( 0 ) {
    cout<<endl<<" ---> noosc fracCOV"<<endl<<endl;// matrix_eff_newworld_pred


    int rows = osc_test->matrix_eff_newworld_abs_syst_total.GetNrows();
    
    TMatrixD matrix_osc_fracCOV(rows, rows);

    for(int ibin=1; ibin<=rows; ibin++) {
      for(int jbin=1; jbin<=rows; jbin++) {
	double cv_i = osc_test->matrix_eff_newworld_pred(0, ibin-1);
	double cv_j = osc_test->matrix_eff_newworld_pred(0, jbin-1);
	double cov_ij = osc_test->matrix_eff_newworld_abs_syst_total(ibin-1, jbin-1);

	double fracCOV = 0;
	if( cv_i!=0 && cv_j!=0 ) {
	  fracCOV = cov_ij/cv_i/cv_j;
	}

	matrix_osc_fracCOV(ibin-1, jbin-1) = fracCOV;
      }
      //cout<<ibin<<"\t"<<sqrt(matrix_osc_fracCOV(ibin-1, ibin-1))<<endl;
    }

    TFile *file = new TFile("file_osc_fracCOV.root", "recreate");
    matrix_osc_fracCOV.Write("matrix_osc_fracCOV");
    //file->Close();


    ///////////////
    
    TMatrixD matrix_user_trans(rows, 14);
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*0 + ibin-1, 0) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*1 + ibin-1, 1) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*2 + ibin-1, 2) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*3 + ibin-1, 3) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*4 + ibin-1, 4) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*5 + ibin-1, 5) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*6 + ibin-1, 6) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*7 + ibin-1, 7) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*8 + ibin-1, 8) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*9 + ibin-1, 9) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*10 + ibin-1, 10) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*11 + ibin-1, 11) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*12 + ibin-1, 12) = 1;
    for(int ibin=1; ibin<=26; ibin++) matrix_user_trans(26*13 + ibin-1, 13) = 1;
    
    TMatrixD matrix_user_trans_T = matrix_user_trans.T(); matrix_user_trans.T(); 

    TMatrixD matrix_rate_COV = matrix_user_trans_T * (osc_test->matrix_eff_newworld_abs_syst_total) * matrix_user_trans;
    TMatrixD matrix_rate_pred = (osc_test->matrix_eff_newworld_pred) * matrix_user_trans;


    int ccrr = 14;
    TMatrixD matrix_osc_fracCOV_rate(ccrr, ccrr);
    for(int ibin=1; ibin<=ccrr; ibin++) {
      for(int jbin=1; jbin<=ccrr; jbin++) {
	double cv_i = matrix_rate_pred(0, ibin-1);
	double cv_j = matrix_rate_pred(0, jbin-1);
	double cov_ij = matrix_rate_COV(ibin-1, jbin-1);

	double fracCOV = 0;
	if( cv_i!=0 && cv_j!=0 ) {
	  fracCOV = cov_ij/cv_i/cv_j;
	}

	matrix_osc_fracCOV_rate(ibin-1, jbin-1) = fracCOV;
      }
      cout<<ibin<<"\t"<<sqrt(matrix_osc_fracCOV_rate(ibin-1, ibin-1))<<endl;
    }


    matrix_osc_fracCOV_rate.Write("matrix_osc_fracCOV_rate");
    file->Close();
    
  }
  
  ///////
  
  if( 0 ) {
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0;
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();
    osc_test->Set_asimov2fitdata();

    double pars_4v[4] = {7.3, 0.26, 0, 0};
    double chi2_4v = osc_test->FCN( pars_4v );
    cout<<endl<<" sensitivity ---> "<<chi2_4v<<endl<<endl;
  }
  
  ///////  
  if( 0 ) {
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.36;
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    
    roostr = TString::Format("sub_fit_%06d.root", ifile);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    int    min_status             = 10;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;
    double min_dm2_41_err         = 0;
    double min_sin2_2theta_14_err = 0;
    double min_sin2_theta_24_err  = 0;
    double min_sin2_theta_34_err  = 0;
    double chi2_3v                = 0;
    double chi2_4v                = 0;

    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );
    tree->Branch( "min_dm2_41_err",         &min_dm2_41_err,         "min_dm2_41_err/D" );
    tree->Branch( "min_sin2_2theta_14_err", &min_sin2_2theta_14_err, "min_sin2_2theta_14_err/D" );
    tree->Branch( "min_sin2_theta_24_err",  &min_sin2_theta_24_err,  "min_sin2_theta_24_err/D" );
    tree->Branch( "min_sin2_theta_34_err",  &min_sin2_theta_34_err,  "min_sin2_theta_34_err/D" );
    tree->Branch( "chi2_3v",                &chi2_3v,                "chi2_3v/D" );
    tree->Branch( "chi2_4v",                &chi2_4v,                "chi2_4v/D" );

    int ntoys = 5;
    osc_test->Set_toy_variations(ntoys);
    for(int idx=1; idx<=ntoys; idx++) {      
      osc_test->Set_toy2fitdata(idx);
      osc_test->Minimization_OscPars_FullCov(6.0, 0.2, 0, 0, "str_flag_fixpar");
      
      min_status            = osc_test->minimization_status;
      min_chi2              = osc_test->minimization_chi2;
      min_dm2_41_val        = osc_test->minimization_dm2_41_val;
      min_sin2_2theta_14_val= osc_test->minimization_sin2_2theta_14_val;
      min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
      min_sin2_theta_34_val = osc_test->minimization_sin2_theta_34_val;
      min_dm2_41_err        = osc_test->minimization_dm2_41_err;
      min_sin2_2theta_14_err= osc_test->minimization_sin2_2theta_14_err;
      min_sin2_theta_24_err = osc_test->minimization_sin2_theta_24_err;
      min_sin2_theta_34_err = osc_test->minimization_sin2_theta_34_err;

      double pars_3v[4] = {0};
      chi2_3v = osc_test->FCN( pars_3v );
      
      double pars_4v[4] = {val_dm2_41, val_sin2_2theta_14, 0, 0};
      chi2_4v = osc_test->FCN( pars_4v );
      
      tree->Fill();
    }
 
    tree->Write();
    subroofile->Close();
  }
    
  
  /////////////////////////////////////////////////////////// exclusion

  if( 0 ) {
    
    cout<<endl;
    cout<<" ---> Exclusion processing"<<endl;

    osc_test->Set_oscillation_pars(0, 0, 0, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    osc_test->Set_asimov2noosc();

    ///////
  
    int bins_theta = 60;
    int bins_dm2   = 60;
  
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);

    int low_theta = 1;
    int hgh_theta = bins_theta;
    int low_dm2   = 1;
    int hgh_dm2   = bins_dm2;
    int low_t24   = 1;
    int hgh_t24   = bins_theta;
    
    if( it14!=0 ) {
      low_theta = it14; hgh_theta = it14;
    }
    if( idm2!=0 ) {
      low_dm2 = idm2; hgh_dm2 = idm2;
    }
    if( it24!=0) {
      low_t24 = it24; hgh_t24 = it24;
    }
    
    //double user_sin2_theta_24 = 0;// sscan
	
    for(int ibin=low_theta; ibin<=hgh_theta; ibin++) {
      for(int jbin=low_dm2; jbin<=hgh_dm2; jbin++) {
	
	roostr = TString::Format("outfile_theta_%03d_dm2_%03d_t24.txt", ibin, jbin);
	ofstream outfile(roostr, ios::out|ios::trunc);
	
	for(int kbin=low_t24; kbin<=hgh_t24; kbin++) {
	  cout<<TString::Format(" ---> processing theta,dm2,t24: %3d %3d %3d", ibin, jbin, kbin)<<endl;    
	  
	  double xcenter = h2_space->GetXaxis()->GetBinCenter(ibin);// angle 14
	  double ycenter = h2_space->GetYaxis()->GetBinCenter(jbin);// dm2
	  double kcenter = h2_space->GetXaxis()->GetBinCenter(kbin);// angle 24
	  
	  double grid_sin2_2theta_14 = pow( 10, xcenter );
	  double grid_dm2_41         = pow( 10, ycenter );
	  double grid_sin2_theta_24  = pow( 10, kcenter );

	  //grid_sin2_2theta_14 = 0.36;
	  //grid_dm2_41         = 7.3;

	  double chi2_4v_on_4vAsimov(0), chi2_3v_on_4vAsimov(0);
	  double chi2_4v_on_3vAsimov(0), chi2_3v_on_3vAsimov(0);
	  double chi2_4v_on_data(0), chi2_3v_on_data(0);

	  double pars_4v[4] = {grid_dm2_41, grid_sin2_2theta_14, grid_sin2_theta_24, 0};
	  double pars_3v[4] = {0};
      
	  /////// 4v Asimov      
	  osc_test->Set_oscillation_pars(pars_4v[0], pars_4v[1], pars_4v[2], pars_4v[3]);
	  osc_test->Apply_oscillation();
	  osc_test->Set_apply_POT();
	  osc_test->Set_asimov2fitdata();
	  chi2_3v_on_4vAsimov = osc_test->FCN( pars_3v );

	  /////// 3v Asimov
	  osc_test->Set_noosc2fitdata();
	  chi2_4v_on_3vAsimov = osc_test->FCN( pars_4v );

	  /////// data
	  osc_test->Set_meas2fitdata();
	  chi2_4v_on_data = osc_test->FCN( pars_4v );
	  chi2_3v_on_data = osc_test->FCN( pars_3v );

	  ///////
	  double dchi2_4vAsimov = chi2_4v_on_4vAsimov - chi2_3v_on_4vAsimov;
	  double dchi2_3vAsimov = chi2_4v_on_3vAsimov - chi2_3v_on_3vAsimov;
	  double dchi2_data     = chi2_4v_on_data - chi2_3v_on_data;
      
	  double delta_4v = dchi2_4vAsimov;
	  double delta_3v = dchi2_3vAsimov;
	  double delta_dd = dchi2_data;
      
	  double data_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_dd) * 100.;
	  double pred_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v) * 100.;
	  double pred_CL_1sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_1sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_2sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-2*(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_2sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+2*(2*sqrt(fabs(delta_3v))) ) * 100.;

	  outfile<<TString::Format("%3d %3d %3d %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
				   ibin, jbin, kbin,
				   chi2_4v_on_4vAsimov, chi2_3v_on_4vAsimov,
				   chi2_4v_on_3vAsimov, chi2_3v_on_3vAsimov,
				   chi2_4v_on_data, chi2_3v_on_data,
				   data_CL, pred_CL,
				   pred_CL_1sigma_plus, pred_CL_1sigma_minus,
				   pred_CL_2sigma_plus, pred_CL_2sigma_minus
				   )<<endl;	  

	  cout<<TString::Format(" ---> pred p-value %10.8f", (100-pred_CL)/100)<<endl;

	}// for(int kbin=low_theta; kbin<=hgh_theta; kbin++)

	outfile.close();
	
      }// for(int jbin=1; jbin<=bins_dm2; jbin++)
    }// for(int ibin=1; ibin<=bins_theta; ibin++)
   
      
  }


  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, theta14(-it14) %d, dm2(-idm2) %d, theta24(-it24) %d",
			display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI, it14, idm2, it24)<<endl;
  
  cout<<endl;
  cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
  cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
  cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
  cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
  cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
  cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

  cout<<endl;
  cout<<" ---> Finished sucessfully"<<endl;
  
  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }
  
  return 0;
}
