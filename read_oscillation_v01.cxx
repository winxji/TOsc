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

  osc_test->flag_NuMI_numuCC_from_overlaynumu = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue  = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue      = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu     = Configure_Osc::flag_NuMI_numuCC_from_appnumu;

  osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;

  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;

  osc_test->flag_BNB_numuCC_from_overlaynumu = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue  = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue      = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu     = Configure_Osc::flag_BNB_numuCC_from_appnumu;

  osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  
  /////// set only one time
  
  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
			       Configure_Osc::default_dirtadd_file,
			       Configure_Osc::default_mcstat_file,
			       Configure_Osc::default_fluxXs_dir,
			       Configure_Osc::default_detector_dir);

  osc_test->Set_oscillation_base();
  
  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
  double val_dm2_41 = 7.2;
  double val_sin2_2theta_14 = 0.26;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  val_dm2_41 = 7.2;
  val_sin2_2theta_14 = 0.26;
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  //osc_test->Set_meas2fitdata();
  osc_test->Set_asimov2fitdata();
  
  ///////
  //osc_test->Plot_user();
  //osc_test->Minimization_OscPars_FullCov(8.0, 0.4, 0, 0, "str_flag_fixpar");

  /////////////////////////////////////////////////////////// exclusion

  cout<<endl;
  cout<<" ---> Exclusion processing"<<endl;

  osc_test->Set_oscillation_pars(0, 0, 0, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  osc_test->Set_asimov2noosc();

  ///////
  
  int bins_theta = 20;
  int bins_dm2   = 20;
  
  /////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  /////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
  TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  for(int ibin=1; ibin<=bins_theta; ibin++) {      
    cout<<TString::Format(" ---> processing %4d/%4d", ibin, bins_theta)<<endl;    
    for(int jbin=1; jbin<=bins_dm2; jbin++) {
      
      double xcenter = h2_space->GetXaxis()->GetBinCenter(ibin);
      double ycenter = h2_space->GetYaxis()->GetBinCenter(jbin);
      
      double grid_sin2_2theta_14 = pow( 10, xcenter );
      double grid_dm2_41         = pow( 10, ycenter );

      double pars_4v[4] = {grid_dm2_41, grid_sin2_2theta_14, 0, 0};
      double pars_3v[4] = {0};
      
      ///////      
      double chi2_4v_on_4vAsimov = 0; double chi2_3v_on_4vAsimov = 0; double dchi2_4vAsimov = 0;
      double chi2_4v_on_3vAsimov = 0; double chi2_3v_on_3vAsimov = 0; double dchi2_3vAsimov = 0;
      double chi2_4v_on_data = 0;     double chi2_3v_on_data = 0;     double dchi2_data = 0;
      
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
      dchi2_4vAsimov = chi2_4v_on_4vAsimov - chi2_3v_on_4vAsimov;
      dchi2_3vAsimov = chi2_4v_on_3vAsimov - chi2_3v_on_3vAsimov;
      dchi2_data     = chi2_4v_on_data - chi2_3v_on_data;
      
      double delta_4v = dchi2_4vAsimov;
      double delta_3v = dchi2_3vAsimov;
      double delta_dd = dchi2_data;
      
      double data_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_dd) * 100.;
      double pred_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v) * 100.;
      double pred_CL_1sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-(2*sqrt(fabs(delta_3v))) ) * 100.;
      double pred_CL_1sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+(2*sqrt(fabs(delta_3v))) ) * 100.;
      double pred_CL_2sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-2*(2*sqrt(fabs(delta_3v))) ) * 100.;
      double pred_CL_2sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+2*(2*sqrt(fabs(delta_3v))) ) * 100.;

      roostr = TString::Format("outfile_theta_%03d_dm2_%03d.txt", ibin, jbin);
      ofstream outfile(roostr, ios::out|ios::trunc);
      outfile<<TString::Format("%3d %3d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
			       ibin, jbin,
			       chi2_4v_on_4vAsimov, chi2_3v_on_4vAsimov,
			       chi2_4v_on_3vAsimov, chi2_3v_on_3vAsimov,
			       chi2_4v_on_data, chi2_3v_on_data,
			       data_CL, pred_CL,
			       pred_CL_1sigma_plus, pred_CL_1sigma_minus,
			       pred_CL_2sigma_plus, pred_CL_2sigma_minus
			       )<<endl;
      outfile.close();
      
    }// for(int jbin=1; jbin<=bins_dm2; jbin++)
  }// for(int ibin=1; ibin<=bins_theta; ibin++)
  
  
  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f",
			display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI)<<endl;
  
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
