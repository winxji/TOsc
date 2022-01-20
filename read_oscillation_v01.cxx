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
  
  cout<<endl<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f",
			      display, ifile, scaleF_POT_BNB, scaleF_POT_NuMI)<<endl<<endl;

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
  
  osc_test->flag_NuMI_nue2nue   = Configure_Osc::flag_NuMI_nue2nue;
  osc_test->flag_NuMI_numu2numu = Configure_Osc::flag_NuMI_numu2numu;
  osc_test->flag_NuMI_numu2nue  = Configure_Osc::flag_NuMI_numu2nue;  
  osc_test->flag_NuMI_nue2numu  = Configure_Osc::flag_NuMI_nue2numu;
  osc_test->flag_NuMI_NC_1minus_nue2sterile  = Configure_Osc::flag_NuMI_NC_1minus_nue2sterile;
  osc_test->flag_NuMI_NC_1minus_numu2sterile = Configure_Osc::flag_NuMI_NC_1minus_numu2sterile;
  
  osc_test->flag_BNB_nue2nue   = Configure_Osc::flag_BNB_nue2nue;
  osc_test->flag_BNB_numu2numu = Configure_Osc::flag_BNB_numu2numu;
  osc_test->flag_BNB_numu2nue  = Configure_Osc::flag_BNB_numu2nue;  
  osc_test->flag_BNB_nue2numu  = Configure_Osc::flag_BNB_nue2numu;
  osc_test->flag_BNB_NC_1minus_nue2sterile  = Configure_Osc::flag_BNB_NC_1minus_nue2sterile;
  osc_test->flag_BNB_NC_1minus_numu2sterile = Configure_Osc::flag_BNB_NC_1minus_numu2sterile;

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
  osc_test->Set_meas2fitdata();

  ///
  osc_test->Minimization_OscPars_FullCov(8.0, 0.4, 0, 0, "str_flag_fixpar");
    
  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ----------------------------> check information"<<endl;
  cout<<endl;
  cout<<" osc_test->flag_NuMI_nue2nue   " <<osc_test->flag_NuMI_nue2nue<<endl;
  cout<<" osc_test->flag_NuMI_numu2numu " <<osc_test->flag_NuMI_numu2numu<<endl;
  cout<<" osc_test->flag_NuMI_numu2nue  " <<osc_test->flag_NuMI_numu2nue<<endl;  
  cout<<" osc_test->flag_NuMI_nue2numu  " <<osc_test->flag_NuMI_nue2numu<<endl;
  cout<<" osc_test->flag_NuMI_NC_1minus_nue2sterile  " <<osc_test->flag_NuMI_NC_1minus_nue2sterile<<endl;
  cout<<" osc_test->flag_NuMI_NC_1minus_numu2sterile " <<osc_test->flag_NuMI_NC_1minus_numu2sterile<<endl;
  cout<<endl;
  cout<<" osc_test->flag_BNB_nue2nue   " <<osc_test->flag_BNB_nue2nue<<endl;
  cout<<" osc_test->flag_BNB_numu2numu " <<osc_test->flag_BNB_numu2numu<<endl;
  cout<<" osc_test->flag_BNB_numu2nue  " <<osc_test->flag_BNB_numu2nue<<endl;  
  cout<<" osc_test->flag_BNB_nue2numu  " <<osc_test->flag_BNB_nue2numu<<endl;
  cout<<" osc_test->flag_BNB_NC_1minus_nue2sterile  " <<osc_test->flag_BNB_NC_1minus_nue2sterile<<endl;
  cout<<" osc_test->flag_BNB_NC_1minus_numu2sterile " <<osc_test->flag_BNB_NC_1minus_numu2sterile<<endl;
  cout<<endl;
  
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f",
			display, ifile, scaleF_POT_BNB, scaleF_POT_NuMI)<<endl;
  cout<<endl;

  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }
  
  return 0;
}
