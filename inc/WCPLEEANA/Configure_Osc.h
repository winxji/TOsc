namespace Configure_Osc
{
  /////////////////////////// default files for spectra and covariance matrixes

  TString default_cv_file      = "./data_inputs/hist_rootfiles_default_noosc/merge.root";
  TString default_mcstat_file  = "./data_inputs/hist_rootfiles_default_noosc/mc_stat/0.log";    
  TString default_fluxXs_dir   = "./data_inputs/hist_rootfiles_default_noosc/XsFlux/";
  TString default_detector_dir = "./data_inputs/hist_rootfiles_default_noosc/DetVar/";

  ///////
  
  bool flag_NuMI_nue2nue   = 1;
  bool flag_NuMI_numu2numu = 0;
  bool flag_NuMI_numu2nue  = 0;  
  bool flag_NuMI_nue2numu  = 0;
  bool flag_NuMI_NC_1minus_nue2sterile  = 0;
  bool flag_NuMI_NC_1minus_numu2sterile = 0;
  
  bool flag_BNB_nue2nue   = 0;
  bool flag_BNB_numu2numu = 0;
  bool flag_BNB_numu2nue  = 0;  
  bool flag_BNB_nue2numu  = 0;
  bool flag_BNB_NC_1minus_nue2sterile  = 0;
  bool flag_BNB_NC_1minus_numu2sterile = 0;
  
}
