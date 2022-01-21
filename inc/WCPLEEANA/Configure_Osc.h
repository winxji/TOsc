namespace Configure_Osc
{
  /////////////////////////// default files for spectra and covariance matrixes

  TString default_cv_file      = "./data_inputs/hist_rootfiles_default_noosc/merge.root";
  TString default_dirtadd_file = "./data_inputs/hist_rootfiles_default_noosc/merge.root";
  TString default_mcstat_file  = "./data_inputs/hist_rootfiles_default_noosc/mc_stat/0.log";    
  TString default_fluxXs_dir   = "./data_inputs/hist_rootfiles_default_noosc/XsFlux/";
  TString default_detector_dir = "./data_inputs/hist_rootfiles_default_noosc/DetVar/";

  ///////
  bool flag_NuMI_nueCC_from_intnue      = 1;
  bool flag_NuMI_nueCC_from_overlaynumu = 0;
  bool flag_NuMI_nueCC_from_appnue      = 0;
  bool flag_NuMI_nueCC_from_appnumu     = 0;
  bool flag_NuMI_nueCC_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 1/557
  bool flag_NuMI_nueCC_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_NuMI_numuCC_from_overlaynumu = 0;
  bool flag_NuMI_numuCC_from_overlaynue  = 0; 
  bool flag_NuMI_numuCC_from_appnue      = 0;
  bool flag_NuMI_numuCC_from_appnumu     = 0;
  bool flag_NuMI_numuCC_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 693.4/129661
  bool flag_NuMI_numuCC_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_NuMI_CCpi0_from_overlaynumu = 0;
  bool flag_NuMI_CCpi0_from_overlaynue  = 0;// approximation: ignore osc-effect. DocDB-36268 (NuMI, pi0-KE): nueCC/data = 4.5/255
  bool flag_NuMI_CCpi0_from_appnue      = 0;
  bool flag_NuMI_CCpi0_from_appnumu     = 0;// approximation: ignore osc-effect. See flag_NuMI_CCpi0_from_overlaynue
  bool flag_NuMI_CCpi0_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 10.4/7953
  bool flag_NuMI_CCpi0_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_NuMI_NCpi0_from_overlaynumu = 0;
  bool flag_NuMI_NCpi0_from_overlaynue  = 0;// approximation: ignore osc-effect. DocDB-36268 (NuMI, pi0-KE): nueCC/data = 34.8/874
  bool flag_NuMI_NCpi0_from_appnue      = 0;
  bool flag_NuMI_NCpi0_from_appnumu     = 0;// approximation: ignore osc-effect. See flag_NuMI_NCpi0_from_overlaynue
  bool flag_NuMI_NCpi0_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 188.3/5936.0
  bool flag_NuMI_NCpi0_from_dirtnumu    = 0;// approximation: ignore osc-effect. 

  ///////
  bool flag_BNB_nueCC_from_intnue      = 0;
  bool flag_BNB_nueCC_from_overlaynumu = 0;
  bool flag_BNB_nueCC_from_appnue      = 0;
  bool flag_BNB_nueCC_from_appnumu     = 0;
  bool flag_BNB_nueCC_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 1/557
  bool flag_BNB_nueCC_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_BNB_numuCC_from_overlaynumu = 0;
  bool flag_BNB_numuCC_from_overlaynue  = 0; 
  bool flag_BNB_numuCC_from_appnue      = 0;
  bool flag_BNB_numuCC_from_appnumu     = 0;
  bool flag_BNB_numuCC_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 693.4/129661
  bool flag_BNB_numuCC_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_BNB_CCpi0_from_overlaynumu = 0;
  bool flag_BNB_CCpi0_from_overlaynue  = 0;// approximation: ignore osc-effect. LEE PRD paper(pi0-KE): nueCC/data =  18.0/7953
  bool flag_BNB_CCpi0_from_appnue      = 0;
  bool flag_BNB_CCpi0_from_appnumu     = 0;// approximation: ignore osc-effect. flag_BNB_CCpi0_from_overlaynue
  bool flag_BNB_CCpi0_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 10.4/7953
  bool flag_BNB_CCpi0_from_dirtnumu    = 0;// approximation: ignore osc-effect.

  bool flag_BNB_NCpi0_from_overlaynumu = 0;
  bool flag_BNB_NCpi0_from_overlaynue  = 0;// approximation: ignore osc-effect. LEE PRD paper(pi0-KE): nueCC/data = 42.2/5936
  bool flag_BNB_NCpi0_from_appnue      = 0;
  bool flag_BNB_NCpi0_from_appnumu     = 0;// approximation: ignore osc-effect. flag_BNB_NCpi0_from_overlaynue
  bool flag_BNB_NCpi0_from_dirtnue     = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 188.3/5936.0
  bool flag_BNB_NCpi0_from_dirtnumu    = 0;// approximation: ignore osc-effect. 


  
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
