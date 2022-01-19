#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

////////////////////////////////////////////////////////////////////////////////////////////////////

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
  cout<<endl;
  
  vec_ratioPOT->push_back( data_pot/mc_pot );// from the output of ./convert_histo.pl
  vec_vec_eventinfo->push_back( vector_eventinfo );
  
}
 

void TOsc::Set_oscillation_base()
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_oscillation_base"<<endl;

  //////////////////////////////////////

  matrix_oscillation_base_NuMI_nueCC.Clear();  matrix_oscillation_base_NuMI_nueCC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_NuMI_numuCC.Clear(); matrix_oscillation_base_NuMI_numuCC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_NuMI_nueNC.Clear();  matrix_oscillation_base_NuMI_nueNC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_NuMI_numuNC.Clear(); matrix_oscillation_base_NuMI_numuNC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  
  matrix_oscillation_base_BNB_nueCC.Clear();  matrix_oscillation_base_BNB_nueCC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_BNB_numuCC.Clear(); matrix_oscillation_base_BNB_numuCC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_BNB_nueNC.Clear();  matrix_oscillation_base_BNB_nueNC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_oscillation_base_BNB_numuNC.Clear(); matrix_oscillation_base_BNB_numuNC.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  
  //////////////////////////////////////

  if( flag_NuMI_nue2nue ) {
    
    {// intrinsic run1     
      TString strfile_mcPOT = "./data_inputs/hist_rootfiles_default_noosc/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root";
      TString strfile_dataPOT = "./data_inputs/hist_rootfiles_default_noosc/run1_data_bnb.root";
      TString strfile_mc_e2e = "./data_inputs/hist_rootfiles_default_noosc/roofile_obj_NuMI_run1_intrinsic_nue.root";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_matrix_oscillation_base_NuMI_nueCC_scaleFPOT, strfile_mc_e2e, &vector_vector_matrix_oscillation_base_NuMI_nueCC);
    }

    TH1D *h1_nue_FC = (TH1D*)map_default_h1d_pred[1]->Clone("h1_nue_FC"); h1_nue_FC->Reset();// define and clear
    TH1D *h1_nue_PC = (TH1D*)map_default_h1d_pred[1]->Clone("h1_nue_PC"); h1_nue_PC->Reset();// define and clear

    int size = vector_vector_matrix_oscillation_base_NuMI_nueCC.size();
    for(int isize=0; isize<size; isize++) {
      double scaleFPOT = vector_matrix_oscillation_base_NuMI_nueCC_scaleFPOT.at(isize);
      int num_events = vector_vector_matrix_oscillation_base_NuMI_nueCC.at(isize).size();

      for(int ievent=0; ievent<num_events; ievent++) {
	EventInfo info = vector_vector_matrix_oscillation_base_NuMI_nueCC.at(isize).at(ievent);
	
      }
    }
    
    
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

////////////////////////////////////////////////////////////////////////////////////////////////////

void TOsc::Set_default_cv_cov(TString default_cv_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir)
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_default_cv_cov"<<endl;

  cout<<endl;
  cout<<TString::Format("      ---> default_cv_file       %50s", default_cv_file.Data())<<endl;
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

    ///////
    matrix_default_oldworld_abs_syst_addi.Clear();   matrix_default_oldworld_abs_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_mcstat.Clear(); matrix_default_oldworld_abs_syst_mcstat.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_flux.Clear();   matrix_default_oldworld_abs_syst_flux.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_geant.Clear();  matrix_default_oldworld_abs_syst_geant.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_Xs.Clear();     matrix_default_oldworld_abs_syst_Xs.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_default_oldworld_abs_syst_det.Clear();    matrix_default_oldworld_abs_syst_det.ResizeTo(default_oldworld_rows, default_oldworld_rows);

    matrix_default_newworld_abs_syst_addi.Clear();   matrix_default_newworld_abs_syst_addi.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_mcstat.Clear(); matrix_default_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_flux.Clear();   matrix_default_newworld_abs_syst_flux.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_geant.Clear();  matrix_default_newworld_abs_syst_geant.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_Xs.Clear();     matrix_default_newworld_abs_syst_Xs.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_default_newworld_abs_syst_det.Clear();    matrix_default_newworld_abs_syst_det.ResizeTo(default_newworld_rows, default_newworld_rows); 
  
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
    matrix_default_newworld_meas.Clear();
    matrix_default_newworld_meas.ResizeTo(1, default_newworld_rows);
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
    matrix_default_oldworld_pred.Clear();
    matrix_default_oldworld_pred.ResizeTo(1, default_oldworld_rows);
    for(int idx=0; idx<(int)(vector_default_oldworld_pred.size()); idx++)  matrix_default_oldworld_pred(0, idx) = vector_default_oldworld_pred.at(idx);
  
    //////////////////////////////////////
  
    cout<<endl;
    cout<<"      ---> Dirt: additional uncertainty"<<endl;
    {
      TMatrixD *matrix_temp = (TMatrixD*)roofile_default_cv_file->Get("cov_mat_add");
      matrix_default_oldworld_abs_syst_addi = (*matrix_temp);
      delete matrix_temp;
    }
  
    ///////
  
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
