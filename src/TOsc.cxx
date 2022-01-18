#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

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
    if(h1f_temp==NULL) break;
    cout<<TString::Format("      %2d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

    roostr = TString::Format("default_h1d_meas_%d", idx);
    int bins = h1f_temp->GetNbinsX();
    map_default_h1d_meas[idx] = new TH1D(roostr, roostr, bins, 0, bins);
    for(int ibin=1; ibin<=bins+1; ibin++) {
      double content = h1f_temp->GetBinContent(ibin);
      map_default_h1d_meas[idx]->SetBinContent(ibin, content );
      vector_newworld_meas.push_back( content );      
    }
    delete h1f_temp;    
  }

  if( default_newworld_rows!=(int)(vector_newworld_meas.size()) ) { cerr<<" ---> ERROR: default_newworld_rows!=vector_newworld_meas"<<endl; exit(1); }
  matrix_newworld_meas.Clear();
  matrix_newworld_meas.ResizeTo(1, default_newworld_rows);
  for(int idx=0; idx<(int)(vector_newworld_meas.size()); idx++)  matrix_newworld_meas(0, idx) = vector_newworld_meas.at(idx);
  
  ///////
  cout<<endl;
  cout<<"      ---> prediction"<<endl;
  for(int idx=1; idx<=10000; idx++) {
    roostr = TString::Format("histo_%d", idx);
    TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
    if(h1f_temp==NULL) break;
    cout<<TString::Format("      %2d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

    roostr = TString::Format("default_h1d_pred_%d", idx);
    int bins = h1f_temp->GetNbinsX();
    map_default_h1d_pred[idx] = new TH1D(roostr, roostr, bins, 0, bins);
    for(int ibin=1; ibin<=bins+1; ibin++) {
      double content = h1f_temp->GetBinContent(ibin);
      map_default_h1d_pred[idx]->SetBinContent(ibin, content );
      vector_oldworld_pred.push_back( content );      
    }
    delete h1f_temp;     
  }

  if( default_oldworld_rows!=(int)(vector_oldworld_pred.size()) ) { cerr<<" ---> ERROR: default_oldworld_rows!=vector_oldworld_pred"<<endl; exit(1); }
  matrix_oldworld_pred.Clear();
  matrix_oldworld_pred.ResizeTo(1, default_oldworld_rows);
  for(int idx=0; idx<(int)(vector_oldworld_pred.size()); idx++)  matrix_oldworld_pred(0, idx) = vector_oldworld_pred.at(idx);
  
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
	double cv_i = matrix_oldworld_pred(0, ibin);
	double cv_j = matrix_oldworld_pred(0, jbin);
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
	double cv_i = matrix_oldworld_pred(0, ibin);
	double cv_j = matrix_oldworld_pred(0, jbin);
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
	double cv_i = matrix_oldworld_pred(0, ibin);
	double cv_j = matrix_oldworld_pred(0, jbin);
	double abs_cov = rel_cov * cv_i * cv_j;
	matrix_default_oldworld_abs_syst_Xs(ibin, jbin) = abs_cov;
      }
    }
  }
    
  //////////////////////////////////////
  
  cout<<endl;
  cout<<"      ---> detector"<<endl;
  {

  }
  
  //////////////////////////////////////
  
  cout<<endl;
}
