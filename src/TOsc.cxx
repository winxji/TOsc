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
  matrix_transform.Clear();
  matrix_transform.ResizeTo(default_oldworld_rows, default_newworld_rows);
  matrix_transform += (*temp_mat_collapse);
  delete temp_mat_collapse;

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
  
  ///////
  cout<<endl;
  delete roofile_default_cv_file;

}
