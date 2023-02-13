/*
 * Description:
 * ================================================================
 * Time-stamp: "2023-02-13 13:16:04 trottar"
 * ================================================================
 *
 * Author:  Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TObject.h>

using namespace std;

void mergeRootFiles(TString RootPath, TString InputFileName, TString InputTreeName, TString OutputFileName, string StringRunNums) {

  stringstream ss(StringRunNums);

  vector<Int_t> ArrRunNums;

  Int_t runnum;

  while (ss >> runnum) {
    ArrRunNums.push_back(runnum);
  }  
  
  TChain chain(InputTreeName); // change "mytree" to the name of your tree in the root files

  // add root files to the chain, you can add multiple files at once or loop over a list of files
  for (const auto &n : ArrRunNums) {
    cout << n << " ";
    TString str = TString(n);
    chain.Add(RootPath+str+InputFileName+".root"); 
  }

  TFile *outfile = new TFile(RootPath+OutputFileName+".root", "RECREATE");
  TObject *obj;
  
  while ((obj = chain.Get())) {
    obj->Write();
  }
  
  outfile->Close();
  
}
