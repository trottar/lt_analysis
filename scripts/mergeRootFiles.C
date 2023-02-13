/*
 * Description:
 * ================================================================
 * Time-stamp: "2023-02-13 14:39:39 trottar"
 * ================================================================
 *
 * Author:  Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include "mergeRootFiles.h"

using namespace std;

int mergeRootFiles(TString RootPath, TString InputFileName, TString InputTreeName, TString OutputFileName, string StringRunNums) {

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
    TString filepath = RootPath + str + InputFileName + ".root";
    TFile *tempfile = TFile::Open(filepath);
    if (tempfile == nullptr || !tempfile->IsOpen() || tempfile->TestBit(TFile::kRecovered)) {
      cout << "File " << filepath << " not found or not opened or corrupted. Skipping this file." << endl;
      continue;
    }
    chain.Add(filepath);
  }

  TFile *outfile = new TFile(RootPath + OutputFileName + ".root", "RECREATE");
  if (!outfile->IsOpen()) {
    cout << "Output file " << outfile->GetName() << " cannot be opened. Exiting the function." << endl;
    return 0;
  }

  chain.Merge(outfile->GetName());

  outfile->Close();

  return 0;
}
