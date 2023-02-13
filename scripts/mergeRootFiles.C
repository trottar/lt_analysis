/*
 * Description:
 * ================================================================
 * Time-stamp: "2023-02-13 14:17:51 trottar"
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
#include <TObject.h>

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

/*
int main(int argc, char* argv[]) {

  if (argc != 6) {
    cout << "Error: Incorrect number of arguments provided." << endl;
    cout << "Usage: ./<executable> <RootPath> <InputFileName> <InputTreeName> <OutputFileName> <StringRunNums>" << endl;
    return 0;
  }
  
  TString RootPath = argv[1];
  TString InputFileName = argv[2];
  TString InputTreeName = argv[3];
  TString OutputFileName = argv[4];
  string StringRunNums = argv[5];
  
  mergeRootFiles(RootPath, InputFileName, InputTreeName, OutputFileName, StringRunNums);
  
  return 0;
}
*/
