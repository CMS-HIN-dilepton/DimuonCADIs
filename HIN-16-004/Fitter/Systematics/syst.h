#ifndef syst_h
#define syst_h

#include "../Macros/Utilities/bin.h"
#include "../Macros/Utilities/resultUtils.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCollection.h"
#include "TList.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#include "RooWorkspace.h"
#include "RooRealVar.h"

struct syst {
   string name;
   double value;
   double value_dR;
   double value_dR_rel;
};

using namespace std;

vector<TString> fileList_syst(const char* token, const char* prependPath="");
map<anabin, syst> readSyst(const char* systfile, const char* workDirName="", const char* path2Fitter="", int mode=0, const char* workDirNameFail=""); // mode=0 (pass), 1 (prompt), 2 (nonprompt)
map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name="Total", bool isdRSyst=false);
map<anabin, syst> readSyst_all(const char* prependPath="", const char* workDirName="", bool doPrintTex=false, const char* texName="Systematics/systs.tex", int mode=0, const char* workDirNameFail="");
map<anabin, syst> readSyst_all_pass(const char* token, const char* prependPath="", const char* workDirName="");
map<anabin, syst> readSyst_all_prompt(const char* token, const char* prependPath="", const char* workDirName="", const char* workDirNameFail="");
map<anabin, syst> readSyst_all_nonprompt(const char* token, const char* prependPath="", const char* workDirName="", const char* workDirNameFail="");
void printTex(vector< map<anabin, syst> > theSysts, const char* texName="Systematics/systs.tex", bool isdRSyst=false);
map<anabin, vector<syst> > vm2mv(vector< map<anabin,syst> > v);
RooWorkspace* getWorkspaceFromBin(anabin thebin, const char* workDirName, const char* token="", const char* path2Fitter="");

#endif // ifndef syst_h
