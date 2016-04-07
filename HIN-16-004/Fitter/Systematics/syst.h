#ifndef syst_h
#define syst_h

#include "../Macros/Utilities/bin.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>

struct syst {
   string name;
   double value;
};

using namespace std;

vector<TString> fileList_syst(const char* token);
map<anabin, syst> readSyst(const char* systfile);
map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name="Total");
map<anabin, syst> readSyst_all(const char* token, bool doPrintTex=false, const char* texName="Systematics/systs.tex");
void printTex(vector< map<anabin, syst> > theSysts, const char* texName="Systematics/systs.tex", bool isLastTotal=false);
map<anabin, vector<syst> > vm2mv(vector< map<anabin,syst> > v);

vector<TString> fileList_syst(const char* token) {
   vector<TString> ans;

   TString basedir("Systematics/csv/");
   TSystemDirectory dir(token,basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".csv") && (TString(token) == "" || fname.Index(token) != kNPOS)) {
            ans.push_back(basedir+fname);
         }
      }
   }

   return ans;
};


map<anabin, syst> readSyst(const char* systfile) {
   map<anabin, syst> ans;

   ifstream file(systfile);
   if (!(file.good())) return ans;

   string systname; getline(file,systname);

   string line;
   double rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;

   while (file.good()) {
      getline(file,line);
      if (line.size()==0) break;
      TString tline(line.c_str());
      TString t; Int_t from = 0, cnt=0;
      while (tline.Tokenize(t, from , ",")) {
         t.Strip(TString::kBoth,' ');
         value = atof(t.Data());
         if (cnt==0) rapmin = atof(t.Data());
         else if (cnt==1) rapmax = value;
         else if (cnt==2) ptmin = value;
         else if (cnt==3) ptmax = value;
         else if (cnt==4) centmin = value;
         else if (cnt==5) centmax = value;
         else if (cnt>6) {
            cout << "Warning, too many fields, I'll take the last one." << endl;
            continue;
         }
         cnt++;
      }
      anabin thebin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);
      syst thesyst; thesyst.value = value; thesyst.name = systname;
      ans[thebin] = thesyst;
   }

   return ans;
};

map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name) {
   map<anabin, syst> ans;

   vector< map<anabin, syst> >::const_iterator it;
   for (it=theSysts.begin(); it!=theSysts.end(); it++) {
      map<anabin, syst>::const_iterator it2;
      for (it2=it->begin(); it2!=it->end(); it2++) {
         anabin thebin = it2->first;
         syst thesyst = it2->second;
         thesyst.name = name;

         // if we already have a syst for this bin, sum quadractically the existing syst and the new syst
         if (ans.find(thebin) != ans.end()) thesyst.value = sqrt(pow(thesyst.value,2) + pow(ans[thebin].value,2));
         ans[thebin] = thesyst;
      }
   }

   return ans;
};

map<anabin, syst> readSyst_all(const char* token, bool doPrintTex, const char* texName) {
   // token should be PP or PbPb

   vector< map<anabin, syst> > systmap_all;
   vector<TString> filelist = fileList_syst(token);

   for (vector<TString>::const_iterator it=filelist.begin(); it!=filelist.end(); it++) {
      cout << "Reading file " << *it << endl;
      map<anabin,syst> systmap = readSyst(it->Data());
      systmap_all.push_back(systmap);
   }

   map<anabin,syst> ans = combineSyst(systmap_all,token);
   systmap_all.push_back(ans);

   if (doPrintTex) printTex(systmap_all, texName, true);

   return ans;
};

void printTex(vector< map<anabin, syst> > theSysts, const char* texName, bool isLastTotal) {
   unsigned int nsyst = theSysts.size();

   ofstream file(texName);
   file << "\\begin{tabular}{|ccc|"; 
   for (unsigned int i=0; i<nsyst; i++) {
      if (i==nsyst-1) file << "|";
      file << "c|";
   }
   file << "}" << endl;
   file << "\\hline" << endl;
   file << "$|y|$ & \\pt & Centrality";
   for (unsigned int i=0; i<nsyst; i++) file << " & " << theSysts[i].begin()->second.name;
   file<< "\\\\" << endl;
   file << "\\hline" << endl;

   anabin oldbin(0,0,0,0,0,0);
   map<anabin, vector<syst> > themap = vm2mv(theSysts);
   map<anabin, vector<syst> >::const_iterator itm;
   for (itm=themap.begin(); itm!=themap.end(); itm++) {
      vector<syst> v = itm->second;
      if (v.size() != nsyst) {
         cout << "Error, not all systematics have the same size. Aborting." << endl;
         file.close();
         return;
      }
      anabin thebin = itm->first;
      if (thebin.rapbin() == oldbin.rapbin()) {
         file << " - ";
      } else {
         if (itm != themap.begin()) file << "\\hline" << endl;
         file.precision(1); file.setf(ios::fixed);
         file << thebin.rapbin().low() << " $< |y| < $ " << thebin.rapbin().high();
      }
      file << " & ";
      if (thebin.ptbin() == oldbin.ptbin()) {
         file << " - ";
      } else {
         file.precision(1); file.setf(ios::fixed);
         file << thebin.ptbin().low() << " $< \\pt < $ " << thebin.ptbin().high();
      }
      file << " & ";
      file.unsetf(ios::fixed);
      file << thebin.centbin().low()/2 << "-" << thebin.centbin().high()/2 << "\\% ";
      file.precision(2);
      file.setf(ios::fixed);

      for (unsigned int i=0; i<nsyst; i++) {
         file << " & " << 100.*v[i].value;
      }
      file << " \\\\" << endl;

      oldbin = thebin;
   }

   file << "\\end{tabular}" << endl;
   file.close();
   cout << "Closed " << texName << endl;
}

map<anabin, vector<syst> > vm2mv(vector< map<anabin,syst> > v) {
   map<anabin, vector<syst> > ans;

   vector<map<anabin,syst> >::const_iterator itv;
   map<anabin,syst>::const_iterator itm;
   for (itv=v.begin(); itv!=v.end(); itv++) {
      for (itm=itv->begin(); itm!=itv->end(); itm++) {
         anabin thebin = itm->first;
         if (ans.find(thebin)==ans.end()) ans[thebin] = vector<syst>();
         ans[thebin].push_back(itm->second);
      }
   }

   return ans;
}

#endif // ifndef syst_h