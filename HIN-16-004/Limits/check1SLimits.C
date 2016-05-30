#ifndef check1SLimits_C
#define check1SLimits_C

#include "TString.h"
#include "TH1D.h"
#include "TLine.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"

#include <vector>
#include <map>

#include "../Fitter/Macros/Utilities/bin.h"
#include "../Fitter/Macros/Utilities/resultUtils.h"
#include "limits.h"

using namespace std;

void check1SLimits(
                   const char* workDirName, // workDirName: usual tag where to look for files in Output
                   const char* lFileName="cLimits_683_NominalABCD_Asym_2SPL_woSyst.csv" // file name to save limits results
)
{
  // list of files
  vector<TString> theFiles_PbPb = fileList(workDirName,"PbPb","DATA","../Fitter");
  vector<TString> theFiles_PP = fileList(workDirName,"PP","DATA","../Fitter");
  if ( (theFiles_PbPb.size() < 1) || (theFiles_PP.size() < 1) )
  {
    cout << "#[Error]: No files found in " << workDirName << endl;
    return;
  }
  
  // bin edges
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  
  // histo for 1sigma limits checks
  TH1* hCL = new TH1D("hOneSigmaCLComparison","",theFiles_PbPb.size(),0,theFiles_PbPb.size());
  hCL->GetYaxis()->SetTitle("CL_{1#sigma}/#sigma");
  hCL->GetYaxis()->SetTitleOffset(1.15);
  hCL->SetStats(0);
  hCL->SetDirectory(0);
  hCL->SetMarkerColor(1);
  hCL->SetMarkerStyle(20);
  hCL->SetMarkerSize(1);
  hCL->SetLineColor(1);
  
  TLine* l1 = new TLine(0.,1.,hCL->GetXaxis()->GetXmax(),1.);
  l1->SetLineWidth(3);
  
  hCL->GetListOfFunctions()->Add(l1);
  
  int cnt=1;
  for (vector<TString>::const_iterator it_PbPb=theFiles_PbPb.begin(); it_PbPb!=theFiles_PbPb.end(); it_PbPb++)
  {
    cout << "Checking 1 sigma limits for analysis bin " << cnt << endl;
    cout << "PbPb workspace " << cnt << " / " << theFiles_PbPb.size() << ": " << *it_PbPb << endl;
    
    TFile *f = TFile::Open(*it_PbPb,"READ");
    if (!f)
    {
      cout << "#[Error]: Unable to read file " << *it_PbPb << endl;
      return;
    }
    
    // Retrieve workspace from file
    RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
    if (!ws)
    {
      cout << "#[Error]: Unable to retrieve PbPb workspace" << endl;
      return;
    }
    
    RooRealVar* singleR_PbPb = ws->var("RFrac2Svs1S_PbPb");
    
    anabin thebin = binFromFile(*it_PbPb);
    ptmin = thebin.ptbin().low();
    ptmax = thebin.ptbin().high();
    ymin = thebin.rapbin().low();
    ymax = thebin.rapbin().high();
    centmin = thebin.centbin().low();
    centmax = thebin.centbin().high();
    
    bool foundPPws = false;
    RooRealVar* singleR_PP(0x0);
    for (vector<TString>::const_iterator it_PP=theFiles_PP.begin(); it_PP!=theFiles_PP.end(); it_PP++)
    {
      if ( !foundPPws && isSameBinPPPbPb(*it_PbPb, *it_PP) )
      {
        foundPPws = true;
        
        cout << "PP workspace " << cnt << " / " << theFiles_PP.size() << ": " << *it_PP << endl;
        
        TFile *fPP = TFile::Open(*it_PP,"READ");
        if (!fPP)
        {
          cout << "#[Error]: Unable to read file " << *it_PP << endl;
          return;
        }
        // Retrieve workspace from file
        RooWorkspace* wsPP = (RooWorkspace*) fPP->Get("workspace");
        if (!wsPP)
        {
          cout << "#[Error]: Unable to retrieve PP workspace" << endl;
          return;
        }
        
        singleR_PP = wsPP->var("RFrac2Svs1S_PP");
        
      }
      else continue;
    }
    if ( !foundPPws ) cout << "# [Error]: No PP workspace found for " << *it_PbPb << endl;
    
    RooRealVar* doubleRatio = ratioVar(singleR_PbPb,singleR_PP,1);
    double sigmaDoubleR = doubleRatio->getError();;
    
    // Get limits from file
    map<anabin,limits> maplim = readLimits(Form("csv/%s",lFileName));
    
    bool foundLims = false;
    map<anabin,limits>::const_iterator it;
    for (it=maplim.begin(); it!=maplim.end(); it++)
    {
      anabin thebinLim = it->first;
      if ( (!foundLims) && (ptmin == thebinLim.ptbin().low()) && (ptmax == thebinLim.ptbin().high()) && (ymin == thebinLim.rapbin().low()) && (ymax == thebinLim.rapbin().high()) && (centmin == thebinLim.centbin().low()) && (centmax == thebinLim.centbin().high()) )
      {
        foundLims = true;
        
        limits lim = it->second;
        
        TString binName(Form("Pt[%.1f,%.1f]-Y[%.1f,%.1f]-C[%.1f,%.1f]",ptmin,ptmax,ymin,ymax,centmin,centmax));
        
        double comp = -1.;
        if ( sigmaDoubleR != 0 ) comp = (lim.val.second-lim.val.first)/(2.*sigmaDoubleR);
        hCL->SetBinContent(cnt,comp);
        hCL->GetXaxis()->SetBinLabel(cnt,binName.Data());
      }
      else continue;
      
    }
    
    if ( !foundLims ) cout << "# [Error]: No CL found for " << *it_PbPb << endl;
    
    cnt++;
  } // loop on the files

  TFile* fSave = new TFile("oneSigmaCLComparison.root","RECREATE");
  
  TCanvas* c = new TCanvas("cOneSigmaCLComparison","",90,116,1265,535);
  c->Range(-3.690909,-0.01066472,33.30606,0.01252061);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.1163896);
  c->SetTopMargin(0.03732809);
  c->SetBottomMargin(0.1630648);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  hCL->Draw("p");
  
  c->Write("cOneSigmaCLComparison", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close(); delete fSave;
  
}

#endif // #ifndef check1SLimits_C