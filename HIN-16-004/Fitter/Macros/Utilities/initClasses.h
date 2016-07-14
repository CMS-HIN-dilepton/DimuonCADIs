#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLegend.h"

#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussModel.h"
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"

#include "../CMS/tdrstyle.C"
#include "../CMS/CMS_lumi.C"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>


using namespace std;
using namespace RooFit;

typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct EvtPar {
  StartEnd RunNb;
  int      TriggerBit;
} EvtPar;

typedef struct DiMuonPar {
  MinMax ctau, ctauErr, ctauTrue, M, Pt, AbsRap;
  string ctauCut;
} DiMuonPar;

typedef struct SiMuonPar {
  MinMax Pt, Eta;
} SiMuonPar;

typedef struct InputOpt {
  int        oniaMode;
  EvtPar     PbPb, pp;
} InputOpt;

typedef struct KinCuts {
  StartEnd   Centrality;
  SiMuonPar  sMuon;
  DiMuonPar  dMuon;
} KinCuts;

bool isEqualKinCuts(struct KinCuts cutA, struct KinCuts cutB, bool isPbPb) 
{
  bool cond = true;

  if (isPbPb) {
    cond = cond && (cutA.Centrality.Start    == cutB.Centrality.Start);
    cond = cond && (cutA.Centrality.End      == cutB.Centrality.End);
  }
  cond = cond && (cutA.sMuon.Pt.Min        == cutB.sMuon.Pt.Min);
  cond = cond && (cutA.sMuon.Pt.Max        == cutB.sMuon.Pt.Max);
  cond = cond && (cutA.sMuon.Eta.Min       == cutB.sMuon.Eta.Min);
  cond = cond && (cutA.sMuon.Eta.Max       == cutB.sMuon.Eta.Max);

  cond = cond && (cutA.dMuon.ctau.Min      == cutB.dMuon.ctau.Min);
  cond = cond && (cutA.dMuon.ctau.Max      == cutB.dMuon.ctau.Max);
  cond = cond && (cutA.dMuon.ctauErr.Min   == cutB.dMuon.ctauErr.Min);
  cond = cond && (cutA.dMuon.ctauErr.Max   == cutB.dMuon.ctauErr.Max);
  cond = cond && (cutA.dMuon.ctauTrue.Min  == cutB.dMuon.ctauTrue.Min);
  cond = cond && (cutA.dMuon.ctauTrue.Max  == cutB.dMuon.ctauTrue.Max);
  cond = cond && (cutA.dMuon.ctauCut       == cutB.dMuon.ctauCut);
  cond = cond && (cutA.dMuon.M.Min         == cutB.dMuon.M.Min);
  cond = cond && (cutA.dMuon.M.Max         == cutB.dMuon.M.Max);
  cond = cond && (cutA.dMuon.Pt.Min        == cutB.dMuon.Pt.Min);
  cond = cond && (cutA.dMuon.Pt.Max        == cutB.dMuon.Pt.Max);
  cond = cond && (cutA.dMuon.AbsRap.Min    == cutB.dMuon.AbsRap.Min);
  cond = cond && (cutA.dMuon.AbsRap.Max    == cutB.dMuon.AbsRap.Max);

  return cond;
}


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z; };
ParticleMass Mass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188};


enum class MassModel 
{
    InvalidModel =0,
    SingleGaussian=1, 
    DoubleGaussian=2, 
    SingleCrystalBall=3, 
    DoubleCrystalBall=4, 
    GaussianAndCrystalBall=6, 
    Uniform=7, 
    Chebychev1=8, 
    Chebychev2=9, 
    Chebychev3=10, 
    Chebychev4=11,
    Chebychev5=12,
    Chebychev6=13,
    ExpChebychev1=14,
    ExpChebychev2=15,
    ExpChebychev3=16,
    ExpChebychev4=17,
    ExpChebychev5=18,
    ExpChebychev6=19,
    Exponential=20
};
map< string , MassModel > MassModelDictionary = {
  {"InvalidModel",            MassModel::InvalidModel},
  {"SingleGaussian",          MassModel::SingleGaussian},
  {"DoubleGaussian",          MassModel::DoubleGaussian},
  {"SingleCrystalBall",       MassModel::SingleCrystalBall},
  {"DoubleCrystalBall",       MassModel::DoubleCrystalBall},
  {"GaussianAndCrystalBall",  MassModel::GaussianAndCrystalBall},
  {"Uniform",                 MassModel::Uniform},
  {"Chebychev1",              MassModel::Chebychev1},
  {"Chebychev2",              MassModel::Chebychev2},
  {"Chebychev3",              MassModel::Chebychev3},
  {"Chebychev4",              MassModel::Chebychev4},
  {"Chebychev5",              MassModel::Chebychev5},
  {"Chebychev6",              MassModel::Chebychev6},
  {"ExpChebychev1",           MassModel::ExpChebychev1},
  {"ExpChebychev2",           MassModel::ExpChebychev2},
  {"ExpChebychev3",           MassModel::ExpChebychev3},
  {"ExpChebychev4",           MassModel::ExpChebychev4},
  {"ExpChebychev5",           MassModel::ExpChebychev5},
  {"ExpChebychev6",           MassModel::ExpChebychev6},
  {"Exponential",             MassModel::Exponential}
};


enum class CtauModel 
{     
    InvalidModel=0,
    DoubleGaussianResolution=1, 
    SingleGaussianResolution=2,
    TripleDecay=3,
    SingleSidedDecay=4, 
    Delta=5
};
map< string , CtauModel > CtauModelDictionary = {
  {"InvalidModel",             CtauModel::InvalidModel},
  {"DoubleGaussianResolution", CtauModel::DoubleGaussianResolution},
  {"SingleGaussianResolution", CtauModel::SingleGaussianResolution},
  {"TripleDecay",              CtauModel::TripleDecay},
  {"SingleSidedDecay",         CtauModel::SingleSidedDecay},
  {"Delta",                    CtauModel::Delta}
};

typedef struct CtauPNP {
  CtauModel    Prompt, NonPrompt;
} CtauPNP;

typedef struct CtauMassModel {
  MassModel  Mass;
  CtauPNP    Ctau;
} FitModel;

typedef struct CharmModel {
  CtauMassModel  Jpsi, Psi2S, Bkg;
  CtauModel CtauRes;
  CtauModel CtauTrue, CtauTrueRes;
} CharmModel;

typedef struct OniaModel {
  CharmModel  PbPb, PP;
} OniaModel;


#endif // #ifndef initClasses_h
