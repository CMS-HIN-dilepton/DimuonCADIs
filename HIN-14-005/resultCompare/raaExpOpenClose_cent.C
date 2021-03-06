#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"

#include "TAttMarker.h"
#include "TFile.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TSystem.h"


#include "dats/data_cms2016.h"
#include "dats/data_others2016.h"

#include "../CMS_lumi.C"
#include "../tdrstyle.C"
#include "../textPosition.h"
#endif
void raaExpOpenClose_cent(const char* inputDir = "../macro_raa/outRoot", // the place where the input root files, with the histograms are
                          bool bSavePlots      = true,
                          bool bDoSameYHighPt  = true 
                          )
{
   // set the style
  gSystem->mkdir(Form("./figs/png"), kTRUE);
  gSystem->mkdir(Form("./figs/pdf"), kTRUE);
  setTDRStyle();
  // gStyle->SetCanvasPreferGL(1);
  
  // read CMS graphs
  TFile *pfRaaCms_cent   = new TFile(Form("%s/makeRaa_cent.root",inputDir));
  // ##################### HIGH PT ############################
  
  TGraphErrors *pgCms     = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsi");
  TGraphErrors *pgCmsP    = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiP");
  TGraphErrors *pgCmsSyst = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiSyst");
  TBox *lumi = (TBox*)pfRaaCms_cent->Get("lumi_pr_y024_pt6530");
  lumi->SetFillColor(kRed-9);
  lumi->SetFillStyle(1001);
  pgCms->SetFillColorAlpha(kRed-9,0.5);
  pgCmsSyst->SetFillColorAlpha(kRed-9,0.5);

  if(bDoSameYHighPt)
  {
    pgCms     = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsi_pt6530y012");
    pgCmsP    = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiP_pt6530y012");
    pgCmsSyst = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiSyst_pt6530y012");
    lumi      = (TBox*)pfRaaCms_cent->Get("lumi_pr_y012_pt6530");
    lumi->SetFillColor(kAzure-9);
    lumi->SetFillStyle(1001);
    pgCms->SetFillColorAlpha(kAzure-9,0.5);
    pgCmsSyst->SetFillColorAlpha(kAzure-9,0.5);
  }
  pgCms->SetName("pgCms");
  lumi->SetX1(375); lumi->SetX2(400);
  //-------------------------------------------------------------------- 
  // *********** alice points:
  //  TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  TGraphAsymmErrors *pgAlice_dpt612         = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                    npartBins5_alice, raaAlice_cent_dpt612,
                                                                    npartBinsErr5_alice, npartBinsErr5_alice, 
                                                                    raaAliceStat_cent_dpt612,raaAliceStat_cent_dpt612);
  TGraphAsymmErrors *pgAliceSyst_dpt612     = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                    npartBins5_alice, raaAlice_cent_dpt612, 
                                                                    npartSystX5_alice, npartSystX5_alice,
                                                                    raaAliceSystLow_cent_dpt612,raaAliceSystHigh_cent_dpt612);
  TGraphAsymmErrors *pgAliceSystCorr_dpt612 = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                    npartBins5_alice, raaAlice_cent_dpt612,
                                                                    npartSystX5_alice, npartSystX5_alice,
                                                                    raaAliceSystCorrLow_cent_dpt612,raaAliceSystCorrHigh_cent_dpt612);
  
  pgAlice_dpt612->SetName("pgAlice_dpt612");
  pgAlice_dpt612->SetMarkerStyle(25);
  pgAlice_dpt612->SetMarkerSize(1.2);
  // systmatic boxes
  // uncorrelated, open
  pgAliceSyst_dpt612->SetLineColor(1);
  pgAliceSyst_dpt612->SetFillColor(0);
  pgAliceSyst_dpt612->SetFillStyle(0);

  // correlated, filled
  pgAliceSystCorr_dpt612->SetLineColor(1);
  pgAliceSystCorr_dpt612->SetFillColor(kGray);
  // pgAliceSystCorr_dpt816->SetFillStyle(0);

   //---------------------------------------------------------
  TLine *line = new TLine(0.,1,400,1);
  line->SetLineStyle(1);
  line->SetLineWidth(1);

  TCanvas *pc = new TCanvas("pc","pc");

  TF1 *f4 = new TF1("f4","1",0,400);
  //  f4->SetLineWidth(1);
  f4->GetXaxis()->SetTitle("<N_{part}>");
  f4->GetYaxis()->SetTitle("R_{AA}");
  f4->GetYaxis()->SetRangeUser(0.0,1.5);
  f4->GetXaxis()->CenterTitle(kTRUE);

  f4->Draw();
  lumi->Draw();
  pc->Update();

  pgAliceSystCorr_dpt612->Draw("2");
  pgCmsSyst->Draw("2");
  pgAliceSyst_dpt612->Draw("2");

  
  pgCmsP->Draw("P");
  pgCms->Draw("P");
  
  pgAlice_dpt612->Draw("P");

  // additional info
  CMS_lumi(pc,12003000,0);

  TLegend *leg_cent = new TLegend(0.2,0.78,0.7,0.88,NULL,"brNDC");
  leg_cent->SetBorderSize(0);
  leg_cent->SetTextFont(132);
  leg_cent->SetTextSize(0.03);
  leg_cent->SetLineColor(1);
  leg_cent->SetLineStyle(1);
  leg_cent->SetLineWidth(1);
  leg_cent->SetFillColor(19);
  leg_cent->SetFillStyle(0);

  TLegendEntry *entry_cent;
  entry_cent=leg_cent->AddEntry("raab","Hidden charm: prompt J/#psi","");
  entry_cent->SetTextSize(ltxSetTextSize3);
  entry_cent->SetFillStyle(1001);
  entry_cent->SetLineColor(1);
  entry_cent->SetLineStyle(1);
  entry_cent->SetLineWidth(1);

  if(bDoSameYHighPt)
  {
    entry_cent->SetMarkerStyle(20);
    entry_cent->SetMarkerColor(kAzure+7);
    entry_cent->SetMarkerSize(1.2);
    entry_cent=leg_cent->AddEntry("pgCms", "6.5 < p_{T} < 30 GeV/c, |y| < 1.2","p");
  }
  else
  {
    entry_cent->SetMarkerStyle(21);
    entry_cent->SetMarkerColor(kRed+2);
    entry_cent->SetMarkerSize(1.);
    entry_cent=leg_cent->AddEntry("pgCms", "6.5 < p_{T} < 30 GeV/c, |y| < 2.4","p");
  }
  entry_cent->SetTextFont(42);
  entry_cent->SetTextSize(entrySize);


  TLegend *leg_alice_cent = new TLegend(0.2,0.67,0.7,0.77,NULL,"brNDC");
  leg_alice_cent->SetBorderSize(0);
  leg_alice_cent->SetTextFont(132);
  leg_alice_cent->SetLineColor(1);
  leg_alice_cent->SetLineStyle(1);
  leg_alice_cent->SetLineWidth(1);
  leg_alice_cent->SetFillColor(19);
  leg_alice_cent->SetFillStyle(0);
  TLegendEntry *entry_alice_cent=leg_alice_cent->AddEntry("pgAlice_dpt612","Open charm: prompt D (ALICE)","");
  entry_alice_cent->SetFillStyle(1001);
  entry_alice_cent->SetLineColor(1);
  entry_alice_cent->SetLineStyle(1);
  entry_alice_cent->SetLineWidth(1);
  entry_alice_cent->SetTextSize(ltxSetTextSize3);
  entry_alice_cent->SetMarkerStyle(25);
  entry_alice_cent->SetMarkerSize(1);
  entry_alice_cent=leg_alice_cent->AddEntry("pgAlice_dpt612","6 < p_{T} < 12 GeV/c, |y| < 0.5","p");
  entry_alice_cent->SetTextFont(42);
  entry_alice_cent->SetTextSize(entrySize);
  

  leg_cent->Draw();
  leg_alice_cent->Draw();
  line->Draw();
  gPad->RedrawAxis();

  if(bSavePlots)
  {
    pc->SaveAs(Form("figs/pdf/raaExpOpenClose_cent_sameY%d.pdf",bDoSameYHighPt));
    pc->SaveAs(Form("figs/png/raaExpOpenClose_cent_sameY%d.png",bDoSameYHighPt));
  }
  // ##################### LOW PT ############################
  // ##################### LOW PT ############################
  // ##################### LOW PT ############################
  // ##################### LOW PT ############################
  // ##################### LOW PT ############################
  
  TGraphErrors *pgCms_low     = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsi_pt365y1624");
  TGraphErrors *pgCmsP_low    = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiP_pt365y1624");
  TGraphErrors *pgCmsSyst_low = (TGraphErrors *)pfRaaCms_cent->Get("gPrJpsiSyst_pt365y1624");
  pgCmsSyst_low->SetFillColorAlpha(kViolet-9, 0.4);
  pgCms_low->SetName("pgCms_low");

  TBox *lumi_low = (TBox*)pfRaaCms_cent->Get("lumi_pr_y1624_pt365");
  lumi_low->SetFillColor(kViolet-9);
  lumi_low->SetFillStyle(1001);
  

  //-------------------------------------------------------------------- 
  // *********** alice points:
  //  TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  TGraphAsymmErrors *pgAlice_dpt25        = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                  npartBins5_alice, raaAlice_cent_dpt25,
                                                                  npartBinsErr5_alice, npartBinsErr5_alice, 
                                                                  raaAliceStat_cent_dpt25,raaAliceStat_cent_dpt25);
  TGraphAsymmErrors *pgAliceSyst_dpt25     = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                    npartBins5_alice, raaAlice_cent_dpt25, 
                                                                    npartSystX5_alice, npartSystX5_alice,
                                                                    raaAliceSystLow_cent_dpt612,raaAliceSystHigh_cent_dpt25);
  TGraphAsymmErrors *pgAliceSystCorr_dpt25 = new TGraphAsymmErrors(nNpartBins5_alice,
                                                                    npartBins5_alice, raaAlice_cent_dpt25,
                                                                    npartSystX5_alice, npartSystX5_alice,
                                                                    raaAliceSystCorrLow_cent_dpt25,raaAliceSystCorrHigh_cent_dpt25);
  
  pgAlice_dpt25->SetName("pgAlice_dpt25");
  pgAlice_dpt25->SetMarkerStyle(25);
  pgAlice_dpt25->SetMarkerSize(1.2);

  // systmatic boxes
  // uncorrelated, open
  pgAliceSyst_dpt25->SetLineColor(1);
  pgAliceSyst_dpt25->SetFillColor(0);
  pgAliceSyst_dpt25->SetFillStyle(0);

  // correlated, filled
  pgAliceSystCorr_dpt25->SetLineColor(1);
  pgAliceSystCorr_dpt25->SetFillColor(kGray);
 
  //---------------------------------------------------------
  TCanvas *pc_low = new TCanvas("pc_low","pc_low");

  f4->Draw();
  lumi_low->Draw();
  pc_low->Update();

  pgAliceSystCorr_dpt25->Draw("2");
  pgCmsSyst_low->Draw("2");
  pgAliceSyst_dpt25->Draw("2");
  
  pgCmsP_low->Draw("P");
  pgCms_low->Draw("P");
  
  pgAlice_dpt25->Draw("P");

  // additional info
  // additional info
  CMS_lumi(pc_low,12003000,0);

  TLegend *leg2_cent = new TLegend(0.2,0.78,0.7,0.88,NULL,"brNDC");
  leg2_cent->SetBorderSize(0);
  leg2_cent->SetTextFont(132);
  leg2_cent->SetTextSize(0.03);
  leg2_cent->SetLineColor(1);
  leg2_cent->SetLineStyle(1);
  leg2_cent->SetLineWidth(1);
  leg2_cent->SetFillColor(19);
  leg2_cent->SetFillStyle(0);

  TLegendEntry *entry2_cent;
  entry2_cent=leg2_cent->AddEntry("raalow","Hidden charm: prompt J/#psi","");
  entry2_cent->SetTextFont(132);
  entry2_cent->SetTextSize(ltxSetTextSize3);
  entry2_cent=leg2_cent->AddEntry("pgCms_low", "3 < p_{T} < 6.5 GeV/c, 1.6 < |y| < 2.4","p");
  entry2_cent->SetTextFont(42);
  entry2_cent->SetTextSize(entrySize);


  TLegend *leg2_alice_cent = new TLegend(0.2,0.67,0.7,0.77,NULL,"brNDC");
  leg2_alice_cent->SetBorderSize(0);
  leg2_alice_cent->SetTextFont(132);
  leg2_alice_cent->SetTextSize(0.03);
  leg2_alice_cent->SetLineColor(1);
  leg2_alice_cent->SetLineStyle(1);
  leg2_alice_cent->SetLineWidth(1);
  leg2_alice_cent->SetFillColor(19);
  leg2_alice_cent->SetFillStyle(0);
  TLegendEntry *entry2_alice_cent=leg2_alice_cent->AddEntry("pgAlice_dpt25","Open charm: prompt D^{0} (ALICE)","");
  entry2_alice_cent->SetTextFont(132);
  entry2_alice_cent->SetTextSize(ltxSetTextSize3);
  entry2_alice_cent=leg2_alice_cent->AddEntry("pgAlice_dpt25","2 < p_{T} < 5 GeV/c, |y| < 0.5","p");
  entry2_alice_cent->SetTextFont(42);
  entry2_alice_cent->SetTextSize(entrySize);

  leg2_cent->Draw();
  leg2_alice_cent->Draw();
  line->Draw();

  gPad->RedrawAxis();

  if(bSavePlots)
    {   
      pc_low->SaveAs("figs/pdf/raaExpOpenClosePt25_cent.pdf");
      pc_low->SaveAs("figs/png/raaExpOpenClosePt25_cent.png");
    }
  
  
}


