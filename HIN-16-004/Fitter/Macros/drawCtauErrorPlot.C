#ifndef drawCtauErrorPlot_C
#define drawCtauErrorPlot_C

#include "Utilities/initClasses.h"

void setCtauErrorRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale);

void drawCtauErrorPlot(RooWorkspace& myws,   // Local workspace
                       string outputDir,     // Output directory
                       struct InputOpt opt,  // Variable with run information (kept for legacy purpose)
                       struct KinCuts cut,   // Variable with current kinematic cuts
                       // Select the type of datasets to fit
                       string DSTAG,         // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                       bool isPbp,          // Define if it is Pbp (True) or PP (False)
                       // Select the type of object to fit
                       bool incJpsi,         // Includes Jpsi model
                       bool incPsi2S,        // Includes Psi(2S) model
                       bool incBkg,          // Includes Background model        
                       bool incPrompt,       // Includes Prompt ctau model       
                       bool incNonPrompt,    // Includes Non-Prompt ctau model
                       // Select the fitting options
                       bool cutCtau,         // Apply prompt ctau cuts
                       bool plotPureSMC,     // Flag to indicate if we want to fit pure signal MC
                       // Select the drawing options
                       bool setLogScale,     // Draw plot with log scale
                       bool incSS,           // Include Same Sign data
                       int  nBins            // Number of bins used for plotting
                       ) 
{
  

  RooMsgService::instance().getStream(0).removeTopic(Caching);  
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  TH1::AddDirectory(kFALSE);

  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbp?"Pbp":"PP"));
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbp?"Pbp":"PP"));
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), (isPbp?"Pbp":"PP"));

  // Create the main plot of the fit
  RooPlot*   frame     = myws.var("ctauErr")->frame(Bins(nBins), Range(cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max));
 
  if (incJpsi) {
    string pdfName = ""; string dsName = "";
    if (incPrompt)  { pdfName = Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"));      dsName = "binSubtractedSIG"; }
    myws.data(dsName.c_str())->plotOn(frame, Name("JPSIDATA"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlue+3), MarkerSize(1.0));
    if ( myws.pdf(pdfName.c_str()) ) { 
      myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI"), LineStyle(1), NumCPU(32), LineColor(kBlue+3), Precision(1e-4) );
    }
  }

  if (incPsi2S) {
    string pdfName = ""; string dsName = "";
    if (incPrompt)  { pdfName = Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"));      dsName = "binSubtractedSIG"; }
    myws.data(dsName.c_str())->plotOn(frame, Name("PSI2SDATA"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlue+3), MarkerSize(1.0));
    if ( myws.pdf(pdfName.c_str()) ) { 
      myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S"), LineStyle(1), NumCPU(32), LineColor(kBlue+3), Precision(1e-4) );
    }
  }
  if (incBkg) {

    string pdfName = ""; string dsName = "";
    if (incPrompt)  { pdfName = Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"));     dsName = "binDataCtErrSB"; }
    //if (!incJpsi && !incPsi2S) { dsName = Form("hCtauErr_Tot_%s", (isPbp?"Pbp":"PP")); } 
    myws.data(dsName.c_str())->plotOn(frame, Name("BKGDATA"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kGreen+3));
    if ( myws.pdf(pdfName.c_str()) ) { 
      myws.pdf(pdfName.c_str())->plotOn(frame,Name("BKG"), LineStyle(1), NumCPU(32), LineColor(kGreen+3), Precision(1e-4) );
    }
  }
 
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  // set the CMS style
  setTDRStyle();
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cCtauErrFig_%s", (isPbp?"Pbp":"PP")), "cCtauErrFig",1000,1000);
  TLine   *pline = new TLine(cut.dMuon.ctauErr.Min, 0.0, cut.dMuon.ctauErr.Max, 0.0);
  
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleOffset(1);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.7);
  frame->GetYaxis()->SetTitleFont(42);
  setCtauErrorRange(myws, frame, dsOSName, setLogScale);
 
  cFig->cd();
  //plot fit
  frame->Draw();

  cFig->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, "Soft Muon ID"); dy+=0.045;
  if (cutCtau) { t->DrawLatex(0.21, 0.86-dy, "#font[12]{l}_{J/#psi} cuts applied"); dy+=0.045; }
  if (isPbp) {
    t->DrawLatex(0.21, 0.86-dy, "HLT PAL1DoubleMuOpen v1"); dy+=0.045;
  } else {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } 
  
  //if (isPbp) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}
  float rapMin = cut.dMuon.AbsRap.Min;
  float rapMax = cut.dMuon.AbsRap.Max;
  if (isPbp) {cut.dMuon.AbsRap.Min = -1*(rapMax + 0.47); cut.dMuon.AbsRap.Max = -1*(rapMin + 0.47);}
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.2f #leq y_{CM} < %.2f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=1.5*0.045;
  
  // Drawing the Legend
  double ymin = 0.7802;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if(incJpsi && frame->findObject("JPSI")) { leg->AddEntry(frame->findObject("JPSI"),"J/#psi PDF","l"); }
  if(incPsi2S && frame->findObject("PSI2S")) { leg->AddEntry(frame->findObject("PSI2S"),"#psi(2S) PDF","l"); }
  if(incBkg && frame->findObject("BKG")) { leg->AddEntry(frame->findObject("BKG"),"Background","l"); }
  leg->Draw("same");

  //Drawing the title
  TString label;
  if (isPbp) {
    if (opt.Pbp.RunNb.Start==opt.Pbp.RunNb.End){
      label = Form("Pbp Run %d", opt.Pbp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "Pbp", "HIOniaL1DoubleMu0", opt.Pbp.RunNb.Start, opt.Pbp.RunNb.End);
    }
  } else {
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  }
  
  CMS_lumi(cFig, isPbp ? 108 : 107, 33, "");
  gStyle->SetTitleFontSize(0.05);
  
  cFig->Update();
  cFig->cd(); 

 
  bool SB = (incBkg&&(!incPsi2S&&!incJpsi));
  // Save the plot in different formats
  gSystem->mkdir(Form("%splot/%s/ctauErr%s/root/", outputDir.c_str(), DSTAG.c_str(), SB?"SB":""), kTRUE); 
  cFig->SaveAs(Form("%splot/%s/ctauErr%s/root/%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), DSTAG.c_str(), SB?"SB":"", DSTAG.c_str(),  "Psi2SJpsi", (isPbp?"Pbp":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%splot/%s/ctauErr%s/png/", outputDir.c_str(), DSTAG.c_str(), SB?"SB":""), kTRUE);
  cFig->SaveAs(Form("%splot/%s/ctauErr%s/png/%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), DSTAG.c_str(), SB?"SB":"", DSTAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%splot/%s/ctauErr%s/pdf/", outputDir.c_str(), DSTAG.c_str(), SB?"SB":""), kTRUE);
  cFig->SaveAs(Form("%splot/%s/ctauErr%s/pdf/%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), DSTAG.c_str(), SB?"SB":"", DSTAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  
  cFig->Clear();
  cFig->Close();

}

void setCtauErrorRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("ctauErr"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 0.0;

  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Ydown = max(0.1, YMin);
    Yup = YMax*100.0;
  }
  else
  {
    Ydown = YMin;
    Yup = YMax*1.4;
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;

}


#endif // #ifndef drawCtauErrorPlot_C
