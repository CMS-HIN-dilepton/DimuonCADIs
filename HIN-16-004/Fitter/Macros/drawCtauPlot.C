#ifndef drawCtauPlot_C
#define drawCtauPlot_C

#include "Utilities/initClasses.h"

void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, int nBins, bool setLogScale, double dMuonYmin = -1.);
void printCtauParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted);
void printCtauChi2(RooWorkspace& myws, TPad* Pad, RooPlot* frame, string varLabel, string dataLabel, string pdfLabel, int nBins, bool isWeighted); 

void drawCtauPlot(RooWorkspace& myws,   // Local workspace
                  string outputDir,     // Output directory
		  struct InputOpt opt,  // Variable with run information (kept for legacy purpose)
                  struct KinCuts cut,   // Variable with current kinematic cuts
                  string plotLabel,     // The label used to define the output file name
                  // Select the type of datasets to fit
                  string DSTAG,         // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                  bool isPbPb,          // Define if it is PbPb (True) or PP (False)
                  // Select the type of object to fit
                  bool fitMass,         // Indicates if the mass is fitted
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


  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string pdfName  = "pdfCTAU";
  if (fitMass) { pdfName = "pdfCTAUMASS"; }
  string pdfTotName  = Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"));
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
    
  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();

  // Create the main plot of the fit
  RooPlot*   frame     = myws.var("ctau")->frame(Bins(nBins), Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max));
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  Double_t norm = 1.0; //myws.data(dsOSName.c_str())->sumEntries();
    
  if (incJpsi) {
    if (incBkg) {
      if (incNonPrompt && incPrompt) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSINOPR"), Components( RooArgSet(
                                                                                           *myws.pdf(Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))),
                                                                                           *myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(),  (isPbPb?"PbPb":"PP")))
                                                                                           )
                                                                                 ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                                             );
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSIPR"),Components( RooArgSet(
                                                                                        *myws.pdf(Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))),
                                                                                        *myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                        )
                                                                              ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32), 
                                             LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                                             );
      }
    } else {
      if (incNonPrompt && incPrompt) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSINOPR"), Components( RooArgSet(
                                                                                           *myws.pdf(Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                           )
                                                                                 ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                                             );
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSIPR"),Components( RooArgSet(
                                                                                        *myws.pdf(Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                        )
                                                                              ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kOrange+2), LineStyle(1), Precision(1e-4), Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max)
                                             );
      }
    }
  }
  if (incPsi2S) {
    if (incBkg) {
      if (incNonPrompt && incPrompt) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SNOPR"), Components( RooArgSet(
                                                                                           *myws.pdf(Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))),
                                                                                           *myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(),  (isPbPb?"PbPb":"PP")))
                                                                                           )
                                                                                 ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                                             );
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SPR"),Components( RooArgSet(
                                                                                         *myws.pdf(Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))),
                                                                                         *myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                         )
                                                                               ),
                                             Normalization(norm, RooAbsReal::NumEvent), LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                                             );
      }
    } else {
      if (incNonPrompt && incPrompt) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SNOPR"), Components( RooArgSet(
                                                                                           *myws.pdf(Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                           )
                                                                                 ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                                             );
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SPR"),Components( RooArgSet(
                                                                                         *myws.pdf(Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))
                                                                                         )
                                                                               ),
                                             Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32),
                                             LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                                             );
      }
    }
  }
  if (incBkg && (!incJpsi && !incPsi2S)) {
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))))),
                                         Normalization(myws.data(dsOSName.c_str())->reduce("invMass<2.8 ||invMass>4.0")->sumEntries(), RooAbsReal::NumEvent), 
                                         Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max), FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue)
                                         );
  } 
  if (incBkg && (incJpsi || incPsi2S)) {
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"))))),
                                         Normalization(norm, RooAbsReal::NumEvent),
                                         Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max), FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue)
                                         );
  } 
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  if (incBkg && (!incJpsi && !incPsi2S)) {
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(norm, RooAbsReal::NumEvent), 
                                         ConditionalObservables(*myws.var("ctauErr")), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32), 
                                         LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max));
  } else {
    myws.pdf(pdfTotName.c_str())->plotOn(frame, Name("PDF"),  Normalization(norm, RooAbsReal::NumEvent), ProjWData(*((RooAbsData*)myws.data(dsOSName.c_str()))), NumCPU(32), 
                                         LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max)
                                         );
  }
  /*
  // Create the pull distribution of the fit
  RooHist *hpull = frame->pullHist(0, 0, true);
  hpull->SetName("hpull");
  */
  RooPlot* frame2 = myws.var("ctau")->frame(Title("Pull Distribution"), Bins(nBins), Range(cut.dMuon.ctau.Min, cut.dMuon.ctau.Max));
  /*
  frame2->addPlotable(hpull, "PX"); 			
  */
  // set the CMS style
  setTDRStyle();
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cCtauFig_%s", (isPbPb?"PbPb":"PP")), "cCtauFig",800,800);
  TPad    *pad1  = new TPad(Form("pad1_%s", (isPbPb?"PbPb":"PP")),"",0,0.23,1,1);
  TPad    *pad2  = new TPad(Form("pad2_%s", (isPbPb?"PbPb":"PP")),"",0,0,1,.228);
  TLine   *pline = new TLine(cut.dMuon.ctau.Min, 0.0, cut.dMuon.ctau.Max, 0.0);
  
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetTitleOffset(3);
  frame->GetXaxis()->SetLabelOffset(3);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.7);
  frame->GetYaxis()->SetTitleFont(42);
  setCtauRange(myws, frame, dsOSName, nBins, setLogScale, cut.dMuon.AbsRap.Min);
 
  cFig->cd();
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000); 
  pad1->SetBottomMargin(0.015); 
  //plot fit
  pad1->Draw();
  pad1->cd(); 
  frame->Draw();

  printCtauParameters(myws, pad1, isPbPb, pdfTotName, isWeighted);
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, "2015 HI Soft Muon ID"); dy+=0.045;
  if (cutCtau) { t->DrawLatex(0.21, 0.86-dy, "c#tau^{J/#psi} cuts applied"); dy+=0.045; }
  if (isPbPb) {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } else {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } 
  if (isPbPb) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=1.5*0.045;

  // Drawing the Legend
  double ymin = 0.7802;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if(frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),"Total fit","l"); }
  if(frame->findObject("BKG")) { leg->AddEntry(frame->findObject("BKG"),"Background","fl"); }
  leg->Draw("same");

  //Drawing the title
  TString label;
  if (isPbPb) {
    if (opt.PbPb.RunNb.Start==opt.PbPb.RunNb.End){
      label = Form("PbPb Run %d", opt.PbPb.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PbPb", "HIOniaL1DoubleMu0", opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End);
    }
  } else {
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  }
  
  CMS_lumi(pad1, isPbPb ? 105 : 104, 33, label);
  gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  //---plot pull
  pad2->Draw();
  pad2->cd();
    
  frame2->SetTitle("");
  frame2->GetYaxis()->CenterTitle(kTRUE);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.1);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetXaxis()->CenterTitle(kTRUE);
  frame2->GetXaxis()->SetTitleOffset(1);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetXaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetTitle("c#tau_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

  frame2->Draw(); 
  
  // *** Print chi2/ndof 
  //printCtauChi2(myws, pad2, frame, "ctau", dsOSName.c_str(), pdfTotName.c_str(), nBins, isWeighted);
  
  pline->Draw("same");
  pad2->Update();
  
  // Save the plot in different formats
  gSystem->mkdir(Form("%splot/%s/ctau/root/", outputDir.c_str(), DSTAG.c_str()), kTRUE); 
  cFig->SaveAs(Form("%splot/%s/ctau/root/%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), DSTAG.c_str(), DSTAG.c_str(),  "Psi2SJpsi", (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%splot/%s/ctau/png/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/%s/ctau/png/%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), DSTAG.c_str(), DSTAG.c_str(), "Psi2SJpsi", (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%splot/%s/ctau/pdf/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/%s/ctau/pdf/%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), DSTAG.c_str(), DSTAG.c_str(), "Psi2SJpsi", (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  
  cFig->Clear();
  cFig->Close();

}

#endif // #ifndef drawCtauPlot_C


void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, int nBins, bool setLogScale, double dMuonYmin)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("ctau"), Binning(nBins));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) isMC = true;
    
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    if (isMC) Ydown = YMin*0.3;
    else Ydown = max(1.0, YMin/(TMath::Power((YMax/YMin), 0.1)));
    Yup = YMax*TMath::Power((YMax/YMin), 0.5);
  }
  else
  {
    Ydown = max(YMin-(YMax-YMin)*0.2,0.0);
    Yup = YMax+(YMax-YMin)*0.5;
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;

}


void printCtauParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.025; 
  RooArgSet* Parameters =  myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")));
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="invMass"){continue;} else if(s1=="ctau"){continue;} else if(s1=="MassRatio"){continue;}   
    else if(s1=="One"){continue;} else if(s1=="mMin"){continue;} else if(s1=="mMax"){continue;}
    else if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){
      s1=Form("#%s",s1.c_str());
    }

    if(s2=="CtauRes")  { s2="Res";   } 
    else if(s2=="JpsiNoPR")  { s2="J/#psi[NoPR]";   } 
    else if(s2=="JpsiPR")  { s2="J/#psi[PR]";   } 
    else if(s2=="Jpsi" && s1=="N")  { s2="J/#psi";   } 
    else if(s2=="Psi2SNoPR") { s2="#psi(2S)[NoPR]"; }
    else if(s2=="Psi2SPR") { s2="#psi(2S)[PR]"; } 
    else if(s2=="Psi2S" && s1=="N")  { s2="#psi(2S)";   } 
    else if(s2=="BkgNoPR")   { s2="bkg[NoPR]";      } 
    else if(s2=="BkgPR")   { s2="bkg[PR]";      } 
    else if(s2=="Bkg" && s1=="N")   { s2="bkg";      }
    else {continue;}
    if(s3!=""){
      label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());
    }
    // Print the parameter's results
    if(s1=="N"){ 
      t->DrawLatex(0.69, 0.75-dy, Form((isWeighted?"%s = %.6f#pm%.6f ":"%s = %.0f#pm%.0f "), label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("sigma")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f GeV/c^{2}", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("lambda")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("ctau")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f cm", (label.insert(1, string("#"))).c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else { 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
  }
}


void printCtauChi2(RooWorkspace& myws, TPad* Pad, RooPlot* frame, string varLabel, string dataLabel, string pdfLabel, int nBins, bool isWeighted) 
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1); 
  unsigned int nFitPar = myws.pdf(pdfLabel.c_str())->getParameters(*myws.data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize(); 
  TH1 *hdatact = myws.data(dataLabel.c_str())->createHistogram("hdatact", *myws.var(varLabel.c_str()), Binning(nBins));
//  RooHist *hpull = frame->pullHist("hdatact",pdfLabel.c_str(), true);
  RooHist *hpull = frame->pullHist(0,0, true);
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i+1) > 0.0) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  //chi2 = myws.pdf(pdfLabel.c_str())->createChi2(*((RooDataSet*)myws.data(dataLabel.c_str())))->getVal(); 
  //chi2 = frame->chiSquare(nFitPar)*ndof;
//  RooDataHist dummy("dummy", "dummy", *myws.var("invMass"), hdatact);
//  if (isWeighted) {
//    chi2 = RooChi2Var("chi2", "chi2", *myws.pdf(pdfLabel.c_str()), dummy, kFALSE, 0, 0, 8, RooFit::Interleave, kFALSE, kFALSE, RooDataHist::SumW2).getVal();
//  }  
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  delete hdatact; 
  delete hpull;
};