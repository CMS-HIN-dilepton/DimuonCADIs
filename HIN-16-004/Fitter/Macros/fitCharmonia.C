#include "Utilities/initClasses.h"
#include "buildCharmoniaMassModel.C"
#include "buildCharmoniaCtauModel.C"
#include "buildCharmoniaCtauTrueModel.C"
#include "drawMassPlot.C"
#include "drawCtauMassPlot.C"
#include "drawCtauPlot.C"
#include "drawCtauErrorPlot.C"
#include "drawCtauTruePlot.C"

#include <algorithm>

void setCtauCuts(struct KinCuts& cut, bool isPbp);
int  importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts& cut, string label, bool incJpsi, bool incPsi2S, bool cutSideBand=false);
void setGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, bool fitCtau, bool fitCtauTrue, bool incJpsi, bool incPsi2S, bool incBkg);
bool setMassModel(struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg);
bool setCtauModel(struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt);
bool setCtauTrueModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp );
void setOptions(struct InputOpt* opt);
void setFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp=true, bool cutSideBand=false);
void setNameFile(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool isPbp=true, bool cutSideBand=false);
bool saveWorkSpace(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp=true, bool cutSideBand=false);
 RooDataHist* subtractSidebands(RooWorkspace ws, RooDataHist* subtrData, RooDataHist* all, RooDataHist* side, float scalefactor, string varName);

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp=true, bool cutSideBand=false);
bool loadPreviousFitResult(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp=true, bool cutSideBand=false);
bool loadFinalFitResult(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue,  bool incJpsi, bool incPsi2S, bool isPbp=true, bool cutSideBand=false);
bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2);
bool makeCtauErrPdf(RooWorkspace& ws, string dsName, bool incPsi2S, bool isPbp, bool fitMass);
bool makeMCCtauErrPdf(RooWorkspace& ws, string dsName, bool incPsi2S, bool isPbp);
bool fitCharmonia( RooWorkspace&  inputWorkspace, // Workspace with all the input RooDatasets
		   struct KinCuts cut,            // Variable containing all kinematic cuts
		   map<string, string>  parIni,   // Variable containing all initial parameters
		   string outputDir,              // Path to output directory
                   // Select the type of datasets to fit
		   string DSTAG,                  // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
		   bool isPbp      = false,      // isPbp = false for pp, true for Pbp
                   // Select the type of object to fit                   
                   bool fitMass      = true,      // Fit mass distribution
                   bool fitCtau      = false,     // Fit ctau distribution
                   bool fitCtauTrue  = false,     // Fit ctau truth MC distribution
                   bool incJpsi      = true,      // Includes Jpsi model
                   bool incPsi2S     = true,      // Includes Psi(2S) model
                   bool incBkg       = true,      // Includes Background model
                   bool incPrompt    = true,      // Includes Prompt ctau model
                   bool incNonPrompt = false,     // Includes NonPrompt ctau model
                   // Select the fitting options
		   bool isErrPdfData = false,        // ctauErr pdf from data or MC
                   bool cutCtau      = false,     // Apply prompt ctau cuts
                   bool doSimulFit   = false,     // Do simultaneous fit
                   bool wantPureSMC  = false,     // Flag to indicate if we want to fit pure signal MC
                   int  numCores     = 32,         // Number of cores used for fitting
                   // Select the drawing options
                   bool  setLogScale  = true,      // Draw plot with log scale
                   bool  incSS        = false,     // Include Same Sign data
                   bool  zoomPsi      = false,     // Zoom Psi(2S) peak on extra pad
                   float binWidth     = 0.05,      // Bin width used for plotting
                   bool  getMeanPT    = false      // Compute the mean PT (NEED TO FIX)
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
  
  //RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-8) ;
  //RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
  // Define the mass range
  //if (cut.dMuon.M.Max==5 && cut.dMuon.M.Min==2) { 
  // Default mass values, means that the user did not specify a mass range
  if ( incJpsi && !incPsi2S) {
    cut.dMuon.M.Min = 2.6;
    cut.dMuon.M.Max = 3.5;
    
    //cut.dMuon.ctauErr.Min = 0.02; 
    //cut.dMuon.ctauErr.Max = 0.15;
    cut.dMuon.ctau.Min = -1.5;
    cut.dMuon.ctau.Max = 3.0;

  }
  else if ( !incJpsi && incPsi2S) {
    cut.dMuon.M.Min = 3.4;
    cut.dMuon.M.Max = 4.2;
    cut.dMuon.ctau.Min = -1.5;
    cut.dMuon.ctau.Max = 3.0;
    
  }
  else {
    //cut.dMuon.M.Min = 2.6;
    //cut.dMuon.M.Max = 3.5;
    cut.dMuon.M.Min = 3.4;
    cut.dMuon.M.Max = 4.2;
    
    //cut.dMuon.ctauErr.Min = 0.01; 
    //cut.dMuon.ctauErr.Max = 0.3;
    cut.dMuon.ctau.Min = -1.5;
    cut.dMuon.ctau.Max = 3.0;
  }
  
  parIni["invMassNorm"] = Form("RooFormulaVar::%s('( -1.0 + 2.0*( @0 - @1 )/( @2 - @1) )', {%s, mMin[%.6f], mMax[%.6f]})", "invMassNorm", "invMass", cut.dMuon.M.Min, cut.dMuon.M.Max );
  bool cutSideBand = false;
  if (incBkg && (!incPsi2S && !incJpsi)) { cutSideBand = true; }
  
  // Define the mass range
  if (fitCtau) {
    if (cut.dMuon.ctau.Min==-50. && cut.dMuon.ctau.Max==50.) { 
      // Default ctau values, means that the user did not specify a ctau range
      if (incNonPrompt) {
        cut.dMuon.ctau.Min = -1.5;
        cut.dMuon.ctau.Max = 3.0;
      } else {
        cut.dMuon.ctau.Min = -1.5;
        cut.dMuon.ctau.Max = 3.0;
      }
    }
    if (cut.dMuon.ctauErr.Min==-50. && cut.dMuon.ctauErr.Max==50.) { 
      // Default ctau values, means that the user did not specify a ctau Err range
      cut.dMuon.ctauErr.Min = 0.01; 
      cut.dMuon.ctauErr.Max = 0.3;
    }
  }
  // Define the mass range
  if (fitCtauTrue) {
    if (cut.dMuon.ctauTrue.Min==-50.0 && cut.dMuon.ctauTrue.Max==50.0) { 
      // Default ctau values, means that the user did not specify a ctau True range
      cut.dMuon.ctauTrue.Min = -1.0;
      cut.dMuon.ctauTrue.Max = 10.0;
    }
  }
  // Apply the ctau cuts to reject non-prompt charmonia
  if (cutCtau) { setCtauCuts(cut, isPbp); }  

  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    if (incJpsi && incPsi2S) { 
      cout << "[ERROR] We can not fit Jpsi and Psi(2S) simultaneously using MC" << endl; return false; 
    }
    if (incPrompt && incNonPrompt) { 
      cout << "[ERROR] We can not fit Prompt and Non-Prompt simultaneously using MC" << endl; return false; 
    }
    isMC = true;
  }
  if (isMC && wantPureSMC) wantPureSMC=true;
  else wantPureSMC=false;

  struct InputOpt opt; setOptions(&opt);
  
  string plotMassLabelPbp,  plotMassLabelPP, plotCtauLabelPbp,  plotCtauLabelPP, plotCtauTrueLabelPbp,  plotCtauTrueLabelPP;

  struct OniaModel model;
  RooWorkspace   myws("workspace", "local workspace");
  bool doFit = true;

  if (!isPbp) {
    
    // Set models based on initial parameters
    if (fitMass)     { if (!setMassModel(model, parIni, false, incJpsi, incPsi2S, incBkg)) { return false; } }
    if (fitCtau)     { if (!setCtauModel(model, parIni, false, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt)) { return false; } }
    if (fitCtauTrue) { if (!setCtauTrueModel(model, parIni, false)) { return false; } }
    
    // Import the local datasets
    string label = Form("%s_%s", DSTAG.c_str(), "PP");
    if (wantPureSMC) label = Form("%s_%s_NoBkg", DSTAG.c_str(), "PP");
    string dsName = Form("dOS_%s", label.c_str());
    int importID = importDataset(myws, inputWorkspace, cut, label, incJpsi, incPsi2S, cutSideBand);
    if (importID<0) { return false; }
    else if (importID==0) { doFit = false; }
    setGlobalParameterRange(myws, parIni, cut, label, fitCtau, fitCtauTrue, incJpsi, incPsi2S, incBkg);
    

    // Build the Fit Model    
    double numEntries = myws.data(dsName.c_str())->sumEntries();
    if (fitMass)     { if (!buildCharmoniaMassModel(myws, model.PP, parIni, false, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries))  { return false; } }
    if (fitCtauTrue) { if (!buildCharmoniaCtauTrueModel(myws, model.PP, parIni, false, incJpsi, incPsi2S, numEntries))  { return false; } }

    if (fitMass) {
      if (incJpsi)  { plotMassLabelPP = plotMassLabelPP + Form("_Jpsi_%s", parIni["Model_Jpsi_PP"].c_str());   } 
      if (incPsi2S) { plotMassLabelPP = plotMassLabelPP + Form("_Psi2S_%s", parIni["Model_Psi2S_PP"].c_str()); }
      if (incBkg)   { plotMassLabelPP = plotMassLabelPP + Form("_Bkg_%s", parIni["Model_Bkg_PP"].c_str());     }
    }
    
    if (fitCtau) {
      if (incJpsi  && incPrompt)    { plotCtauLabelPP = plotCtauLabelPP + Form("_JpsiPR_%s", parIni["Model_JpsiPR_PP"].c_str());      } 
      if (incJpsi  && incNonPrompt) { plotCtauLabelPP = plotCtauLabelPP + Form("_JpsiNoPR_%s", parIni["Model_JpsiNoPR_PP"].c_str());  }
      if (incPsi2S && incPrompt)    { plotCtauLabelPP = plotCtauLabelPP + Form("_Psi2SPR_%s", parIni["Model_Psi2SPR_PP"].c_str());    } 
      if (incPsi2S && incNonPrompt) { plotCtauLabelPP = plotCtauLabelPP + Form("_Psi2SNoPR_%s", parIni["Model_Psi2SNoPR_PP"].c_str());}
      if (incBkg && incPrompt)      { plotCtauLabelPP = plotCtauLabelPP + Form("_BkgPR_%s", parIni["Model_BkgPR_PP"].c_str());        } 
      if (incBkg && incNonPrompt)   { plotCtauLabelPP = plotCtauLabelPP + Form("_BkgNoPR_%s", parIni["Model_BkgNoPR_PP"].c_str());    }
    }
    
    if (fitCtauTrue) {
      plotCtauTrueLabelPP = string("") + Form("_CtauTrueRes_%s", parIni["Model_CtauTrueRes_PP"].c_str()) + Form("_CtauTrue_%s", parIni["Model_CtauTrue_PP"].c_str());
    }
  }
  if (isPbp) {
    
    // Set models based on initial parameters
    if (fitMass)     { if (!setMassModel(model, parIni, true, incJpsi, incPsi2S, incBkg)) { return false; } }
    if (fitCtau)     { if (!setCtauModel(model, parIni, true, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt)) { return false; } }
    if (fitCtauTrue) { if (!setCtauTrueModel(model, parIni, true)) { return false; } }
    
    // Import the local datasets
    string label = Form("%s_%s", DSTAG.c_str(), "Pbp");
    if (wantPureSMC) label = Form("%s_%s_NoBkg", DSTAG.c_str(), "Pbp");
    string dsName = Form("dOS_%s", label.c_str());
    int importID = importDataset(myws, inputWorkspace, cut, label, incJpsi, incPsi2S, cutSideBand);
    if (importID<0) { return false; }
    else if (importID==0) { doFit = false; }
    setGlobalParameterRange(myws, parIni, cut, label, fitCtau, fitCtauTrue, incJpsi, incPsi2S, incBkg);
    
    // Build the Fit Model
    double    numEntries = myws.data(dsName.c_str())->sumEntries();
    if (fitMass) { if (!buildCharmoniaMassModel(myws, model.Pbp, parIni, true, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries)) { return false; } }
    if (fitCtauTrue) { if (!buildCharmoniaCtauTrueModel(myws, model.Pbp, parIni, true, incJpsi, incPsi2S, numEntries))  { return false; } }

    if (fitMass) {
      if (incJpsi)  { plotMassLabelPbp = plotMassLabelPbp + Form("_Jpsi_%s", parIni["Model_Jpsi_Pbp"].c_str());   } 
      if (incPsi2S) { plotMassLabelPbp = plotMassLabelPbp + Form("_Psi2S_%s", parIni["Model_Psi2S_Pbp"].c_str()); }
      if (incBkg)   { plotMassLabelPbp = plotMassLabelPbp + Form("_Bkg_%s", parIni["Model_Bkg_Pbp"].c_str());     }
    }
    if (fitCtau) {
      if (incJpsi  && incPrompt)    { plotCtauLabelPbp = plotCtauLabelPbp + Form("_JpsiPR_%s", parIni["Model_JpsiPR_Pbp"].c_str());      } 
      if (incJpsi  && incNonPrompt) { plotCtauLabelPbp = plotCtauLabelPbp + Form("_JpsiNoPR_%s", parIni["Model_JpsiNoPR_Pbp"].c_str());  }
      if (incPsi2S && incPrompt)    { plotCtauLabelPbp = plotCtauLabelPbp + Form("_Psi2SPR_%s", parIni["Model_Psi2SPR_Pbp"].c_str());    } 
      if (incPsi2S && incNonPrompt) { plotCtauLabelPbp = plotCtauLabelPbp + Form("_Psi2SNoPR_%s", parIni["Model_Psi2SNoPR_Pbp"].c_str());}
      if (incBkg && incPrompt)      { plotCtauLabelPbp = plotCtauLabelPbp + Form("_BkgPR_%s", parIni["Model_BkgPR_Pbp"].c_str());        } 
      if (incBkg && incNonPrompt)   { plotCtauLabelPbp = plotCtauLabelPbp + Form("_BkgNoPR_%s", parIni["Model_BkgNoPR_Pbp"].c_str());    }
    }
    if (fitCtauTrue) {
      plotCtauTrueLabelPbp = string("") + Form("_CtauTrueRes_%s", parIni["Model_CtauTrueRes_Pbp"].c_str()) + Form("_CtauTrue_%s", parIni["Model_CtauTrue_Pbp"].c_str());
    }
  }
  
  if (doFit) {
    if (isPbp) {
      string plotLabel = ""; string pdfName = "";
      string dsName = Form("dOS_%s_Pbp", DSTAG.c_str());
      if (wantPureSMC) dsName = Form("dOS_%s_Pbp_NoBkg", DSTAG.c_str());
      
      if (fitCtauTrue) {
	pdfName = "pdfCTAUTRUE_Tot_Pbp"; plotLabel = plotCtauTrueLabelPbp;
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauTrue")));
	bool found =  true; bool skipFit = false;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, false, false, fitCtauTrue, false, true);
	if (found) {
	  if (fitCtau) {
	    cout << "[INFO] This ctauTrue fit was already done, so I'll load the fit results." << endl;
	    if ( loadPreviousFitResult(myws, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, false, false, fitCtauTrue, false, true, false) ) {
	      skipFit = true;
	    }
	  } else {
	    cout << "[INFO] This ctauTrue fit was already done, so I'll just go to the next one." << endl;
	    return true;
	  }
	}
	if (skipFit==false) {
	  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("FullWindow"), NumCPU(numCores), Save());
	  fitResult->Print();
	  // Draw the mass plot
	  int nBins = min(int( round((cut.dMuon.ctauTrue.Max - cut.dMuon.ctauTrue.Min)/binWidth) ), 1000);
	  drawCtauTruePlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, true, incJpsi, incPsi2S, cutCtau, wantPureSMC, setLogScale, incSS, nBins); 
	  saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, false, false, fitCtauTrue, false, true, cutSideBand);
	  myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
	}
	if (fitCtau) {
	  myws.pdf("pdfCTAUTRUE_Tot_Pbp")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
	}
      }
      if (fitMass && !fitCtau) {
	pdfName = "pdfMASS_Tot_Pbp"; plotLabel = plotMassLabelPbp;
        
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass")));
	bool found =  true; bool skipFit = false;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, fitMass, false, false, false, true, cutSideBand);
	if (found) {
	  if (fitCtau) {
	    cout << "[INFO] This mass fit was already done, so I'll load the fit results." << endl;
	    string plotMassLabel = ""; 
	    if (incJpsi)  { plotMassLabel = plotMassLabel + Form("_Jpsi_%s", parIni["Model_Jpsi_Pbp"].c_str());   } 
	    if (incPsi2S) { plotMassLabel = plotMassLabel + Form("_Psi2S_%s", parIni["Model_Psi2S_Pbp"].c_str()); }
	    if (incBkg)   { plotMassLabel = plotMassLabel + Form("_Bkg_%s", parIni["Model_Bkg_Pbp"].c_str());     }
	    if ( loadPreviousFitResult(myws, outputDir, plotMassLabel, DSTAG, cut, fitMass, false, false, false, true, cutSideBand) ) {
	      skipFit = true;
	    }
	  } else {
	    cout << "[INFO] This mass fit was already done, so I'll just go to the next one." << endl;
	    return true;
	  }
	}
	if (skipFit==false) {
	  RooFitResult* fitResult = myws.pdf("pdfMASS_Tot_Pbp")->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores), Range((incJpsi||incPsi2S)?"FullWindow":parIni["BkgMassRange_FULL_Label"].c_str()), Save());  
	  //RooFitResult* fitResult = myws.pdf("pdfMASS_Tot_Pbp")->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores), Save());  
	  fitResult->Print();
	  // Draw the mass plot
	  int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
	  drawMassPlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, true, incJpsi, incPsi2S, incBkg, cutCtau, wantPureSMC, setLogScale, incSS, zoomPsi, nBins, getMeanPT);

	  if(isErrPdfData)
	    if( !makeCtauErrPdf(myws, dsName, incPsi2S ,isPbp, fitMass) ) { return false; }
	  
	  saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, fitMass, false, false, false, true, cutSideBand);
	  myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
	}
	
	if (fitCtau) {
	  myws.pdf("pdfMASS_Tot_Pbp")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
	}
      }

      if (isMC && incPrompt &&(!isErrPdfData)){
	if( !makeMCCtauErrPdf(myws, dsName, incPsi2S, isPbp) ) { return false; }
      }
      
      //myws.Print("v");
      if (fitCtau &&(isMC || fitMass) &&isErrPdfData){
	string plotMassLabel;
	if (incJpsi) plotMassLabel = "_Jpsi_GaussianAndCrystalBall_Bkg_Exponential";
	if (incPsi2S) plotMassLabel = "_Psi2S_SingleCrystalBall_Bkg_Exponential";
	loadFinalFitResult(myws, outputDir, plotMassLabel, "DATA", cut, fitMass, fitCtau, false,  incJpsi, incPsi2S, true, cutSideBand);
	cout<<"outputDir :"<<outputDir<<endl;
      }

      
      if (fitCtau) { 
	if (!fitMass && fitCtau) { pdfName = "pdfCTAU_Tot_Pbp";     plotLabel = plotCtauLabelPbp;                   }
	if (fitMass &&  fitCtau) { pdfName = "pdfCTAUMASS_Tot_Pbp"; plotLabel = plotMassLabelPbp+plotCtauLabelPbp; }

	cout << "====================ok12========================================" << endl;
	
	if (fitMass) {
	  string plotMassLabel;
	  if (incJpsi)  plotMassLabel = "_Jpsi_GaussianAndCrystalBall_Bkg_Exponential";
	  if (incPsi2S) plotMassLabel = "_Psi2S_SingleCrystalBall_Bkg_Exponential"; 
	  //plotMassLabel = plotMassLabel + Form("_Jpsi_%s", parIni["Model_Jpsi_Pbp"].c_str());
	  //plotMassLabel = plotMassLabel + Form("_Psi2S_%s", parIni["Model_Psi2S_Pbp"].c_str());
	  //plotMassLabel = plotMassLabel + Form("_Bkg_%s", parIni["Model_Bkg_Pbp"].c_str());
	  loadPreviousFitResult(myws, outputDir, plotMassLabel, DSTAG, cut, true, false, false, false, true, cutSideBand);
	  myws.var("invMass")->setRange("InclusiveMassRange", cut.dMuon.M.Min, cut.dMuon.M.Max);
	  
	  if (myws.pdf("pdfMASS_Tot_Pbp")) {
	    myws.pdf("pdfMASS_Tot_Pbp")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE); 
	  }

	  cout << "====================ok12a========================================" << endl;
	  
	  /*if (myws.var("N_Bkg_Pbp")) { 
	    myws.var("N_Bkg_Pbp")->setConstant(kFALSE); 
	    }

	  if (myws.var("N_Jpsi_Pbp")) { 
	    myws.var("N_Jpsi_Pbp")->setConstant(kFALSE); 
	  }

	  if (myws.var("N_Psi2S_Pbp")) { 
	    myws.var("N_Psi2S_Pbp")->setConstant(kFALSE); 
	    }*/
	}

	if (fitCtau && incBkg ){
	  if(!myws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))))
	    if( !makeCtauErrPdf(myws, dsName, incPsi2S ,isPbp, fitMass) ) { return false; }
	}

	cout << "====================ok12b========================================" << endl;
	
	if (fitMass &&  fitCtau){
	  //string plotCtauTrueLabel = string("") + Form("_CtauTrueRes_%s", parIni["Model_CtauTrueRes_Pbp"].c_str()) + Form("_CtauTrue_%s", parIni["Model_CtauTrue_Pbp"].c_str());
	  string plotCtauTrueLabel = "_CtauTrueRes_SingleGaussianResolution_CtauTrue_SingleSidedDecay";
	  if (incJpsi) loadFinalFitResult(myws, outputDir, plotCtauTrueLabel, "MCJPSINOPR", cut, false, false, true, incJpsi, incPsi2S, true, false);
	  if (incPsi2S) loadFinalFitResult(myws, outputDir, plotCtauTrueLabel, "MCPSI2SNOPR", cut, false, false, true, incJpsi, incPsi2S, true, false);
	  
	  //string plotCtauLabel = Form("_BkgPR_%s", (parIni["Model_BkgPR_Pbp"]==""?"Delta":parIni["Model_BkgPR_Pbp"].c_str()));
	  //loadFinalFitResult(myws, outputDir, plotCtauLabel, "DATA", cut, false, true, false, incJpsi, incPsi2S, true, true);
	}                           
	double numEntries = myws.data(dsName.c_str())->sumEntries();
	if (!buildCharmoniaCtauModel(myws, model.Pbp, parIni, dsName, true, incBkg, incJpsi, incPsi2S, incPrompt, incNonPrompt, isErrPdfData, numEntries)) { return false; }
	
	cout << "====================ok13========================================" << endl;

	cout<<" pdfName.c_str() "<<pdfName.c_str()<<endl;
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters( RooArgSet(*myws.var("invMass"), *myws.var("ctau")  ));
	bool found =  true;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, fitMass, fitCtau, false, false, true, cutSideBand);
	if (found) {
	  cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
	  return true;
	}

	cout << "====================ok14========================================" << endl;
	
	if (!isMC){
	  if (incPrompt && fitMass) {
	    if(incJpsi) {
	      string plotCtauLabel = Form("_JpsiPR_%s", (parIni["Model_JpsiPR_Pbp"]==""?"Delta":parIni["Model_JpsiPR_Pbp"].c_str()));
	      loadPreviousFitResult(myws, outputDir, plotCtauLabel, "MCJPSIPR", cut, false, true, false, false, true, false);
	      if(myws.pdf("pdfCTAUCOND_JpsiPR_Pbp")) {
		myws.pdf("pdfCTAUCOND_JpsiPR_Pbp")->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	      }
	      
	      RooArgSet* model_params1 = myws.pdf("pdfCTAUCOND_JpsiPR_Pbp")->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"))) ;
	      //model_params1->Print("v");
	    }
	    
	    if(incPsi2S) {
	      string plotCtauLabel = Form("_Psi2SPR_%s", (parIni["Model_Psi2SPR_Pbp"]==""?"Delta":parIni["Model_Psi2SPR_Pbp"].c_str()));
	      loadPreviousFitResult(myws, outputDir, plotCtauLabel, "MCPSI2SPR", cut, false, true, false, false, true, false);                      
	      if(myws.pdf("pdfCTAUCOND_Psi2SPR_Pbp") || myws.pdf("pdfCTAURES_Psi2SPR_Pbp")) { 
		//myws.pdf("pdfCTAUCOND_Psi2SPR_Pbp")->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	      }
	    }
	  }
	  
	  cout << "====================ok15========================================" << endl;
	  //return  true;
	  
	  if (incNonPrompt && (incJpsi || incPsi2S)) {
	    string plotCtauTrueLabel = string("") + Form("_CtauTrueRes_%s", parIni["Model_CtauTrueRes_Pbp"].c_str()) + Form("_CtauTrue_%s", parIni["Model_CtauTrue_Pbp"].c_str());
	    if(incJpsi) loadPreviousFitResult(myws, outputDir, plotCtauTrueLabel, "MCJPSINOPR", cut, false, false, true, false, true, false);
	    if ( myws.var("sigmaMC_JpsiNoPR_Pbp") )    { myws.var("sigmaMC_JpsiNoPR_Pbp")->setConstant(kTRUE);  }
	    if ( myws.var("lambdaDSS_JpsiNoPR_Pbp") )    { myws.var("lambdaDSS_JpsiNoPR_Pbp")->setConstant(kTRUE);  }
	    if(incPsi2S) loadPreviousFitResult(myws, outputDir, plotCtauTrueLabel, "MCPSI2SNOPR", cut, false, false, true, false, true, false);
	    if ( myws.var("sigmaMC_Psi2SNoPR_Pbp") )   { myws.var("sigmaMC_Psi2SNoPR_Pbp")->setConstant(kTRUE); }
	    if ( myws.var("lambdaDSS_Psi2SNoPR_Pbp") )    { myws.var("lambdaDSS_Psi2SNoPR_Pbp")->setConstant(kTRUE);  }
	    
	    string plotCtauLabel = Form("_BkgPR_%s", parIni["Model_BkgPR_Pbp"].c_str());
	    plotCtauLabel = plotCtauLabel + Form("_BkgNoPR_%s", parIni["Model_BkgNoPR_Pbp"].c_str());
	    loadPreviousFitResult(myws, outputDir, plotCtauLabel, "DATA", cut, false, true, false, false, true, true);
	    if (myws.pdf("pdfCTAUCOND_BkgPR_Pbp")) { 
	      myws.pdf("pdfCTAUCOND_BkgPR_Pbp")->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	    }
	    //if (myws.var("b_Bkg_Pbp")) { myws.var("b_Bkg_Pbp")->setConstant(kTRUE); }
	  }
	  
	  //myws.var(Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")))->setConstant(kTRUE);
	  //myws.var(Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);
	  //myws.var(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))) ->setConstant(kFALSE);
	  //myws.var(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))) ->setConstant(kFALSE);
	  //myws.var("N_Bkg_Pbp")->setConstant(kFALSE);
	  //myws.var("N_Jpsi_Pbp")->setConstant(kFALSE);
	  //myws.var("N_Psi2S_Pbp")->setConstant(kFALSE); 
	  
	}
	
	myws.Print("v");
	

	if (isMC) {
	  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores), ConditionalObservables(*myws.var("ctauErr")), Save());
	  fitResult->Print();
	} else {
	  if (fitMass){
	    cout << "[INFO] Fitting pdf: "<< pdfName << endl;
	    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kFALSE), NumCPU(numCores), Save(), ConditionalObservables(*myws.var("ctauErr")));
	    fitResult->Print("v");
	  } else {
	    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores),  Save(), ConditionalObservables(*myws.var("ctauErr")));
	    cout<<"==========good1=================="<<endl;
	    fitResult->Print("v");
	  }
	}

	cout<<"==========good2=================="<<endl;
	int nBins = min(int( round((cut.dMuon.ctau.Max - cut.dMuon.ctau.Min)/binWidth) ), 1000);
	drawCtauPlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, true, fitMass, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, cutCtau, wantPureSMC, setLogScale, incSS, nBins);
	cout<<"==========good3=================="<<endl;
	
	double width = (myws.var("ctauErr")->getMax() - myws.var("ctauErr")->getMin()) / 100.0;
	nBins = min(int( round( (cut.dMuon.ctauErr.Max - cut.dMuon.ctauErr.Min)/width )), 1000);
	if (!isMC) drawCtauErrorPlot(myws, outputDir, opt, cut, DSTAG, true, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, cutCtau, wantPureSMC, setLogScale, incSS, nBins);
	saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, fitMass, fitCtau, false, false, true, cutSideBand);
	myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
      }
    }
    else {
      cut.Centrality.Start = 0;
      cut.Centrality.End = 200;
      
      string plotLabel = ""; string pdfName = "";
      string dsName = Form("dOS_%s_PP", DSTAG.c_str());
      if (wantPureSMC) dsName = Form("dOS_%s_PP_NoBkg", DSTAG.c_str());
      
      if (fitCtauTrue) {
	
	pdfName = "pdfCTAUTRUE_Tot_PP"; plotLabel = plotCtauTrueLabelPP;
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauTrue")));
	bool found =  true; bool skipFit = false;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, false, false, fitCtauTrue, false, false, cutSideBand);
	if (found) {
	  if (fitCtau) {
	    cout << "[INFO] This ctauTrue fit was already done, so I'll load the fit results." << endl;
	    if ( loadPreviousFitResult(myws, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, false, false, fitCtauTrue, false, false, false) ) {
	      skipFit = true;
	    }
	  } else {
	    cout << "[INFO] This ctauTrue fit was already done, so I'll just go to the next one." << endl;
	    return true;
	  }
	}
	if (skipFit==false) {
	  // Fit the Datasets
	  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("FullWindow"), NumCPU(numCores), Save());
	  fitResult->Print();
	  // Draw the mass plot
	  int nBins = min(int( round((cut.dMuon.ctauTrue.Max - cut.dMuon.ctauTrue.Min)/binWidth) ), 1000);
	  drawCtauTruePlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, false, incJpsi, incPsi2S, cutCtau, wantPureSMC, setLogScale, incSS, nBins); 
	  saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, false, false, fitCtauTrue, false, false, cutSideBand);
	  myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
	}
	if (fitCtau) {
	  myws.pdf("pdfCTAUTRUE_Tot_PP")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
	}
      }
      if (fitMass && !fitCtau) {
	pdfName = "pdfMASS_Tot_PP"; plotLabel = plotMassLabelPP;
        
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass")));
	bool found =  true;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, fitMass, false, false, false, false, cutSideBand);
	if (found) {
	  cout << "[INFO] This mass fit was already done, so I'll just go to the next one." << endl;
	  return true;
	}
      
	RooFitResult* fitResult = myws.pdf("pdfMASS_Tot_PP")->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE),
								    NumCPU(numCores), Range((incJpsi||incPsi2S)?"FullWindow":parIni["BkgMassRange_FULL_Label"].c_str()), Save());
	fitResult->Print();
	
	// Draw the mass plot
	int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
	drawMassPlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, false, incJpsi, incPsi2S, incBkg, cutCtau, wantPureSMC, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
	
	if (fitMass && incBkg ){
	vector<string> rangeColl; vector<string> cutColl; vector<string> objectColl;
	rangeColl.push_back("JpsiWindow"); cutColl.push_back(parIni["JpsiMassRange_Cut"]);  objectColl.push_back("Jpsi");
	rangeColl.push_back(parIni["BkgMassRange_JPSI_Label"]); cutColl.push_back(parIni["BkgMassRange_JPSI_Cut"]);   objectColl.push_back("Bkg");
	//if (incPsi2S) { rangeColl.push_back("Psi2SWindow");                cutColl.push_back(parIni["Psi2SMassRange_Cut"]); objectColl.push_back("Psi2S"); }
	if( !makeCtauErrPdf(myws, dsName, incPsi2S, isPbp, fitMass) ) { return false; }
	}
	
	saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, fitMass, false, false, false, false, cutSideBand);
	myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
      }

      if (fitCtau) { 
	if (!fitMass && fitCtau) { pdfName = "pdfCTAU_Tot_PP";     plotLabel = plotCtauLabelPP;                 }
	if (fitMass &&  fitCtau) { pdfName = "pdfCTAUMASS_Tot_PP"; plotLabel = plotMassLabelPP+plotCtauLabelPP; }
	cout << "====================ok1========================================" << endl;
        
	if (fitMass) {
	  string plotMassLabel = ""; 
	  plotMassLabel = plotMassLabel + Form("_Jpsi_%s", parIni["Model_Jpsi_PP"].c_str());
	  plotMassLabel = plotMassLabel + Form("_Psi2S_%s", parIni["Model_Psi2S_PP"].c_str());
	  plotMassLabel = plotMassLabel + Form("_Bkg_%s", parIni["Model_Bkg_PP"].c_str());
	  loadPreviousFitResult(myws, outputDir, plotMassLabel, DSTAG, cut, true, false, false, false, false, cutSideBand);
	  cout << "====================ok2========================================" << endl;
	  myws.var("invMass")->setRange("InclusiveMassRange",cut.dMuon.M.Min, cut.dMuon.M.Max);
	  if (myws.pdf("pdfMASS_Tot_PP")) {
	    myws.pdf("pdfMASS_Tot_PP")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE); 
	  }
	  if (myws.pdf("pdfMASS_Bkg_PP"))   { 
	    ((RooChebychev*)myws.pdf("pdfMASS_Bkg_PP"))->selectNormalizationRange("InclusiveMassRange", kTRUE);   
	    double NormFactor = 1.0/(myws.pdf("pdfMASS_Bkg_PP")->createIntegral(*myws.var("invMass"), NormSet(*myws.var("invMass")), Range("InclusiveMassRange"))->getValV());
	    myws.var("N_Bkg_PP")->setVal(myws.var("N_Bkg_PP")->getValV()*NormFactor); 
	  } 
	  if (myws.var("N_Bkg_PP")) { 
	    myws.var("N_Bkg_PP")->setConstant(kFALSE); 
	  }
	  if (myws.pdf("pdfMASS_Jpsi_PP"))   { 
	    double NormFactor = 1.0/(myws.pdf("pdfMASS_Jpsi_PP")->createIntegral(*myws.var("invMass"), NormSet(*myws.var("invMass")), Range("InclusiveMassRange"))->getValV());
	    myws.var("N_Jpsi_PP")->setVal(myws.var("N_Jpsi_PP")->getValV()*NormFactor); 
	  } 
	  if (myws.var("N_Jpsi_PP")) { 
	    myws.var("N_Jpsi_PP")->setConstant(kFALSE); 
	  }
	  if (myws.pdf("pdfMASS_Psi2S_PP"))   { 
	    double NormFactor = 1.0/(myws.pdf("pdfMASS_Psi2S_PP")->createIntegral(*myws.var("invMass"), NormSet(*myws.var("invMass")), Range("InclusiveMassRange"))->getValV());
	    myws.var("N_Psi2S_PP")->setVal(myws.var("N_Psi2S_PP")->getValV()*NormFactor); 
	  } 
	  if (myws.var("N_Psi2S_PP")) { 
	    myws.var("N_Psi2S_PP")->setConstant(kFALSE); 
	  }
	}
	
	cout << "====================ok3========================================" << endl;
	double numEntries = myws.data(dsName.c_str())->sumEntries();
	buildCharmoniaCtauModel(myws, model.PP, parIni, dsName, false, incBkg, incJpsi, incPsi2S, incPrompt, incNonPrompt, isErrPdfData, numEntries);      
	//{ return false; }
	cout << "====================ok4========================================" << endl;
	
	// check if we have already done this fit. If yes, do nothing and return true.
	RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctau")));
	bool found =  true;
	found = found && isFitAlreadyFound(newpars, outputDir, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, cut, fitMass, fitCtau, false, false, false, cutSideBand);
	if (found) {
	  cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
	  return true;
	}
	
	if (!isMC){
	  if (incPrompt) {
	    string plotCtauLabel = Form("_JpsiPR_%s", (parIni["Model_JpsiPR_PP"]==""?"Delta":parIni["Model_JpsiPR_PP"].c_str()));
	    loadPreviousFitResult(myws, outputDir, plotCtauLabel, "MCJPSIPR", cut, false, true, false, false, false, false);
	    if (myws.pdf("pdfCTAUCOND_JpsiPR_PP")) { 
	      myws.pdf("pdfCTAUCOND_JpsiPR_PP")->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	    }
	    if (myws.pdf("pdfCTAUCOND_Psi2SPR_PP")) { 
	      myws.pdf("pdfCTAUCOND_Psi2SPR_PP")->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	    }
	    if (myws.pdf("pdfCTAUCOND_BkgPR_PP")) { 
	      myws.pdf("pdfCTAUCOND_BkgPR_PP")->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	    }
	  }
	  
	  if (incNonPrompt && (incJpsi || incPsi2S)) {
	    string plotCtauTrueLabel = string("") + Form("_CtauTrueRes_%s", parIni["Model_CtauTrueRes_PP"].c_str()) + Form("_CtauTrue_%s", parIni["Model_CtauTrue_PP"].c_str());
	    loadPreviousFitResult(myws, outputDir, plotCtauTrueLabel, "MCJPSINOPR", cut, false, false, true, false, false, false);
	    if ( myws.var("sigmaMC_JpsiNoPR_PP") )    { myws.var("sigmaMC_JpsiNoPR_PP")->setConstant(kTRUE);  }
	    if ( myws.var("sigmaMC_Psi2SNoPR_PP") )   { myws.var("sigmaMC_Psi2SNoPR_PP")->setConstant(kTRUE); }
            
	    string plotCtauLabel = Form("_BkgPR_%s", parIni["Model_BkgPR_PP"].c_str());
	    plotCtauLabel = plotCtauLabel + Form("_BkgNoPR_%s", parIni["Model_BkgNoPR_PP"].c_str());
	    loadPreviousFitResult(myws, outputDir, plotCtauLabel, "DATA", cut, false, true, false, false, false, true);
	    if (myws.pdf("pdfCTAUCOND_BkgNoPR_PP")) { 
	      myws.pdf("pdfCTAUCOND_BkgNoPR_PP")->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
	    }
	    if (myws.var("b_Bkg_PP")) { myws.var("b_Bkg_PP")->setConstant(kTRUE); }
	  }
	}
	if (isMC) {
	  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores), ConditionalObservables(*myws.var("ctauErr")), Save());
	  fitResult->Print();
	} else {
	  if (fitMass){
	    cout << "[INFO] Fitting pdf: " << pdfName << endl;
	    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores), Save(), ConditionalObservables(*myws.var("ctauErr")));
	    fitResult->Print("v");
	  } else {
	    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), NumCPU(numCores),  Save(), ConditionalObservables(*myws.var("ctauErr")));
	    fitResult->Print();
	  }
	}
      
	//==============================================================================
        cout << "====================okA========================================" << endl;
        
	int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
	drawCtauMassPlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, false, fitMass, incJpsi, incPsi2S, incBkg, cutCtau, wantPureSMC, setLogScale, incSS, zoomPsi, nBins); 
	//==================================================================================
	
	cout << "====================ok5========================================" << endl;
	nBins = min(int( round((cut.dMuon.ctau.Max - cut.dMuon.ctau.Min)/binWidth) ), 1000);
	drawCtauPlot(myws, outputDir, opt, cut, (wantPureSMC ? (plotLabel+"_NoBkg") : plotLabel), DSTAG, false, fitMass, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, cutCtau, wantPureSMC, setLogScale, incSS, nBins);
	cout << "====================ok6========================================" << endl;
	double width = (myws.var("ctauErr")->getMax() - myws.var("ctauErr")->getMin()) / 100.0;
	nBins = min(int( round( (cut.dMuon.ctauErr.Max - cut.dMuon.ctauErr.Min)/width )), 1000);
	if (!isMC) drawCtauErrorPlot(myws, outputDir, opt, cut, DSTAG, false, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, cutCtau, wantPureSMC, setLogScale, incSS, nBins);
	saveWorkSpace(myws, outputDir, plotLabel, DSTAG, cut, fitMass, fitCtau, false, false, false, cutSideBand);
	myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE) ;
      }
    }   
  }
  return true;
};


void setCtauCuts(struct KinCuts& cut, bool isPbp) 
{
  if (cut.dMuon.AbsRap.Max<=1.6 && isPbp) {
    cut.dMuon.ctauCut = "( ctau < (0.013 + (0.22/pt)) )";
  }
  if (cut.dMuon.AbsRap.Min>=1.6 && isPbp) {
    cut.dMuon.ctauCut = "( ctau < (0.015 + (0.28/pt)) )";
  }
  if (cut.dMuon.AbsRap.Max<=1.6 && !isPbp) {
    cut.dMuon.ctauCut = "( ctau < (0.010 + (0.25/pt)) )";
  }
  if (cut.dMuon.AbsRap.Min>=1.6 && !isPbp) {
    cut.dMuon.ctauCut = "( ctau < (0.013 + (0.29/pt)) )";
  }
};


bool setMassModel( struct OniaModel& model, map<string, string>& parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg )
{
  if (isPbp && incBkg) {
    if (parIni.count("Model_Bkg_Pbp")>0) {
      model.Pbp.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_Pbp"]];
      if (model.Pbp.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incJpsi) {
    if (parIni.count("Model_Jpsi_Pbp")>0) {
      model.Pbp.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_Pbp"]];
      if (model.Pbp.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incPsi2S) {
    if (parIni.count("Model_Psi2S_Pbp")>0) {
      model.Pbp.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_Pbp"]];
      if (model.Pbp.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incBkg) {
    if (parIni.count("Model_Bkg_PP")>0) {
      model.PP.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PP"]];
      if (model.PP.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incJpsi) {
    if (parIni.count("Model_Jpsi_PP")>0) {
      model.PP.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PP"]];
      if (model.PP.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incPsi2S) {
    if (parIni.count("Model_Psi2S_PP")>0) {
      model.PP.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PP"]];
      if (model.PP.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};


bool setCtauModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt )
{
  if (isPbp) {
    if (parIni.count("Model_CtauRes_Pbp")>0) {
      model.Pbp.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_Pbp"]];
      if (model.Pbp.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  } else {
    if (parIni.count("Model_CtauRes_PP")>0) {
      model.PP.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PP"]];
      if (model.PP.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incBkg && incNonPrompt) {
    if (parIni.count("Model_BkgNoPR_Pbp")>0) {
      model.Pbp.Bkg.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_BkgNoPR_Pbp"]];
      if (model.Pbp.Bkg.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The background non-prompt ctau model: " << parIni["Model_BkgNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incBkg && incPrompt) {
    if (parIni.count("Model_BkgPR_Pbp")>0) {
      model.Pbp.Bkg.Ctau.Prompt = CtauModelDictionary[parIni["Model_BkgPR_Pbp"]];
      if (model.Pbp.Bkg.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The background prompt ctau model: " << parIni["Model_BkgPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else {
      parIni["Model_BkgPR_Pbp"] = "Delta";
      model.Pbp.Bkg.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (isPbp && incJpsi && incNonPrompt) {
    if (parIni.count("Model_JpsiNoPR_Pbp")>0) {
      model.Pbp.Jpsi.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_JpsiNoPR_Pbp"]];
      if (model.Pbp.Jpsi.Ctau.NonPrompt==CtauModel(0)) {
	cout << "[ERROR] The Jpsi non-prompt ctau model: " << parIni["Model_JpsiNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incJpsi && incPrompt) {
    if (parIni.count("Model_JpsiPR_Pbp")>0) {
      model.Pbp.Jpsi.Ctau.Prompt = CtauModelDictionary[parIni["Model_JpsiPR_Pbp"]];
      if (model.Pbp.Jpsi.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi prompt ctau model: " << parIni["Model_JpsiPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_JpsiPR_Pbp"] = "Delta";
      model.Pbp.Jpsi.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (isPbp && incPsi2S && incNonPrompt) {
    if (parIni.count("Model_Psi2SNoPR_Pbp")>0) {
      model.Pbp.Psi2S.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_Psi2SNoPR_Pbp"]];
      if (model.Pbp.Psi2S.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) non-prompt ctau model: " << parIni["Model_Psi2SNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incPsi2S && incPrompt) {
    if (parIni.count("Model_Psi2SPR_Pbp")>0) {
      model.Pbp.Psi2S.Ctau.Prompt = CtauModelDictionary[parIni["Model_Psi2SPR_Pbp"]];
      if (model.Pbp.Psi2S.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) prompt ctau model: " << parIni["Model_Psi2SPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_Psi2SPR_Pbp"] = "Delta";
      model.Pbp.Psi2S.Ctau.Prompt=CtauModel::Delta;
    }
  }

  if (!isPbp && incBkg && incNonPrompt) {
    if (parIni.count("Model_BkgNoPR_PP")>0) {
      model.PP.Bkg.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_BkgNoPR_PP"]];
      if (model.PP.Bkg.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The background non-prompt ctau model: " << parIni["Model_BkgNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incBkg && incPrompt) {
    if (parIni.count("Model_BkgPR_PP")>0) {
      model.PP.Bkg.Ctau.Prompt = CtauModelDictionary[parIni["Model_BkgPR_PP"]];
      if (model.PP.Bkg.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The background prompt ctau model: " << parIni["Model_BkgPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_BkgPR_PP"] = "Delta";
      model.PP.Bkg.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (!isPbp && incJpsi && incNonPrompt) {
    if (parIni.count("Model_JpsiNoPR_PP")>0) {
      model.PP.Jpsi.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_JpsiNoPR_PP"]];
      if (model.PP.Jpsi.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi non-prompt ctau model: " << parIni["Model_JpsiNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incJpsi && incPrompt) {
    if (parIni.count("Model_JpsiPR_PP")>0) {
      model.PP.Jpsi.Ctau.Prompt = CtauModelDictionary[parIni["Model_JpsiPR_PP"]];
      if (model.PP.Jpsi.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi prompt ctau model: " << parIni["Model_JpsiPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_JpsiPR_PP"] = "Delta";
      model.PP.Jpsi.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (!isPbp && incPsi2S && incNonPrompt) {
    if (parIni.count("Model_Psi2SNoPR_PP")>0) {
      model.PP.Psi2S.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_Psi2SNoPR_PP"]];
      if (model.PP.Psi2S.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) non-prompt ctau model: " << parIni["Model_Psi2SNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incPsi2S && incPrompt) {
    if (parIni.count("Model_Psi2SPR_PP")>0) {
      model.PP.Psi2S.Ctau.Prompt = CtauModelDictionary[parIni["Model_Psi2SPR_PP"]];
      if (model.PP.Psi2S.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) prompt ctau model: " << parIni["Model_Psi2SPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_Psi2SPR_PP"] = "Delta";
      model.PP.Psi2S.Ctau.Prompt=CtauModel::Delta;
    }
  }

  return true;
};


bool setCtauTrueModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp )
{
  if (isPbp) {
    if (parIni.count("Model_CtauTrueRes_Pbp")>0) {
      model.Pbp.CtauTrueRes = CtauModelDictionary[parIni["Model_CtauTrueRes_Pbp"]];
      if (model.Pbp.CtauTrueRes==CtauModel(0)) {
        cout << "[ERROR] The ctau truth resolution model: " << parIni["Model_CtauTrueRes_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Truth resolution model for Pbp was not found in the initial parameters!" << endl; return false;
    }
    if (parIni.count("Model_CtauTrue_Pbp")>0) {
      model.Pbp.CtauTrue = CtauModelDictionary[parIni["Model_CtauTrue_Pbp"]];
      if (model.Pbp.CtauTrue==CtauModel(0)) {
        cout << "[ERROR] The ctau truth model: " << parIni["Model_CtauTrue_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Truth model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  } else {
    if (parIni.count("Model_CtauTrueRes_PP")>0) {
      model.PP.CtauTrueRes = CtauModelDictionary[parIni["Model_CtauTrueRes_PP"]];
      if (model.PP.CtauTrueRes==CtauModel(0)) {
        cout << "[ERROR] The ctau truth resolution model: " << parIni["Model_CtauTrue_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau truth resolution model for PP was not found in the initial parameters!" << endl; return false;
    }
    if (parIni.count("Model_CtauTrue_PP")>0) {
      model.PP.CtauTrue = CtauModelDictionary[parIni["Model_CtauTrue_PP"]];
      if (model.PP.CtauTrue==CtauModel(0)) {
        cout << "[ERROR] The ctau truth model: " << parIni["Model_CtauTrue_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau truth model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};

   
int importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts& cut, string label, bool incJpsi, bool incPsi2S, bool cutSideBand)
{
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  if (cutSideBand) {
    //indMuonMass =  indMuonMass + "&&" + "((2.0 < invMass && invMass < 2.9) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
    //indMuonMass =  "((2.6 <= invMass && invMass <= 2.9) || (3.3 < invMass && invMass < 3.5))";
    indMuonMass =  "((3.4 <= invMass && invMass <= 3.55) || (3.85 < invMass && invMass < 4.2) )";
  }
  
  //string indMuonRap     = Form("(%.6f <= abs(rap) && abs(rap) < %.6f)",    cut.dMuon.AbsRap.Min,  cut.dMuon.AbsRap.Max);
  //if(runNb >= 211313 && runNb <= 211631){
  //double rapvalHi = (-1)*cut.dMuon.AbsRap.Min;
    //double rapvalLow = (-1)*cut.dMuon.AbsRap.Max;
    //string indMuonRap     = Form("(%.6f <= rap && rap < %.6f)",    rapvalLow,  rapvalHi);
  //}
  
  string indMuonRap     = Form("(%.6f <= rap && rap < %.6f)",    cut.dMuon.AbsRap.Min,  cut.dMuon.AbsRap.Max);
  string indMuonPt      = Form("(%.6f <= pt && pt < %.6f)",                cut.dMuon.Pt.Min,      cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.12f <= ctau && ctau < %.12f)",          cut.dMuon.ctau.Min,    cut.dMuon.ctau.Max); 
  
  string indMuonCtauErr = Form("(%.12f <= ctauErr && ctauErr < %.12f)",    cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  string indMuonCtauErrMC = Form("(%.12f <= ctauErr && ctauErr < %.12f)",  0.005, 0.3);
  string inCentrality   = Form("(%d <= cent && cent < %d)",                cut.Centrality.Start,  cut.Centrality.End);
  string indMuonCtauTrue = Form("(%.12f <= ctauTrue && ctauTrue < %.12f)", cut.dMuon.ctauTrue.Min,    cut.dMuon.ctauTrue.Max);

  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr;
  string strCutMC         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt+"&&"+ indMuonCtau +"&&"+ indMuonCtauErrMC;
  //if (label.find("Pbp")!=std::string::npos){ strCut = strCut +"&&"+ inCentrality; }
  if (label.find("MC")!=std::string::npos && label.find("PR")!=std::string::npos){ strCut =  strCutMC;} 
  if (label.find("MC")!=std::string::npos && label.find("NOPR")!=std::string::npos){ strCut = strCutMC +"&&"+ indMuonCtauTrue; } 

  //Double_t ctauErrMax; Double_t ctauErrMin;
  //dataOS->getRange("ctauErr", ctauErrMin, ctauErrMax);
  ////string indMuonCtauErr = Form("(%.12f <= ctauErr && ctauErr < %.12f)",    cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  //string indMuonCtauErr = Form("(%.12f <= ctauErr && ctauErr < %.12f)",    ctauErrMin, ctauErrMax);
  
  // Reduce and import the datasets
  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return -1;
  }
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  cout<<"---------------------------------------------------------------"<<endl<<endl;
  cout<<" no of entries in Roodata : "<<dataOS->sumEntries()<<endl<<endl;
  cout<<"---------------------------------------------------------------"<<endl<<endl;

  if (dataOS->sumEntries()==0){ 
    cout << "[WARNING] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return 0;
  }
  myws.import(*dataOS);
  delete dataOS;
  
  if (label.find("NoBkg")==std::string::npos) // Don't try to find SS dataset if label contais NoBkg
  {
    if (!(inputWS.data(Form("dSS_%s", label.c_str())))){
      cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
      return -1;
    }
    RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
    if (dataSS->sumEntries()==0){
      cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    }
    myws.import(*dataSS);
    delete dataSS;
  }
  
   cout<<"------------------------------ok1---------------------------------"<<endl<<endl;
  const RooArgSet* rowOS = myws.data(Form("dOS_%s", label.c_str()))->get();
  ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);        
  ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
  ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);            
  ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
  ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);        
  ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);  
  ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
  //if (label.find("PbPb")!=std::string::npos){
  //((RooRealVar*)rowOS->find("cent"))->setMin(cut.Centrality.Start);      
  //((RooRealVar*)rowOS->find("cent"))->setMax(cut.Centrality.End);
  //}
  if (label.find("MC")!=std::string::npos && label.find("NOPR")!=std::string::npos){ 
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);      
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
    }
  
  cout<<"------------------------------ok2---------------------------------"<<endl<<endl;
  
  // Set the range of each global parameter in the local workspace
  myws.var("invMass")->setMin(cut.dMuon.M.Min);        
  myws.var("invMass")->setMax(cut.dMuon.M.Max);
  myws.var("pt")->setMin(cut.dMuon.Pt.Min);            
  myws.var("pt")->setMax(cut.dMuon.Pt.Max);
  myws.var("rap")->setMin(cut.dMuon.AbsRap.Min);       
  myws.var("rap")->setMax(cut.dMuon.AbsRap.Max);
  myws.var("ctau")->setMin(cut.dMuon.ctau.Min);        
  myws.var("ctau")->setMax(cut.dMuon.ctau.Max);
  myws.var("ctauErr")->setMin(cut.dMuon.ctauErr.Min);  
  myws.var("ctauErr")->setMax(cut.dMuon.ctauErr.Max);
  //if (label.find("PbPb")!=std::string::npos){
  //myws.var("cent")->setMin(cut.Centrality.Start);      
  //myws.var("cent")->setMax(cut.Centrality.End);
  //}
  if (label.find("MC")!=std::string::npos && label.find("NOPR")!=std::string::npos){ 
    myws.var("ctauTrue")->setMin(cut.dMuon.ctauTrue.Min);      
    myws.var("ctauTrue")->setMax(cut.dMuon.ctauTrue.Max);
  }
  cout << "[INFO] Analyzing bin: " << Form(
                                           "%.3f < pt < %.3f, %.3f < rap < %.3f", 
                                           cut.dMuon.Pt.Min,
                                           cut.dMuon.Pt.Max,
                                           cut.dMuon.AbsRap.Min,
                                           cut.dMuon.AbsRap.Max
					   ) << endl;
   cout<<"-------------------------ok3--------------------------------------"<<endl<<endl;

  return 1;

};

void setGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, bool fitCtau, bool fitCtauTrue, bool incJpsi, bool incPsi2S, bool incBkg)
{
  if (fitCtau) {
    Double_t ctauErrMax; Double_t ctauErrMin;
    myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctauErr"), ctauErrMin, ctauErrMax);
    //if (cut.dMuon.ctauErr.Min<ctauErrMin && cut.dMuon.ctauErr.Max>ctauErrMin) { cut.dMuon.ctauErr.Min = ctauErrMin; }
    //if (cut.dMuon.ctauErr.Max>ctauErrMax && cut.dMuon.ctauErr.Min<ctauErrMax) { cut.dMuon.ctauErr.Max = ctauErrMax; }
    cout << "Range from data: ctauErrMin: " << ctauErrMin << "  ctauErrMax: " << ctauErrMax << endl;
    cout << "Range from file: ctauErrMin: " << cut.dMuon.ctauErr.Min << "  ctauErrMax: " <<cut.dMuon.ctauErr.Max << endl;

    if (label.find("MC")!=std::string::npos){
      ctauErrMin = cut.dMuon.ctauErr.Min;
      ctauErrMax = cut.dMuon.ctauErr.Max;
    }
    
    myws.var("ctauErr")->setMin(ctauErrMin);  
    myws.var("ctauErr")->setMax(ctauErrMax);
    myws.var("ctauErr")->setRange("FullWindow",  ctauErrMin, ctauErrMax);
   
    if (incBkg) { 
      myws.var("ctauErr")->setRange("SideBandTOP_FULL",   ctauErrMin, ctauErrMax); 
      myws.var("ctauErr")->setRange("SideBandMID_FULL",   ctauErrMin, ctauErrMax);
      myws.var("ctauErr")->setRange("SideBandBOT_FULL",   ctauErrMin, ctauErrMax); 
      if (incJpsi)  { 
        myws.var("ctauErr")->setRange("SideBandMID_JPSI",   ctauErrMin, ctauErrMax);
        myws.var("ctauErr")->setRange("SideBandBOT_JPSI",   ctauErrMin, ctauErrMax);
      }
      if (incPsi2S) { 
        myws.var("ctauErr")->setRange("SideBandTOP_PSI2S",  ctauErrMin, ctauErrMax); 
        myws.var("ctauErr")->setRange("SideBandMID_PSI2S",  ctauErrMin, ctauErrMax);
      }
    }
    if (incJpsi)  { myws.var("ctauErr")->setRange("JpsiWindow",  ctauErrMin, ctauErrMax); }
    if (incPsi2S) { myws.var("ctauErr")->setRange("Psi2SWindow", ctauErrMin, ctauErrMax); }
    parIni["CtauErrRange_Cut"]   = Form("(%.12f <= ctauErr && ctauErr < %.12f)", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);; 

    Double_t ctauMax=cut.dMuon.ctau.Max; Double_t ctauMin=cut.dMuon.ctau.Min;
    cout << "Range from input: ctauMin: " << ctauMin << "  ctauMax: " << ctauMax << endl;
    myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctau"), ctauMin, ctauMax);
    //if (ctauMin<cut.dMuon.ctau.Min) { ctauMin = cut.dMuon.ctau.Min; }
    //if (ctauMax>cut.dMuon.ctau.Max) { ctauMax = cut.dMuon.ctau.Max; }
    //TH1* h = myws.data(Form("dOS_%s", label.c_str()))->createHistogram("hist", *myws.var("ctau"));
    //ctauMin = h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove(1.0));
    //ctauMax = h->GetXaxis()->GetBinUpEdge(h->FindLastBinAbove(1.0));
    //delete h;
    //myws.var("ctau")->setMin(ctauMin);  
    //myws.var("ctau")->setMax(ctauMax);

    //ctauMax=cut.dMuon.ctau.Max;
    //ctauMin=cut.dMuon.ctau.Min;
    
    cout << "Range from data: ctauMin: " << ctauMin << "  ctauMax: " << ctauMax << endl;
    myws.var("ctau")->setRange("FullWindow", ctauMin, ctauMax);
    if (incBkg) { 
      myws.var("ctau")->setRange("SideBandTOP_FULL",  ctauMin, ctauMax); 
      myws.var("ctau")->setRange("SideBandMID_FULL",  ctauMin, ctauMax);
      myws.var("ctau")->setRange("SideBandBOT_FULL",  ctauMin, ctauMax); 
      if (incJpsi)  { 
        myws.var("ctau")->setRange("SideBandMID_JPSI",  ctauMin, ctauMax);
        myws.var("ctau")->setRange("SideBandBOT_JPSI",  ctauMin, ctauMax);
      }
      if (incPsi2S) { 
        myws.var("ctau")->setRange("SideBandTOP_PSI2S", ctauMin, ctauMax); 
        myws.var("ctau")->setRange("SideBandMID_PSI2S", ctauMin, ctauMax);
      }
    }
    if (incJpsi)  { myws.var("ctau")->setRange("JpsiWindow", ctauMin, ctauMax);  }
    if (incPsi2S) { myws.var("ctau")->setRange("Psi2SWindow", ctauMin, ctauMax); }
    parIni["CtauRange_Cut"]   = Form("(%.12f <= ctau && ctau < %.12f)", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  }

  if (fitCtauTrue){ 
    Double_t ctauTrueMax; Double_t ctauTrueMin;
    myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctauTrue"), ctauTrueMin, ctauTrueMax);
    if (ctauTrueMin<cut.dMuon.ctauTrue.Min) { ctauTrueMin = cut.dMuon.ctauTrue.Min; }
    if (ctauTrueMax>cut.dMuon.ctauTrue.Max) { ctauTrueMax = cut.dMuon.ctauTrue.Max; }
    cout << "Range from data: ctauTrueMin: " << ctauTrueMin << "  ctauTrueMax: " << ctauTrueMax << endl;
    myws.var("ctauTrue")->setRange("CtauTrueWindow", ctauTrueMin, ctauTrueMax);
    parIni["CtauTrueRange_Cut"]   = Form("(%.12f <= ctauTrue && ctauTrue < %.12f)", ctauTrueMin, ctauTrueMax);
  }

  if (label.find("MC")!=std::string::npos)
    {
      myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
      parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, cut.dMuon.M.Max);
      /*
      if (incPsi2S)
        {
          if (cut.dMuon.AbsRap.Min >= 1.6) {
            myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, 3.95);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.95);
          }
          else { 
            myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, 3.85);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.85);
          }
        }
      if (incJpsi)
        {
          if (cut.dMuon.AbsRap.Min >= 1.6) {
            myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, 3.32);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.32);
          }
          else {
            myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, 3.26);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.26);
          }
	  }*/
    }
  else {
    myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
    parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, cut.dMuon.M.Max);
    if (incJpsi) {
      myws.var("invMass")->setRange("JpsiWindow", 2.9, 3.3);
      parIni["JpsiMassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 2.9, 3.3);
    }
    if (incPsi2S) {
      myws.var("invMass")->setRange("Psi2SWindow", 3.55, 3.85);
      parIni["Psi2SMassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 3.55, 3.85);
    }
    if (incBkg) {
      //myws.var("invMass")->setRange("SideBandMID_FULL",  cut.dMuon.M.Min, cut.dMuon.M.Max);
      myws.var("invMass")->setRange("SideBandMID_FULL",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.5)?3.5:cut.dMuon.M.Max));
      myws.var("invMass")->setRange("SideBandMID_JPSI",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.4)?3.4:cut.dMuon.M.Max));
      myws.var("invMass")->setRange("SideBandMID_PSI2S", ((cut.dMuon.M.Min<3.4)?3.4:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.5)?3.55:cut.dMuon.M.Max));
      parIni["BkgMassRange_FULL_Label"]  = "SideBandMID_FULL";
      parIni["BkgMassRange_JPSI_Label"]  = "SideBandMID_JPSI";
      parIni["BkgMassRange_PSI2S_Label"] = "SideBandMID_PSI2S";
      if (cut.dMuon.M.Min < 2.9) {
        myws.var("invMass")->setRange("SideBandBOT_FULL", cut.dMuon.M.Min, 2.9);
        myws.var("invMass")->setRange("SideBandBOT_JPSI", ((cut.dMuon.M.Min<2.5)?2.5:cut.dMuon.M.Min), 2.9);
        parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandBOT_FULL";
        parIni["BkgMassRange_JPSI_Label"] = parIni["BkgMassRange_JPSI_Label"] + "," + "SideBandBOT_JPSI";
      }
      if (cut.dMuon.M.Max > 3.9) {
        myws.var("invMass")->setRange("SideBandTOP_FULL", 3.85, cut.dMuon.M.Max);
        myws.var("invMass")->setRange("SideBandTOP_PSI2S", 3.85, ((cut.dMuon.M.Max>4.2)?4.2:cut.dMuon.M.Max));
        parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandTOP_FULL";
        parIni["BkgMassRange_PSI2S_Label"] = parIni["BkgMassRange_PSI2S_Label"] + "," + "SideBandTOP_PSI2S";
      }
      parIni["BkgMassRange_FULL_Cut"]  = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
      //parIni["BkgMassRange_FULL_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "((2.0 < invMass && invMass < 2.9) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
      parIni["BkgMassRange_JPSI_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "((2.5 < invMass && invMass < 3.4))";
      parIni["BkgMassRange_PSI2S_Cut"] = parIni["BkgMassRange_FULL_Cut"] + "&&" + "((3.4 < invMass && invMass < 4.2))";
      parIni["BkgMassRange_FULL_Cut"]  = "("+parIni["BkgMassRange_FULL_Cut"]+")";
      parIni["BkgMassRange_JPSI_Cut"]  = "("+parIni["BkgMassRange_JPSI_Cut"]+")";
      parIni["BkgMassRange_PSI2S_Cut"] = "("+parIni["BkgMassRange_PSI2S_Cut"]+")";
    }
  }

};


void setOptions(struct InputOpt* opt) 
{
  opt->pp.RunNb.Start   = 262157; opt->Pbp.RunNb.Start = 210676;
  opt->pp.RunNb.End     = 262328; opt->Pbp.RunNb.End   = 211631;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; //HLT_PAL1DoubleMuOpen_v1; 
  opt->Pbp.TriggerBit  = (int) HI::HLT_PAL1DoubleMuOpen_v1; //HLT_HIL1DoubleMu0_v1; 
  return;
};


void setNameFile(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut,bool fitMass, bool fitCtau, bool fitCtauTrue, bool isPbp, bool cutSideBand) 
{
  
  if (fitMass || (fitCtau && !fitCtauTrue)) {
    FileName = Form("%sresult/%s/massErrPdf%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"":""), "MASS", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  }
  if (fitCtauTrue) {
    //TAG = "MCJPSINOPR";
    TAG = "MCPSI2SNOPR";
    FileName = Form("%sresult/%s/ctauTrue%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"":""), "CTAUTRUE", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  }
}


void setFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut,bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp, bool cutSideBand) 
{
  if (fitCtauTrue) {
    FileName = Form("%sresult/%s/ctauTrue%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "CTAUTRUE", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  }  else {
    if (doSimulFit) {
      if (fitMass && fitCtau) { 
        FileName = Form("%sresult/%s/ctaumass%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "CTAUMASS", TAG.c_str(), "Psi2SJpsi", "COMB", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      } else if (fitMass) {
        FileName = Form("%sresult/%s/mass%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "MASS", TAG.c_str(), "Psi2SJpsi", "COMB", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      } else if (fitCtau) {
        FileName = Form("%sresult/%s/ctau%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "CTAU", TAG.c_str(), "Psi2SJpsi", "COMB", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      }
    } else {    
      if (fitMass && fitCtau) {
        FileName = Form("%sresult/%s/ctaumass%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "CTAUMASS", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      } else if (fitMass) {
        FileName = Form("%sresult/%s/mass%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "MASS", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      } else if (fitCtau) {
        FileName = Form("%sresult/%s/ctau%s/FIT_%s_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), (cutSideBand?"SB":""), "CTAU", TAG.c_str(), "Psi2SJpsi", (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
      } 
    }
  }
};

bool saveWorkSpace(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp, bool cutSideBand) 
{
  // Save the workspace
  if (fitCtauTrue) { 
    gSystem->mkdir(Form("%sresult/%s/ctauTrue%s", outputDir.c_str(), DSTAG.c_str(), (cutSideBand?"SB":"")), kTRUE); 
  } else {
    if (fitMass && fitCtau) { gSystem->mkdir(Form("%sresult/%s/ctaumass%s", outputDir.c_str(), DSTAG.c_str(), (cutSideBand?"SB":"")), kTRUE); } 
    else if (fitMass) { gSystem->mkdir(Form("%sresult/%s/mass%s", outputDir.c_str(), DSTAG.c_str(), (cutSideBand?"SB":"")), kTRUE); } 
    else if (fitCtau) { gSystem->mkdir(Form("%sresult/%s/ctau%s", outputDir.c_str(), DSTAG.c_str(), (cutSideBand?"SB":"")), kTRUE); }
  }
  string FileName = "";
  setFileName(FileName, outputDir, DSTAG, plotLabel, cut, fitMass, fitCtau, fitCtauTrue, doSimulFit, isPbp, cutSideBand);
  cout << FileName << endl;
  TFile *file =  new TFile(FileName.c_str(), "RECREATE");
  if (!file) { 
    cout << "[ERROR] Output root file with fit results could not be created!" << endl; return false; 
  } else {
    file->cd();    
    myws.Write("workspace"); 
    file->Write(); file->Close(); delete file;
  }
  return true;
};


bool loadFinalFitResult(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool incJpsi, bool incPsi2S, bool isPbp, bool cutSideBand){
  string FileName = "";

  //if(fitCtau && fitMass)
  setNameFile(FileName, outputDir, DSTAG, plotLabel, cut, fitMass, fitCtau, fitCtauTrue, isPbp, cutSideBand);
  //if(fitCtauTrue)
  //setFileName(FileName, outputDir, DSTAG, plotLabel, cut, fitMass, fitCtau, fitCtauTrue, doSimulFit, isPbp, cutSideBand);
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
  
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    file->Close(); delete file;
    return false;
  }

  //ws->Print("v");
  
  string Obj;
  if(incJpsi) { Obj = "Jpsi";}
  if(incPsi2S) {
    Obj = "Psi2S";
    cout << "Filename: " << FileName<<"  "<<Obj << endl;
  }
  
  cout<<"=====================okA1====================="<<endl;
  RooArgSet*  varList1 = NULL;
  RooArgSet*  varList2 = NULL;
  
  if(fitCtau &&!fitMass && !cutSideBand){
    varList1 = new RooArgSet(*(ws->pdf(Form("pdfCTAUERR_%s_%s", Obj.c_str() , (isPbp?"Pbp":"PP")))) );
    myws.import(*ws->data("binDataCtErr")); myws.import(*ws->data("binDataCtErrSIG")); myws.import(*ws->data("binSubtractedSIG")); myws.import(*ws->data(Form("hCtauErr_Tot_%s", (isPbp?"Pbp":"PP"))));
    myws.import(*varList1);
  }
  
  //if(cutSideBand && !fitMass){
  //varList1 =   new RooArgSet(*(ws->pdf(Form("pdfCTAU_Bkg_%s", (isPbp?"Pbp":"PP")))), *(ws->pdf(Form("pdfCTAURES_BkgPR_%s", (isPbp?"Pbp":"PP")))), *(ws->pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP")))) );
  //varList1 =   new RooArgSet(*(ws->pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP")))) );  
  //myws.import(*ws->data("binDataCtErr")); myws.import(*ws->data("binDataCtErrSIG")); myws.import(*ws->data("binSubtractedSIG")); myws.import(*ws->data(Form("hCtauErr_Tot_%s", (isPbp?"Pbp":"PP")))); myws.import(*ws->data("binDataCtErrSB")); myws.import(*ws->data("binScaleBKG"));
  //myws.import(*varList1);
  //}
  
  if(fitCtauTrue){                             
    if(incJpsi)  varList1 =   new RooArgSet(*(ws->var(Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))) );
    if(incPsi2S) varList1 =   new RooArgSet(*(ws->var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))) );
    myws.import(*varList1);
  }
  
  //if(fitCtau && fitMass) {
  //varList2 = new RooArgSet(*(ws->pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP")))), *(ws->pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP")))) );
    //myws.import(*ws->data("binDataCtErrSB")); myws.import(*ws->data("binScaleBKG")); 
    //myws.import(*varList2);
  //}
  
  cout<<"=====================okA2====================="<<endl;
  //myws.Print("v");
  
}


bool loadPreviousFitResult(RooWorkspace& myws, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp, bool cutSideBand)
{
  string FileName = "";
  setFileName(FileName, outputDir, DSTAG, plotLabel, cut, fitMass, fitCtau, fitCtauTrue, doSimulFit, isPbp, cutSideBand);
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
 
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    file->Close(); delete file;
    return false;
  }

  
  /*cout << "[INFO] FileName: " << FileName << endl;
    if(ws->pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP")))) {
    myws.import(*ws->pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))));
    myws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
    RooArgSet* model_params = myws.pdf("pdfCTAUCOND_JpsiPR_Pbp")->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"))) ;
    model_params->Print("v");
    cout << " Fixed parameters of pdfCTAUCOND_JpsiPR_Pbp for 2D fit from the" << FileName << endl;
  }
  
  if (ws->pdf("pdfCTAUCOND_Psi2SPR_Pbp")) {
    myws.pdf("pdfCTAUCOND_Psi2SPR_Pbp")->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
    cout << " Fixed parameters of pdfCTAUCOND_Psi2SPR_Pbp for 2D fit from the" << FileName << endl;
    }*/

  //ok244
    
  if(cutSideBand) {
    RooArgSet*  varList3 = NULL;
    RooArgSet*  varList4 = NULL;
    
    varList3 = new RooArgSet(*(ws->var(Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")))), *(ws->var(Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")))) );
    //myws.import(*varList3);
  }

  cout <<  "[INFO] Loading variables and functions from: " << FileName << endl;
  RooArgSet listVar = ws->allVars();
  TIterator* parIt = listVar.createIterator();
  string print = "[INFO] Variables loaded: ";
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    string name = it->GetName();
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || 
         name=="ctauTrue" || name=="pt" || name=="cent" || 
         name=="rap" || name=="One" || name=="Ntracks" || name=="SumET_HFEta4") continue;
    if ( (DSTAG.find("MC")!=std::string::npos || cutSideBand) && (name.find("N_")!=std::string::npos) ) continue; 
    if (myws.var(name.c_str())) { 
      print = print + Form("  %s: %.5f->%.5f  ", name.c_str(), myws.var(name.c_str())->getValV(), ws->var(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->var(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( ws->var(name.c_str())->getError() );
    } else {
      if ( (name==Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))) && myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))) ) {
        print = print + Form("  %s: %.5f->%.5f", Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
      if ( (name==Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))) && myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))) {
        print = print + Form("  %s: %.5f->%.5f  ", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
    }
  }

  cout << print << endl;
  RooArgSet listFun = ws->allFunctions();
  TIterator* parFunIt = listFun.createIterator();
  string printFun = "[INFO] Functions loaded: ";
  for (RooRealVar* it = (RooRealVar*)parFunIt->Next(); it!=NULL; it = (RooRealVar*)parFunIt->Next() ) {
    string name = it->GetName();
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || 
         name=="ctauTrue" || name=="pt" || name=="cent" || 
         name=="rap" || name=="One" || name=="Ntracks" || name=="SumET_HFEta4") continue;
    if ( (DSTAG.find("MC")!=std::string::npos || cutSideBand) && (name.find("N_")!=std::string::npos) ) continue; 
    if (myws.var(name.c_str())) { 
      printFun = printFun + Form("  %s: %.5f->%.5f  ", name.c_str(), myws.var(name.c_str())->getValV(), ws->function(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->function(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( 0.0 );
    }
  }
  cout << printFun << endl;
  return true;
};

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string DSTAG, struct KinCuts cut, bool fitMass, bool fitCtau, bool fitCtauTrue, bool doSimulFit, bool isPbp, bool cutSideBand) 
{
  string FileName = "";
  setFileName(FileName, outputDir, DSTAG, plotLabel, cut, fitMass, fitCtau, fitCtauTrue, doSimulFit, isPbp, cutSideBand);  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "FileName: " << FileName << " was not found" << endl;
    return false; // File was not found
  }

  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    file->Close(); delete file;
    return false;
  }

  string snapShotName = "";
  if (fitCtauTrue) {
    snapShotName = Form("pdfCTAUTRUE_Tot_%s_parIni", (isPbp?"Pbp":"PP"));
  } else {
    if (doSimulFit) {
      string snapShotName = "simPdf_parIni";
    } else {
      if (fitMass && fitCtau) {
        snapShotName = Form("pdfCTAUMASS_Tot_%s_parIni", (isPbp?"Pbp":"PP"));
      } else if (fitMass) {
        snapShotName = Form("pdfMASS_Tot_%s_parIni", (isPbp?"Pbp":"PP"));
      } else if (fitCtau) {
        snapShotName = Form("pdfCTAU_Tot_%s_parIni", (isPbp?"Pbp":"PP"));
      }    
    }
  }

  const RooArgSet *params = ws->getSnapshot(snapShotName.c_str());
  if (!params) {
    delete ws;
    file->Close(); delete file;
    return false;
  }

  bool result = compareSnapshots(newpars, params);

  delete ws;
  file->Close(); delete file; 

  return result;
};

bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2) {
  TIterator* parIt = pars1->createIterator(); 

  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2->getRealValue(it->GetName(),-1e99);
    if ( strcmp(it->GetName(),"ctauErr")==0 || strcmp(it->GetName(),"ctau")==0 || strcmp(it->GetName(),"ctauTrue")==0 ) continue;
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false; // the parameter was found, but with a different value!
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMin() != it->getMin() ) return false; // the parameter has different lower limit
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMax() != it->getMax() ) return false; // the parameter has different upper limit
  }

  return true;
};


bool makeCtauErrPdf(RooWorkspace& ws, string dsName, bool incPsi2S, bool isPbp, bool fitMass) 
{

  //// *** scaleF to scale down ct err dist in 2.9-3.3 GeV/c2
  float bc = -1.0;
  if(ws.var(Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")))){
    if (fitMass) {
      bc = ws.var(Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")))->getVal();
      if(abs(bc) < 0.001) bc = -0.1;
    }
    else bc = -1.0;
  }
  cout << "====================ok22========================================" << endl;
  float scaleF = 0.01;
  if(!incPsi2S) { scaleF = (exp(2.9*bc)-exp(3.3*bc))/(exp(2.6*bc)-exp(2.9*bc)+exp(3.3*bc)-exp(3.5*bc));}
  if(incPsi2S)  { scaleF = (exp(3.55*bc)-exp(3.85*bc))/(exp(3.4*bc)-exp(3.55*bc)+exp(3.85*bc)-exp(4.2*bc));}
  cout<<" scaleF: "<<scaleF<<endl;
  
  RooDataSet *redDataSIG = NULL;
  RooDataSet*redDataSB = NULL;
  
  if(!incPsi2S) {
    redDataSIG = (RooDataSet*)ws.data(dsName.c_str())->reduce("invMass > 2.9 && invMass < 3.3");
    redDataSB = (RooDataSet*)ws.data(dsName.c_str())->reduce("((2.6 < invMass && invMass < 2.9) || (3.3 < invMass && invMass < 3.5) )");
  }

  if(incPsi2S) {
    redDataSIG = (RooDataSet*)ws.data(dsName.c_str())->reduce("invMass > 3.55 && invMass < 3.85");
    redDataSB = (RooDataSet*)ws.data(dsName.c_str())->reduce("((3.4 < invMass && invMass < 3.55) || (3.85 < invMass && invMass < 4.2) )");
  }

  cout<<" redDataSIG->numEntries()  "<<redDataSIG->numEntries()<<" redDataSB "<< redDataSB->numEntries()<<endl; 
  //RooDataSet* redDataSBL = (RooDataSet*)dataOS->reduce("invMass<2.9");
  //RooDataSet* redDataSBR = (RooDataSet*)dataOS->reduce("invMass>3.3");
  //ws.import(*redDataSIG);
  //ws.import(*redDataSB);
  
  cout << "====================ok221========================================" << endl;
  // Extract the ctau error histograms for each mass range 
  RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr","binDataCtErr",RooArgSet(*(ws.var("ctauErr"))), *ws.data(dsName.c_str()) );
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG","Data ct error distribution for sig",RooArgSet(*(ws.var("ctauErr"))), *redDataSIG);
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB","Data ct error distribution for bkg",RooArgSet(*(ws.var("ctauErr"))), *redDataSB);

  cout<<" binDataCtErr->numEntries()  "<<binDataCtErr->numEntries()<<" binDataCtErrSIG "<< binDataCtErrSIG->numEntries()<<" binDataCtErrSB "<<binDataCtErrSB->numEntries()<<endl; 
  
  ws.import(*binDataCtErr);
  ws.import(*binDataCtErrSIG);
  ws.import(*binDataCtErrSB);
  
  cout<<"dsName  "<<dsName<<endl;
  cout << "====================ok23========================================" << endl;
  
  //// **** (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  RooDataHist *binSubtractedSIG, *binScaledBKG;
  binSubtractedSIG = new RooDataHist("binSubtractedSIG","Subtracted data",RooArgSet(*(ws.var("ctauErr")))); 
  binScaledBKG = subtractSidebands(ws,binSubtractedSIG,binDataCtErrSIG,binDataCtErrSB,scaleF,"ctauErr");
  
  cout<<" binScaledBKG  "<<binScaledBKG->numEntries()<<" binSubtractedSIG  "<<binSubtractedSIG->numEntries()<<endl; 
  if (binScaledBKG==NULL ||binSubtractedSIG ==NULL ) { 
    cout << "[ERROR] binScaledBKG " << binScaledBKG << " is NULL!" << endl;
    return false; 
  } 

  ws.import(*binSubtractedSIG);
  ws.import(*binScaledBKG);
  
  cout << "====================ok24========================================" << endl;
  
  // Transform each RooDataHist to a RooHistPDF and import them in the Workspace
  //if (data==NULL || data->sumEntries()==0) { 
  //cout << "[ERROR] DataHist used to create ctau error pdf for " << objectColl.at(j) << " is empty!" << endl; return false; 
  //} 
  string Obj;
  if(incPsi2S) cout<<"Psi2s  "<<endl;
  if(!incPsi2S) {Obj = "Jpsi";}
  if(incPsi2S) {Obj = "Psi2S";}
  
  string pdfNameSIG = Form("pdfCTAUERR_%s_%s", Obj.c_str(), (isPbp?"Pbp":"PP"));
  string pdfNameBkg = Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"));
  RooHistPdf errPdfSig(pdfNameSIG.c_str(), "Error PDF signal", RooArgSet(*(ws.var("ctauErr"))) , *binSubtractedSIG);
  RooHistPdf errPdfBkg(pdfNameBkg.c_str(), "Error PDF Bkg", RooArgSet(*(ws.var("ctauErr"))), *binScaledBKG);
  
  //if ((errPdfSig==NULL) || (errPdfBkg==NULL) ) { 
  //cout << "[ERROR] RooKeysPDF " << pdfNameSIG<<"  "<<pdfNameBkg << " is NULL!" << endl;
  //return false; 
  //} 
  
  cout << "====================ok242========================================" << endl;
  
  ws.import(errPdfSig);
  ws.import(errPdfBkg);
  
  cout << "====================ok25========================================" << endl;
  
  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("HIST_CtauErr_Tot_%s", (isPbp?"Pbp":"PP")), *ws.var("ctauErr"))->Clone(); 
  cout << "====================ok251========================================" << endl;
  RooDataHist* dataFULL = new RooDataHist(Form("hCtauErr_Tot_%s", (isPbp?"Pbp":"PP")), "", *ws.var("ctauErr"), hTot);
  ws.import(*dataFULL);
  hTot->Delete();
  
  cout << "====================ok252========================================" << endl;
  string pdfName = Form("pdfCTAUERR_Tot_%s", (isPbp?"Pbp":"PP"));
  RooHistPdf* pdf = new RooHistPdf( pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *dataFULL);
  ws.import(*pdf);
  cout << "====================ok253========================================" << endl;
  
  delete pdf;
  cout << "====================ok254========================================" << endl;
  delete dataFULL;
  //hTot->Clear();
  //hTot->Delete();
  
  cout << "====================ok255========================================" << endl;

  //ws.Print("v");
  
  return true;
  
};


bool makeMCCtauErrPdf(RooWorkspace& ws, string dsName, bool incPsi2S, bool isPbp) 
{

  RooDataSet *redMCSIG = NULL;
  if(!incPsi2S) 
    redMCSIG = (RooDataSet*)ws.data(dsName.c_str())->reduce("invMass > 2.9 && invMass < 3.3");
  
  if(incPsi2S)
    redMCSIG = (RooDataSet*)ws.data(dsName.c_str())->reduce("invMass > 3.55 && invMass < 3.85");
  
  // Extract the ctau error histograms for each mass range 
  RooDataHist *binMCCtErr = new RooDataHist("binMCCtErr","binMCCtErr",RooArgSet(*(ws.var("ctauErr"))), *ws.data(dsName.c_str()) );
  RooDataHist *binMCCtErrSIG = new RooDataHist("binMCCtErrSIG","MC ct error distribution for sig",RooArgSet(*(ws.var("ctauErr"))), *redMCSIG);
  ws.import(*binMCCtErr);
  ws.import(*binMCCtErrSIG);
  
  string Obj;
  if(!incPsi2S) {Obj = "Jpsi";}
  if(incPsi2S) {Obj = "Psi2S";}
  
  string pdfNameMCSIG = Form("pdfMCCTAUERR_%s_%s", Obj.c_str(), (isPbp?"Pbp":"PP"));
  RooHistPdf errPdfMCSig(pdfNameMCSIG.c_str(), "Error PDF MC signal", RooArgSet(*(ws.var("ctauErr"))) , *binMCCtErrSIG);
  
  ws.import(errPdfMCSig);
  
  string pdfName = Form("pdfMCCTAUERR_Tot_%s", (isPbp?"Pbp":"PP"));
  RooHistPdf* MCpdf = new RooHistPdf( pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *binMCCtErr);
  ws.import(*MCpdf);
  
  delete MCpdf;
  return true;
  
};




RooDataHist* subtractSidebands(RooWorkspace ws, RooDataHist* binSubtrSIG, RooDataHist* binSIG, RooDataHist* binSB, float scalefactor, string varName = "ctauErr") {
 
  if (binSIG->numEntries() != binSB->numEntries()) {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }
  RooDataHist* binScaleBKG = new RooDataHist("binScaleBKG","scaled SB",RooArgSet(*(ws.var(varName.c_str())))); 

  //// **** bin-by-bin scaling
  const RooArgSet* argSIG;
  const RooArgSet* argSB;
  for (Int_t i=0; i<binSIG->numEntries(); i++) {
    argSIG = binSIG->get(i);
    argSB = binSB->get(i);
    RooRealVar* thisVar = (RooRealVar*)argSIG->find(varName.c_str());
    ws.var(varName.c_str())->setVal(thisVar->getVal());
    //// *** set minimum as 0.1 to prevent empty PDF
    float wBkg = binSB->weight(*argSB,0,false);
    if (wBkg <= 0.1) wBkg = 0.1;
    binScaleBKG->add(RooArgSet(*(ws.var(varName.c_str()))),wBkg);
    float newWeight = binSIG->weight(*argSIG,0,false) - scalefactor*binSB->weight(*argSB,0,false);
    if (newWeight <= 0.1) newWeight = 0.1;
    binSubtrSIG->add(RooArgSet(*(ws.var(varName.c_str()))),newWeight);
  }
  return binScaleBKG;
}
