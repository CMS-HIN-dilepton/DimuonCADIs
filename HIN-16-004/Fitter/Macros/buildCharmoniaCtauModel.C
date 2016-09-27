#include "Utilities/initClasses.h"

void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbp);
void setCtauDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries);
bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp, bool usePromptCtauRes); 
bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp); 
bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp);
bool createCtauErrPDF(RooWorkspace& ws, string dsName, string pdfName, string histName, string cut="");

bool buildCharmoniaCtauModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, string dsName,
                             bool isPbp,                 // Determine if we are working with Pbp (True) or PP (False)
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             bool incPsi2S,               // Include Psi(2S) model
                             bool incPrompt,              // Include Prompt models
                             bool incNonPrompt,           // Include NonPrompt models
			     bool isErrPdfData,           // ctauErr pdf from data or MC
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{
  // If the initial parameters are empty, set defaul parameter values
  setCtauDefaultParameters(parIni, isPbp, numEntries);
  
  // Fix all psi2S parameters to jpsi
  if (incJpsi && incPsi2S) {
    fixCtauParPsi2StoJpsi(parIni, isPbp);
  }

  
  // C r e a t e   m o d e l 
     
  bool fitMass = false;
  if ( ws.pdf(Form("pdfMASS_Tot_%s", (isPbp?"Pbp":"PP"))) ) { fitMass = true; } 

  string pdfName     = "pdfCTAU";
  if (fitMass) { pdfName = "pdfCTAUMASS"; }
  bool isMC = (dsName.find("MC")!=std::string::npos);
  bool incCtauErrPDF = ( !isMC );

  cout << "====================ok1========================================" << endl;
  TH1::AddDirectory(kFALSE); 
  
  if (incJpsi) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "JpsiPR", model.CtauRes, parIni, isPbp, true)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "JpsiPR", model.Jpsi.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Jpsi Ctau Model failed" << endl; return false; }
    }
    
    cout << "====================ok4========================================" << endl;
    if (incNonPrompt) {
      if(!defineCtauResolModel(ws, "JpsiNoPR", model.CtauRes, parIni, isPbp, false)) { cout << "[ERROR] Defining the Non-Prompt Ctau Resolution Model failed" << endl; return false; }
      cout << "====================ok5========================================" << endl;
      if(!addSignalCtauModel(ws, "JpsiNoPR", model.Jpsi.Ctau.NonPrompt, parIni, isPbp)) { cout << "[ERROR] Adding NonPrompt Jpsi Ctau Model failed" << endl; return false; }
    }
    
    cout << "====================ok6========================================" << endl;
    
    if(!ws.pdf(Form("pdf%sCTAUERR_Jpsi_%s", (isErrPdfData ? "" : "MC"), (isPbp?"Pbp":"PP")))) {cout<<" not pdfMCCTAUERR_Jpsi_%s "<<endl;}

    //ws.Print("v");
    if (incPrompt) {
      RooProdPdf pdfJpsi(Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdf%sCTAUERR_Jpsi_%s", (isErrPdfData ? "" : "MC"), (isPbp?"Pbp":"PP"))),
			 Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau"))
				      )); 
      ws.import(pdfJpsi);
      //else {
      //ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")),
      //                Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))
      //                ));
      //}
    }
    
    cout << "====================ok7========================================" << endl;
    if (incPrompt && incNonPrompt) {
      
      if( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
                        ));
	}

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
                      ));
      
      if(!ws.var(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))) {cout<<" not sigmaMC_JpsiNoP || sigma2_CtauRes "<<endl;}
      //// **** convolution NPMCTrue with PRRes
	RooFormulaVar sigmaNPResN("sigmaNPResN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws.var(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP")))), *(ws.var("ctauErr")), *(ws.var(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))) )));  ws.import(sigmaNPResN);
      RooFormulaVar sigmaNPResW("sigmaNPResW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws.var(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP")))), *(ws.var("ctauErr")), *(ws.var(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))))) );  ws.import(sigmaNPResW);
      
      ws.factory(Form("GaussModel::GN_NPRes(ctau, %s,sigmaNPResN)", Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")) ));
      ws.factory(Form("GaussModel::GW_NPRes(ctau, %s,sigmaNPResW)", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")) ));
      ws.factory("AddModel::CtNPRes({GW_NPRes,GN_NPRes},{fracRes[0.01,0.001,0.999]})");
      //// **** final CtNP 
      float coefExpNPTrueVal = ws.var(Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))->getVal();
      ws.factory(Form("Decay::CtNPTot(ctau, %s,CtNPRes,RooDecay::SingleSided)",Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")) ));
      
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))) ){
	ws.factory(Form("PROD::%s(%s, CtNPTot)", Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
			Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
			));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		      Form("%s_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
		      ));
      
      
      if (!ws.pdf(Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) )) {
	cout<<"Total CtauMass_JpsiPR not found "<<Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))<<endl;
      }
      
      if(!ws.pdf(Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))) ) {
	cout<<"Total CtauMass_JpsiNoPR not found "<<Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))<<endl;
      }
      
      RooProdPdf pdfMassCtPR(Form("pdfMassCtPR_Jpsi_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
			     Conditional( *ws.pdf(Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ),
					  RooArgList(*ws.var("ctau"), *ws.var("invMass") ))
			     ); 
      ws.import(pdfMassCtPR);
      
      RooProdPdf pdfMassCtNoPR(Form("pdfMassCtNoPR_Jpsi_%s", (isPbp?"Pbp":"PP")), "",
			       *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
			       Conditional( *ws.pdf(Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ),
					    RooArgList(*ws.var("ctau"), *ws.var("invMass") ))
			       ); 
      ws.import(pdfMassCtNoPR);

      ws.factory( parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Jpsi_%s", (isPbp?"Pbp":"PP")),
		      Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")),
		      Form("pdfMassCtNoPR_Jpsi_%s", (isPbp?"Pbp":"PP")),
		      Form("pdfMassCtPR_Jpsi_%s", (isPbp?"Pbp":"PP"))
		      ));
      
      //ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Jpsi_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
      //              Form("%s_Jpsi_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
      //              Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
      //	      ));
    }
  }
  
  if (incPsi2S) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "Psi2SPR", model.CtauRes, parIni, isPbp, true)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "Psi2SPR", model.Psi2S.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Psi2S Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!defineCtauResolModel(ws, "Psi2SNoPR", model.CtauRes, parIni, isPbp, false)) { cout << "[ERROR] Defining the Non-Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "Psi2SNoPR", model.Psi2S.Ctau.NonPrompt, parIni, isPbp)) { cout << "[ERROR] Adding NonPrompt Psi2S Ctau Model failed" << endl; return false; }
    }
      
    cout << "====================ok6========================================" << endl;
    
    if (incPrompt) {
      RooProdPdf pdfPsi2S(Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdf%sCTAUERR_Psi2S_%s", (isErrPdfData ? "" : "MC"), (isPbp?"Pbp":"PP"))),
			 Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
			 ); 
      ws.import(pdfPsi2S);
      //else {
      //ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
      //                Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
      //                ));
      //}
    }
    
    cout << "====================ok7========================================" << endl;
    if (incPrompt && incNonPrompt) {
      
      if( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
                        ));
	}

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))
                      ));
      
      if(!ws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))) {cout<<" not sigmaMC_Psi2SNoP || sigma2_CtauRes "<<endl;}
      //// **** convolution NPMCTrue with PRRes
	RooFormulaVar sigmaNPResN("sigmaNPResN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws.var(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP")))), *(ws.var("ctauErr")), *(ws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))) )));  ws.import(sigmaNPResN);
      RooFormulaVar sigmaNPResW("sigmaNPResW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws.var(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP")))), *(ws.var("ctauErr")), *(ws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))))) );  ws.import(sigmaNPResW);
      
      ws.factory(Form("GaussModel::GN_NPRes(ctau, %s,sigmaNPResN)", Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")) ));
      ws.factory(Form("GaussModel::GW_NPRes(ctau, %s,sigmaNPResW)", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")) ));
      ws.factory("AddModel::CtNPRes({GW_NPRes,GN_NPRes},{fracRes[0.01,0.001,0.999]})");
      //// **** final CtNP 
      float coefExpNPTrueVal = ws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->getVal();
      ws.factory(Form("Decay::CtNPTot(ctau, %s,CtNPRes,RooDecay::SingleSided)",Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")) ));
      
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP"))) ){
	ws.factory(Form("PROD::%s(%s, CtNPTot)", Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
			Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
			));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		      Form("%s_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		      Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))
		      ));
      
      
      if (!ws.pdf(Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) )) {
	cout<<"Total CtauMass_Psi2SPR not found "<<Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))<<endl;
      }
      
      if(!ws.pdf(Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))) ) {
	cout<<"Total CtauMass_Psi2SNoPR not found "<<Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP"))<<endl;
      }
      
      RooProdPdf pdfMassCtPR(Form("pdfMassCtPR_Psi2S_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
			     Conditional( *ws.pdf(Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ),
					  RooArgList(*ws.var("ctau"), *ws.var("invMass") ))
			     ); 
      ws.import(pdfMassCtPR);
      
      RooProdPdf pdfMassCtNoPR(Form("pdfMassCtNoPR_Psi2S_%s", (isPbp?"Pbp":"PP")), "",
			       *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
			       Conditional( *ws.pdf(Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ),
					    RooArgList(*ws.var("ctau"), *ws.var("invMass") ))
			       ); 
      ws.import(pdfMassCtNoPR);

      ws.factory( parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Psi2S_%s", (isPbp?"Pbp":"PP")),
		      Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")),
		      Form("pdfMassCtNoPR_Psi2S_%s", (isPbp?"Pbp":"PP")),
		      Form("pdfMassCtPR_Psi2S_%s", (isPbp?"Pbp":"PP"))
		      ));
      
      //ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2S_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
      //              Form("%s_Psi2S_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
      //              Form("N_JPsi2S_%s", (isPbp?"Pbp":"PP"))
      //	      );)
    }
  }
  
  if (incBkg) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "BkgPR", model.CtauRes, parIni, isPbp, true)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addBackgroundCtauModel(ws, "BkgPR", model.Bkg.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Bkg Ctau Model failed" << endl; return false; }
    }

    if (incPrompt && incNonPrompt) {
      RooProdPdf pdfBkg(Form("pdfCTAU_Bkg_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
			Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
			);
      ws.import(pdfBkg);
      
      if (fitMass){
	if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))) ){
	  ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgPR_%s", (isPbp?"Pbp":"PP")),
			  Form("pdfCTAU_Bkg_%s", (isPbp?"Pbp":"PP")),
			  Form("pdfMASS_Bkg_%s",(isPbp?"Pbp":"PP"))
			  ));
	}
	//ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_BkgPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
	//		Form("%s_BkgPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
	//		Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))
	//		));
	
	cout << "====================okA4========================================" << endl;
	RooProdPdf pdfMassCtPRBkg(Form("pdfMassCtPR_Bkg_%s", (isPbp?"Pbp":"PP")), "",
				  *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
				  Conditional( *ws.pdf(Form("pdfCTAUMASS_BkgPR_%s", (isPbp?"Pbp":"PP")) ),
					       RooArgList(*ws.var("ctau"), *ws.var("invMass")) )
				  ); 
	ws.import(pdfMassCtPRBkg);


	cout << "====================okA5========================================" << endl;
	RooFormulaVar fracBkg("fracBkg","@0/(@0+@1)",RooArgList(*(ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))),*(ws.var(Form("N_%s_%s", (incJpsi ? "Jpsi" : "Psi2S") ,(isPbp?"Pbp":"PP"))  ))));
	ws.import(fracBkg);

	ws.factory( parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))].c_str() );
	if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))) ){
	  ws.factory(Form("PROD::%s(%s)", Form("pdfCTAUMASS_Bkg_%s", (isPbp?"Pbp":"PP")),
			  Form("pdfMassCtPR_Bkg_%s", (isPbp?"Pbp":"PP"))
			  ));
	}
	//ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
	//		Form("%s_Bkg_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
	//		Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))
	//		));
	cout << "====================okA6========================================" << endl;
      }
    }
  }
  
  ///ws.Print("v");
  
  // Total PDF
  RooAbsPdf *themodel = NULL;
  string tag = "";
  if (incPrompt  && !incNonPrompt) { tag = "PR"; }
  if (!incPrompt && incNonPrompt)  { tag = "NoPR";}
  if (incPrompt && incNonPrompt)  { tag = "";   }

  if (incJpsi && incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))) 
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))), 
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))),
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))) 
                                         )
                             );
  }
  if (incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))))
                                         )
                             );
  }
  if (incJpsi && !incPsi2S && incBkg) {
    /*themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))), 
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
					)
			     );*/
    ws.factory(Form("RSUM::%s(fracBkg*%s, %s*%s, %s)", Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		    Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")),
		    Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")), Form("pdfMassCtNoPR_Jpsi_%s", (isPbp?"Pbp":"PP")),
		    Form("pdfMassCtPR_Jpsi_%s", (isPbp?"Pbp":"PP")) 
		    ));
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ); 
    //ws->factory("RSUM::totPDF_PEE(fracBkg*MassCtBkg_PEE,Bfrac[0.25,0.0,1.]*MassCtNP_PEE,MassCtPR_PEE)");
  }
  
  if (!incJpsi && incPsi2S && incBkg) {
    /*themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))), *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
			     ));*/
    
    ws.factory(Form("RSUM::%s(fracBkg*%s, %s*%s, %s)", Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
		    Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")),
		    Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("pdfMassCtNoPR_Psi2S_%s", (isPbp?"Pbp":"PP")),
		    Form("pdfMassCtPR_Psi2S_%s", (isPbp?"Pbp":"PP")) 
		    ));
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")) ); 
    
  }
  if (incJpsi && !incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (!incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(),tag.c_str(),(isPbp?"Pbp":"PP")))
                                         ),  
                             RooArgList( 
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }

  //ws.import(*themodel);
  
  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("ctau"), *ws.var("invMass"), *ws.var("ctauErr"))) ;
  pdfName+=Form("_Tot_%s_parIni", (isPbp?"Pbp":"PP"));
  ws.saveSnapshot(pdfName.c_str(),*params,kTRUE) ;
  
  //ws.Print();
  return true;
};

void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbp)
{
  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi" << endl;
  parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP") ));
  parIni[Form("lambdaDSS_Psi2SNoPR_Pbp")] = Form("RooFormulaVar::%s('@0',{%s})", Form("lambdaDSS_Psi2SNoPR_Pbp"), Form("lambdaDSS_JpsiNoPR_Pbp"));
  //parIni[Form("b_Psi2S_Pbp")] = Form("RooFormulaVar::%s('@0',{%s})", Form("b_Psi2S_Pbp"), Form("b_Jpsi_Pbp"));
};


bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp, bool usePromptCtauRes) 
{ 

  if (ws.pdf(Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Ctau Resolution Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Ctau Resolution Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (CtauModel::SingleGaussianResolution):  
      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Single Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }

      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );  }
      if (!ws.var(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP")))) { ws.factory( parIni[Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }

      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, ctauErr)", Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))
                      ));
      
      cout << Form("[INFO] %s Single Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
      
    case (CtauModel::DoubleGaussianResolution):  
      
      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Double Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );      }
      
      // create the two PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, ctauErr)", Form("pdfCTAURES1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, ctauErr)", Form("pdfCTAURES2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))
                      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("pdfCTAURES1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("pdfCTAURES2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));

      //((RooGaussModel*)ws.pdf(Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))))->advertiseFlatScaleFactorIntegral(kTRUE);
      cout << Form("[INFO] %s Double Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
 
    default :
      cout<< "[ERROR] Selected Ctau Resolution Model has not been implemented"<< endl; return false;

    }
 
  return true;

};

bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{
  if (ws.pdf(Form("pdfCTAUTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Background Ctau Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Background Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   }

  switch(model) 
    {  
    case (CtauModel::TripleDecay): 
      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Background Triple Decay Ctau Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) && !ws.function(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      } 
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()  );
      ws.factory( parIni[Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()     ); 
      ws.factory( parIni[Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()     );

      // create the three PDFs
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      
      // combine the three PDFs
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAU1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
                              
      cout << Form("[INFO] %s Background Triple Decay Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    case (CtauModel::Delta): 
      
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) && !ws.function(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }

      // create the PDF
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));

      cout << Form("[INFO] %s Background Delta Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
      
    default :

      cout<< "[ERROR] Selected Background Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};


bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{
  if (ws.pdf(Form("pdfCTAUTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Signal Ctau Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Signal Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   } 

  switch(model) 
    {  
    case (CtauModel::SingleSidedDecay): 

      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Signal Single Sided Decay Ctau Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) && !ws.function(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
            
      // create the PDF
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
 
      cout << Form("[INFO] %s Signal Single Sided Decay Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    case (CtauModel::Delta):
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) && !ws.function(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }

      // create the three PDFs 
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
 		      
      cout << Form("[INFO] %s Signal Delta Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    default :

      cout<< "[ERROR] Selected Signal Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};


void setCtauDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;
  
  cout << "====================ok1a========================================" << endl;

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Bkg_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")), 0.25, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")), 0.2, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Bkg_%s", (isPbp?"Pbp":"PP")), 0.2, 0.0, 1.0);
  }

  cout << "====================ok2a========================================" << endl;
 // CTAU FIT PARAMETERS

 // Resolution Ctau Model
  if (parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("sigma1_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 2.0);
  }
  if (parIni.count(Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("sigma2_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 0.001, 60.0);
  }
  if (parIni.count(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.000001, 0.0000001, 1.0);
  }
  if (parIni.count(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.000001, 0.0000001, 1.0);
  }
  if (parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.8, 0.0, 1.0);
  }

  cout << "====================ok3a========================================" << endl;
  // Signal Ctau Model
  if (parIni.count(Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.01, 2.0);
  }
  if (parIni.count(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.01, 2.0);
  }

  cout << "====================ok4a========================================" << endl;
  // Background Ctau Model
  if (parIni.count(Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.5, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.45, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.30, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.06, 0.0001, 5.0);
  }
  if (parIni.count(Form("f_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.42, 0.05, 1.5);
  }
  if (parIni.count(Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 1.5);
  }
  if (parIni.count(Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.7, 0.001, 5.0);
  }
cout << "====================ok5a========================================" << endl;
};

