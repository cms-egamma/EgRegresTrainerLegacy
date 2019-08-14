{
  //#include <iomanip.h>

  gSystem->AddIncludePath("-I$PWD/include");
  gSystem->AddIncludePath("-I$ROOFITSYS/include");

  
  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("./libs/$SCRAM_ARCH/libUtility.so");
  gSystem->Load("./libs/$SCRAM_ARCH/libGBRLikelihood.so");
  gSystem->Load("./libs/$SCRAM_ARCH/libRegresTrainer.so");
  gSystem->Load("./libs/$SCRAM_ARCH/libResAnalysis.so");
    
  gROOT->SetStyle("Plain");
  
  gStyle->SetCanvasColor(0);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptFit(1111);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0); 
  gStyle->SetStripDecimals(false);
  
  gStyle->SetTitleXOffset(.9);
  gStyle->SetTitleXSize(0.047);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleYSize(0.047);
  
  gStyle->SetTitleX(0.04);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleW(0.88);
  gStyle->SetTitleH(0.06);
  
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

  gStyle->SetPadTopMargin(0.09);

  // gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadBottomMargin(0.105);

  gStyle->SetPadLeftMargin(0.13);

  gStyle->SetPadRightMargin(0.07);

  gStyle->SetStripDecimals(false);
  gROOT->ProcessLine(".x style2014.C");
  //  setTDRStyle();
//gStyle->SetTitleBorderSize(0);

  //gStyle->SetStatFont(62);
  //gStyle->SetStatFontSize(0.025);
  //gStyle->SetTitleFont(62,"t");
  //gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleFont(62, "XYZ");
  //  gStyle->SetLabelFont(62, "XYZ");
  // gStyle->SetLegendFont(62);


}
