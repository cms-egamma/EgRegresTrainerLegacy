void style2014()
{

  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetTitleFont(42,"t");
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.457589);
  gStyle->SetStatY(0.312937);
  gStyle->SetStatW(0.29241/2+0.0185);
  gStyle->SetStatH(0.169580+0.05);
  gStyle->SetStatFontSize(0.0402098);
 
  gStyle->SetFitFormat("5.3g");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFontSize(0.040209);
  gStyle->SetStatFontSize(0.035209);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(900); //Width of canvas
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);

  gStyle->SetTitleX(0.1);

  gStyle->SetTitleXOffset(.9);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleYSize(0.055);
  gStyle->SetLabelSize(0.047,"XYZ");
}
