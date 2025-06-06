{

	// Use times new roman, precision 2 
	Int_t jpFont        = 42;  // Old LHCb style: 62;
	// Line thickness
	Double_t jpWidth    = 2.00; // Old LHCb style: 3.00;
	// Text size
	Double_t jpTSize    = 0.06;

	gROOT->SetStyle("Plain");
	TStyle *jpStyle= new TStyle("jpStyle","Jinping plots style");

	//jpStyle->SetFillColor(1);
	//jpStyle->SetFillStyle(1001);   // solid
	jpStyle->SetFrameFillColor(0);
	jpStyle->SetFrameBorderMode(0);
	jpStyle->SetPadBorderMode(0);
	jpStyle->SetPadColor(0);
	jpStyle->SetCanvasBorderMode(0);
	jpStyle->SetCanvasColor(0);
	jpStyle->SetStatColor(0);
	jpStyle->SetLegendBorderSize(0);



	jpStyle->SetPadTopMargin(0.135);
	jpStyle->SetPadRightMargin(0.08); // increase for colz plots
	jpStyle->SetPadBottomMargin(0.16);
	jpStyle->SetPadLeftMargin(0.15);

	// use large fonts
	jpStyle->SetTextFont(jpFont);
	jpStyle->SetTextSize(jpTSize);
	jpStyle->SetLabelFont(jpFont,"x");
	jpStyle->SetLabelFont(jpFont,"y");
	jpStyle->SetLabelFont(jpFont,"z");
	jpStyle->SetLabelSize(jpTSize,"x");
	jpStyle->SetLabelSize(jpTSize,"y");
	jpStyle->SetLabelSize(jpTSize,"z");
	jpStyle->SetTitleFont(jpFont);
	jpStyle->SetTitleFont(jpFont,"x");
	jpStyle->SetTitleFont(jpFont,"y");
	jpStyle->SetTitleFont(jpFont,"z");
	jpStyle->SetTitleSize(1.2*jpTSize,"x");
	jpStyle->SetTitleSize(1.2*jpTSize,"y");
	jpStyle->SetTitleSize(1.2*jpTSize,"z");
	jpStyle->SetTitleSize(1.5*jpTSize,"");

  //jpStyle->SetLegendTextSize(28);
  //jpStyle->SetLegendBorderSize(28);

	// histogram divisions: only 5 in x to avoid label overlaps
	jpStyle->SetNdivisions(505,"x");
	jpStyle->SetNdivisions(505,"y");
	jpStyle->SetNdivisions(505,"z");

	jpStyle->SetLabelOffset(0.010,"X");
	jpStyle->SetLabelOffset(0.010,"Y");
	jpStyle->SetLabelOffset(0.010,"Z");


	TGaxis::SetMaxDigits(3);

	// by default, do not display histogram decorations:
	jpStyle->SetOptStat(0);
	jpStyle->SetOptTitle(1);

	jpStyle->SetOptFit(1);

	//titles
	jpStyle->SetTitleOffset(0.95,"X");
	jpStyle->SetTitleOffset(0.95,"Y");
	jpStyle->SetTitleOffset(0.85,"Z");
  

  jpStyle->SetTitleFillColor(0);
  jpStyle->SetTitleStyle(0);
  jpStyle->SetTitleBorderSize(0);
  jpStyle->SetTitleFont(jpFont,"title");
  jpStyle->SetTitleX(0.1);
  jpStyle->SetTitleY(0.985); 
  //jpStyle->SetTitleW(1.0);
  //jpStyle->SetTitleH(0.10);



	gROOT->SetStyle("jpStyle");
	gROOT->ForceStyle();





}
